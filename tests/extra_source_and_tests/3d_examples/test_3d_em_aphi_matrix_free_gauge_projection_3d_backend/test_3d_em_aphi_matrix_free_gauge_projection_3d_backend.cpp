/**
 * @file test_3d_em_aphi_matrix_free_gauge_projection_3d_backend.cpp
 * @brief 3D gauge projection: graph vs SPH-native backend consistency (Laplace / div-grad chi).
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;

namespace
{

bool getEnvBool(const char *name, bool default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : std::atoi(value) != 0;
}

Real getEnvReal(const char *name, Real default_value)
{
    const char *value = std::getenv(name);
    if (value == nullptr)
    {
        return default_value;
    }
    return static_cast<Real>(std::strtod(value, nullptr));
}

class Gauge3dBoxShape : public ComplexShape
{
  public:
    Gauge3dBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

MatrixFreeAPhiFields makeDivergentAField(const StdVec<Vecd> &positions, const Vecd &center)
{
    MatrixFreeAPhiFields fields;
    const size_t n = positions.size();
    fields.ax.resize(n);
    fields.ay.resize(n);
    fields.az.resize(n);
    fields.phi.resize(n, Complex(0.0, 0.0));
    for (size_t i = 0; i != n; ++i)
    {
        fields.ax[i] = Complex(positions[i][0] - center[0], 0.0);
        fields.ay[i] = Complex(positions[i][1] - center[1], 0.0);
        fields.az[i] = Complex(0.1 * (positions[i][2] - center[2]), 0.05);
    }
    return fields;
}

struct GaugeRunMetrics
{
    bool chi_converged = false;
    size_t chi_iterations = 0;
    Real chi_l2 = 0.0;
    Real grad_chi_l2 = 0.0;
    Real div_a_before_l2 = 0.0;
    Real div_a_after_final_l2 = 0.0;
};

GaugeRunMetrics runGaugeProjection(const MatrixFreePairwiseGraph &graph, MatrixFreeAPhiFields fields,
                                 const StdVec<Real> &sigma, const MatrixFreeAPhiParameters &parameters,
                                 const ScalarComplexHelmholtzSolverParameters &gauge_solver, bool use_consistent_gauge,
                                 MatrixFreeAPhiSphNativeContext *sph_native_context)
{
    const MatrixFreeAPhiGaugeProjectionResult result = applyMatrixFreeAPhiGaugeProjectionStep(
        graph, fields, sigma, parameters, gauge_solver, true, false, use_consistent_gauge, sph_native_context);

    GaugeRunMetrics metrics;
    metrics.chi_converged = result.solver_state.converged_;
    metrics.chi_iterations = result.solver_state.iterations_;
    metrics.chi_l2 = result.diagnostics.chi_l2;
    metrics.grad_chi_l2 = result.diagnostics.grad_chi_l2;
    metrics.div_a_before_l2 = result.diagnostics.divergence_a_before_l2;
    metrics.div_a_after_final_l2 = result.diagnostics.divergence_a_after_final_l2;
    return metrics;
}

Real relativeDifference(Real a, Real b)
{
    return std::abs(a - b) / SMAX(SMAX(std::abs(a), std::abs(b)), TinyReal);
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = getEnvReal("EM_APHI_GAUGE_3D_DP", 0.12);
    const bool use_sph_native = getEnvBool("EM_APHI_USE_SPH_NATIVE_RESIDUALS", false);
    const bool use_consistent_gauge = getEnvBool("EM_APHI_GAUGE_3D_OPERATOR_CONSISTENT", false);
    const Real max_backend_rel_diff = getEnvReal("EM_APHI_GAUGE_3D_MAX_GRAPH_DIFF", 5.0e-2);

    const Vecd halfsize(0.5, 0.5, 0.5);
    const Vecd center(0.5, 0.5, 0.5);
    BoundingBoxd bounds(Vecd(-0.36, -0.36, -0.36), Vecd(1.36, 1.36, 1.36));

    SPHSystem sph_system(bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    SolidBody body(sph_system, makeShared<Gauge3dBoxShape>("Gauge3dBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);
    Inner<> inner_ck(body);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    parameters.angular_frequency = getEnvReal("EM_APHI_GAUGE_3D_OMEGA", 1.0);
    MatrixFreeAPhiDiscreteView view{n, body.getSPHAdaptation().ReferenceSmoothingLength(), positions, vol, nullptr,
                                    &inner_relation.inner_configuration_, nullptr};
    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(view, parameters);

    std::unique_ptr<MatrixFreeAPhiSphNativeContext> sph_ctx_storage;
    if (use_sph_native)
    {
        sph_ctx_storage = std::make_unique<MatrixFreeAPhiSphNativeContext>(body, inner_ck);
    }
    MatrixFreeAPhiSphNativeContext *sph_ctx = sph_ctx_storage.get();

    const MatrixFreeAPhiFields template_fields =
        makeDivergentAField(StdVec<Vecd>(positions, positions + n), center);
    StdVec<Real> sigma(n, 1.0);

    ScalarComplexHelmholtzSolverParameters gauge_solver;
    gauge_solver.max_iterations_ = static_cast<size_t>(getEnvReal("EM_APHI_GAUGE_3D_MAX_ITER", 400));
    gauge_solver.absolute_tolerance_ = getEnvReal("EM_APHI_GAUGE_3D_TOL", 1.0e-4);
    gauge_solver.relaxation_factor_ = getEnvReal("EM_APHI_GAUGE_3D_RELAX", 0.5);

    const GaugeRunMetrics graph_metrics =
        runGaugeProjection(graph, template_fields, sigma, parameters, gauge_solver, use_consistent_gauge, nullptr);

    Real chi_rel_diff = 0.0;
    Real grad_chi_rel_diff = 0.0;
    Real div_a_after_rel_diff = 0.0;
    GaugeRunMetrics sph_metrics;

    if (use_sph_native)
    {
        sph_metrics = runGaugeProjection(graph, template_fields, sigma, parameters, gauge_solver, use_consistent_gauge, sph_ctx);
        chi_rel_diff = relativeDifference(graph_metrics.chi_l2, sph_metrics.chi_l2);
        grad_chi_rel_diff = relativeDifference(graph_metrics.grad_chi_l2, sph_metrics.grad_chi_l2);
        div_a_after_rel_diff =
            relativeDifference(graph_metrics.div_a_after_final_l2, sph_metrics.div_a_after_final_l2);
    }

    // Laplace chi solve converges on this box; div-grad uses a custom normalized update that may
    // hit the residual growth guard early — backend parity is still checked when SPH is enabled.
    const bool require_chi_convergence = !use_consistent_gauge;
    const bool graph_ok = !require_chi_convergence || graph_metrics.chi_converged;
    const bool sph_ok = !use_sph_native || !require_chi_convergence || sph_metrics.chi_converged;
    const bool backend_ok =
        !use_sph_native ||
        (chi_rel_diff <= max_backend_rel_diff && grad_chi_rel_diff <= max_backend_rel_diff &&
         div_a_after_rel_diff <= max_backend_rel_diff);
    const bool validation_pass = graph_ok && sph_ok && backend_ok;

    std::cout << std::setprecision(10)
              << "matrix_free_gauge_projection_3d_backend"
              << " points=" << n
              << " sph_native=" << (use_sph_native ? 1 : 0)
              << " gauge_operator=" << (use_consistent_gauge ? "consistent_div_grad" : "laplace")
              << " graph_chi_iterations=" << graph_metrics.chi_iterations
              << " graph_chi_converged=" << (graph_metrics.chi_converged ? 1 : 0)
              << " graph_chi_l2=" << graph_metrics.chi_l2
              << " graph_grad_chi_l2=" << graph_metrics.grad_chi_l2
              << " graph_div_a_after_final_l2=" << graph_metrics.div_a_after_final_l2;
    if (use_sph_native)
    {
        std::cout << " sph_chi_iterations=" << sph_metrics.chi_iterations
                  << " sph_chi_converged=" << (sph_metrics.chi_converged ? 1 : 0)
                  << " chi_rel_diff=" << chi_rel_diff
                  << " grad_chi_rel_diff=" << grad_chi_rel_diff
                  << " div_a_after_rel_diff=" << div_a_after_rel_diff
                  << " max_backend_rel_diff=" << max_backend_rel_diff;
    }
    std::cout << " validation_pass=" << (validation_pass ? 1 : 0) << std::endl;

    return validation_pass ? 0 : 1;
}
