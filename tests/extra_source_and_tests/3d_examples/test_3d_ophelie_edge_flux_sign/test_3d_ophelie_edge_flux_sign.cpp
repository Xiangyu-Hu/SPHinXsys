/**
 * @file test_3d_ophelie_edge_flux_sign.cpp
 * @brief Edge-flux sign MMS for imag and real scalar chains (edge_drop -> 0).
 */
#include "electromagnetic_ophelie.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline bool checkEdgeDropSignCase(const OphelieEdgeFluxEdgeDropMetrics &drop_metrics, Real phi_scale,
                                  const char *chain_name)
{
    const Real edge_drop_linf_rel = drop_metrics.edge_drop_linf / (phi_scale + TinyReal);
    const Real edge_drop_l2_rel = drop_metrics.edge_drop_l2 / (phi_scale + TinyReal);
    const bool passed = edge_drop_linf_rel < Real(1e-4) && edge_drop_l2_rel < Real(1e-4);
    std::cout << "test_3d_ophelie_edge_flux_sign " << chain_name << " edge_drop_linf=" << drop_metrics.edge_drop_linf
              << " edge_drop_l2=" << drop_metrics.edge_drop_l2 << " edge_drop_linf_rel=" << edge_drop_linf_rel
              << " edge_drop_l2_rel=" << edge_drop_l2_rel << " edge_drop_mean_abs=" << drop_metrics.edge_drop_mean_abs
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed;
}

} // namespace

int main(int, char *[])
{
    OphelieParameters params;
    params.frequency_ = 50.0e3;
    params.sigma_glass_ = 16.0;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    OphelieTestCliOptions cli_options;
    finalizeOphelieCurrentFormConfiguration(params, cli_options);

    const Vecd a_uniform(1.5, 0.0, 0.0);
    const Real omega = params.omega();

    const Real dp = 0.04;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names,
                                                                               params.sigma_glass_);
    assign_sigma.exec();

    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    Real *phi_real = particles.getVariableDataByName<Real>(glass_names.phi_real);
    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = a_uniform;
        a_coil_imag[i] = Vecd::Zero();
        phi_imag[i] = -omega * a_uniform.dot(pos[i]);
        phi_real[i] = Real(0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);

    const OphelieEdgeFluxQAntisymMetrics q_antisym = evaluateOphelieEdgeFluxQAntisymmetryForComponent<MainExecutionPolicy>(
        glass_body, *glass_inner, glass_names, makeOphelieEdgeFluxImagComponent(glass_names, params), params);
    const OphelieEdgeFluxResidualMetrics imag_residual_metrics =
        evaluateOphelieEdgeFluxResidual<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    const OphelieEdgeFluxEdgeDropMetrics imag_drop_metrics =
        evaluateOphelieEdgeFluxImagEdgeDropMetrics<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    const Real phi_scale = std::abs(omega * a_uniform.norm() * halfsize[0]);
    const Real edge_res_rel = imag_residual_metrics.edge_res_l2 / (phi_scale + TinyReal);
    const bool imag_passed = n > 0 && edge_res_rel < 0.05 &&
                             checkEdgeDropSignCase(imag_drop_metrics, phi_scale, "imag_chain") &&
                             q_antisym.q_nonfinite_count == 0 && q_antisym.q_antisym_rel_l2 < Real(1e-5);

    for (size_t i = 0; i < n; ++i)
    {
        a_coil_real[i] = Vecd::Zero();
        a_coil_imag[i] = a_uniform;
        phi_imag[i] = Real(0);
        phi_real[i] = omega * a_uniform.dot(pos[i]);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_coil_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);
    syncVariableToDevice<Real>(particles, glass_names.phi_real);

    const OphelieEdgeFluxEdgeDropMetrics real_drop_metrics =
        evaluateOphelieEdgeFluxRealEdgeDropMetrics<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    const OphelieEdgeFluxQAntisymMetrics real_q_antisym = evaluateOphelieEdgeFluxQAntisymmetryForComponent<
        MainExecutionPolicy>(glass_body, *glass_inner, glass_names,
                             makeOphelieEdgeFluxRealComponent(glass_names, params), params);
    const bool real_passed =
        checkEdgeDropSignCase(real_drop_metrics, phi_scale, "real_chain") && real_q_antisym.q_nonfinite_count == 0 &&
        real_q_antisym.q_antisym_rel_l2 < Real(1e-5);
    const bool passed = imag_passed && real_passed;

    std::cout << "test_3d_ophelie_edge_flux_sign n=" << n << " edge_res_l2=" << imag_residual_metrics.edge_res_l2
              << " edge_res_rel=" << edge_res_rel               << " q_antisym_rel_l2=" << q_antisym.q_antisym_rel_l2
              << " real_q_antisym_rel_l2=" << real_q_antisym.q_antisym_rel_l2
              << " q_nonfinite_count=" << q_antisym.q_nonfinite_count << " phi_scale=" << phi_scale
              << " edge_res_integral=" << imag_residual_metrics.edge_res_volume_integral
              << " imag_chain_passed=" << (imag_passed ? 1 : 0) << " real_chain_passed=" << (real_passed ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
