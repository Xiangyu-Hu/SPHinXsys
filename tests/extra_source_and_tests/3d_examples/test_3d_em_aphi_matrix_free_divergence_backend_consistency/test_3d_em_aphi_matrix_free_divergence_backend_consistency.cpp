/**
 * @file test_3d_em_aphi_matrix_free_divergence_backend_consistency.cpp
 * @brief Backend consistency for matrix-free vector divergence and sigma-grad-phi assembly.
 *
 * Validates:
 *   1) div(constant vector) = 0 (standard and harmonic-weighted)
 *   2) sigma * grad(phi) assembly via harmonic gradient
 *   3) SPH-native CK backends match graph CPU reference within tolerance
 *
 * Optional:
 *   EM_APHI_USE_SPH_NATIVE_DIVERGENCE=1
 *   EM_APHI_SPH_NATIVE_DIVERGENCE_MAX_GRAPH_DIFF (default 1e-2)
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_assembly.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;
using namespace SPH::electromagnetics::sph;

namespace
{

Real getEnvRealLocal(const std::string &name, Real default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    char *end_ptr = nullptr;
    const Real parsed = static_cast<Real>(std::strtod(value, &end_ptr));
    return end_ptr == value || !std::isfinite(parsed) ? default_value : parsed;
}

bool getEnvBoolLocal(const std::string &name, bool default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    return std::atoi(value) != 0;
}

class DivergenceBoxShape : public ComplexShape
{
  public:
    DivergenceBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Real complexAbs(const Complex &value)
{
    return std::abs(value);
}

Real complexAbsError(const Complex &lhs, const Complex &rhs)
{
    return std::abs(lhs - rhs);
}

Complex constantComponentValue(int component)
{
    const Real offsets[3] = {-0.2, 0.35, -0.15};
    const Real imags[3] = {0.4, -0.25, 0.5};
    return Complex(offsets[component], imags[component]);
}

Complex linearComponentValue(int component, const Vecd &position, const Vecd &center)
{
    const Complex slopes[3] = {Complex(0.6, 0.1), Complex(-0.2, 0.3), Complex(0.15, -0.12)};
    const Complex offsets[3] = {Complex(0.1, -0.05), Complex(-0.08, 0.2), Complex(0.12, 0.04)};
    return offsets[component] + slopes[component] * (position[component] - center[component]);
}

struct Summary
{
    Real constant_standard_max_abs = 0.0;
    Real constant_harmonic_max_abs = 0.0;
    Real standard_max_graph_diff = 0.0;
    Real harmonic_max_graph_diff = 0.0;
    Real sigma_grad_constant_max_abs = 0.0;
    Real sigma_grad_max_graph_diff = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_DP", 0.10);
    const Real body_length = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_LENGTH", 1.0);
    const Real body_height = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_HEIGHT", 1.0);
    const Real body_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_WIDTH", 1.0);
    const Real boundary_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Real sigma_left = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_SIGMA_LEFT", 1.0);
    const Real sigma_right = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_SIGMA_RIGHT", 100.0);
    const Real validation_constant_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_DIVERGENCE_MAX_CONSTANT_ABS", 1.0e-8);
    const Real validation_max_graph_diff = getEnvRealLocal("EM_APHI_SPH_NATIVE_DIVERGENCE_MAX_GRAPH_DIFF", 2.0e-2);
    const bool use_sph_native_divergence = getEnvBoolLocal("EM_APHI_USE_SPH_NATIVE_DIVERGENCE", false);

    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Real interface_x = body_center[0];
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system, makeShared<DivergenceBoxShape>("MatrixFreeDivergenceBody", body_center, body_halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);
    Inner<> inner_ck(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    parameters.max_change_rate = MaxReal;
    parameters.interface_contrast_threshold = 10.0;
    MatrixFreeAPhiDiscreteView discrete_view{number_of_particles,
                                             body.getSPHAdaptation().ReferenceSmoothingLength(),
                                             positions,
                                             volumetric_measure,
                                             nullptr,
                                             &inner_relation.inner_configuration_,
                                             nullptr};
    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);

    StdVec<Real> sigma(number_of_particles, sigma_left);
    StdVec<Complex> constant_x(number_of_particles);
    StdVec<Complex> constant_y(number_of_particles);
    StdVec<Complex> constant_z(number_of_particles);
    StdVec<Complex> linear_x(number_of_particles);
    StdVec<Complex> linear_y(number_of_particles);
    StdVec<Complex> linear_z(number_of_particles);
    StdVec<Complex> constant_phi(number_of_particles, constantComponentValue(0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        sigma[i] = positions[i][0] < interface_x ? sigma_left : sigma_right;
        constant_x[i] = constantComponentValue(0);
        constant_y[i] = constantComponentValue(1);
        constant_z[i] = constantComponentValue(2);
        linear_x[i] = linearComponentValue(0, positions[i], body_center);
        linear_y[i] = linearComponentValue(1, positions[i], body_center);
        linear_z[i] = linearComponentValue(2, positions[i], body_center);
    }

    const StdVec<Complex> graph_div_constant =
        computeGradientDivergenceOfVectorField(graph, constant_x, constant_y, constant_z);
    const StdVec<Complex> graph_div_harmonic_constant =
        computeHarmonicDivergenceOfVectorField(graph, constant_x, constant_y, constant_z, sigma);
    const StdVec<Complex> graph_div_linear =
        computeGradientDivergenceOfVectorField(graph, linear_x, linear_y, linear_z);
    const StdVec<Complex> graph_div_harmonic_linear =
        computeHarmonicDivergenceOfVectorField(graph, linear_x, linear_y, linear_z, sigma);
    const StdVec<Vec3c> graph_sigma_grad_constant =
        applyMatrixFreeHarmonicWeightedGradient(graph, constant_phi, sigma);
    const StdVec<Vec3c> graph_sigma_grad_linear =
        applyMatrixFreeHarmonicWeightedGradient(graph, linear_x, sigma);

    std::unique_ptr<SPHAPhiCoupledOperators> sph_native_operators;
    if (use_sph_native_divergence)
    {
        sph_native_operators = std::make_unique<SPHAPhiCoupledOperators>(body, inner_ck);
    }

    const auto compute_standard_divergence = [&](const StdVec<Complex> &fx, const StdVec<Complex> &fy,
                                                 const StdVec<Complex> &fz) -> StdVec<Complex>
    {
        if (use_sph_native_divergence)
        {
            return sph_native_operators->computeGradientDivergenceOfVectorField(fx, fy, fz);
        }
        return computeGradientDivergenceOfVectorField(graph, fx, fy, fz);
    };

    const auto compute_harmonic_divergence = [&](const StdVec<Complex> &fx, const StdVec<Complex> &fy,
                                                 const StdVec<Complex> &fz) -> StdVec<Complex>
    {
        if (use_sph_native_divergence)
        {
            return sph_native_operators->computeHarmonicDivergenceOfVectorField(fx, fy, fz, sigma);
        }
        return computeHarmonicDivergenceOfVectorField(graph, fx, fy, fz, sigma);
    };

    const auto compute_sigma_grad_phi = [&](const StdVec<Complex> &phi) -> StdVec<Vec3c>
    {
        if (use_sph_native_divergence)
        {
            return sph_native_operators->computeSigmaGradPhi(phi, sigma);
        }
        return applyMatrixFreeHarmonicWeightedGradient(graph, phi, sigma);
    };

    const StdVec<Complex> backend_div_constant = compute_standard_divergence(constant_x, constant_y, constant_z);
    const StdVec<Complex> backend_div_harmonic_constant = compute_harmonic_divergence(constant_x, constant_y, constant_z);
    const StdVec<Complex> backend_div_linear = compute_standard_divergence(linear_x, linear_y, linear_z);
    const StdVec<Complex> backend_div_harmonic_linear = compute_harmonic_divergence(linear_x, linear_y, linear_z);
    const StdVec<Vec3c> backend_sigma_grad_constant = compute_sigma_grad_phi(constant_phi);
    const StdVec<Vec3c> backend_sigma_grad_linear = compute_sigma_grad_phi(linear_x);

    Summary summary;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        summary.constant_standard_max_abs = SMAX(summary.constant_standard_max_abs, complexAbs(backend_div_constant[i]));
        summary.constant_harmonic_max_abs =
            SMAX(summary.constant_harmonic_max_abs, complexAbs(backend_div_harmonic_constant[i]));
        summary.standard_max_graph_diff =
            SMAX(summary.standard_max_graph_diff, complexAbsError(backend_div_linear[i], graph_div_linear[i]));
        summary.harmonic_max_graph_diff = SMAX(summary.harmonic_max_graph_diff,
                                                 complexAbsError(backend_div_harmonic_linear[i], graph_div_harmonic_linear[i]));
        Real sigma_grad_constant_abs = 0.0;
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            sigma_grad_constant_abs += std::norm(backend_sigma_grad_constant[i][axis]);
            summary.sigma_grad_max_graph_diff =
                SMAX(summary.sigma_grad_max_graph_diff, std::abs(backend_sigma_grad_linear[i][axis] -
                                                                graph_sigma_grad_linear[i][axis]));
        }
        summary.sigma_grad_constant_max_abs = SMAX(summary.sigma_grad_constant_max_abs, std::sqrt(sigma_grad_constant_abs));
    }

    if (!use_sph_native_divergence)
    {
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            summary.standard_max_graph_diff =
                SMAX(summary.standard_max_graph_diff, complexAbsError(graph_div_constant[i], backend_div_constant[i]));
            summary.harmonic_max_graph_diff = SMAX(summary.harmonic_max_graph_diff,
                                                   complexAbsError(graph_div_harmonic_constant[i], backend_div_harmonic_constant[i]));
        }
    }

    bool validation_pass = true;
    validation_pass = validation_pass && summary.constant_standard_max_abs <= validation_constant_max;
    validation_pass = validation_pass && summary.constant_harmonic_max_abs <= validation_constant_max;
    validation_pass = validation_pass && summary.sigma_grad_constant_max_abs <= validation_constant_max;
    if (use_sph_native_divergence)
    {
        validation_pass = validation_pass && summary.standard_max_graph_diff <= validation_max_graph_diff;
        validation_pass = validation_pass && summary.harmonic_max_graph_diff <= validation_max_graph_diff;
        validation_pass = validation_pass && summary.sigma_grad_max_graph_diff <= validation_max_graph_diff;
    }

    std::cout << std::setprecision(12)
              << "matrix_free_divergence_backend_consistency"
              << " points=" << number_of_particles
              << " dp=" << dp_0
              << " sph_native_divergence=" << (use_sph_native_divergence ? 1 : 0)
              << " divergence_backend=" << (use_sph_native_divergence ? "sph_native_ck" : "graph_cpu")
              << " constant_standard_max_abs=" << summary.constant_standard_max_abs
              << " constant_harmonic_max_abs=" << summary.constant_harmonic_max_abs
              << " standard_max_graph_diff=" << summary.standard_max_graph_diff
              << " harmonic_max_graph_diff=" << summary.harmonic_max_graph_diff
              << " sigma_grad_constant_max_abs=" << summary.sigma_grad_constant_max_abs
              << " sigma_grad_max_graph_diff=" << summary.sigma_grad_max_graph_diff
              << " validation_max_constant_abs=" << validation_constant_max
              << " validation_max_graph_diff=" << validation_max_graph_diff
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << std::endl;

    return validation_pass ? 0 : 1;
}
