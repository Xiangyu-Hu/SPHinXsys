/**
 * @file test_3d_em_aphi_matrix_free_harmonic_gradient_variable_sigma.cpp
 * @brief Manufactured verification for matrix-free harmonic-weighted gradient operator.
 *
 * Validates:
 *   1) harmonic-gradient(constant field, variable sigma) = 0
 *   2) for a linear field, the operator matches sigma * grad(phi) away from the
 *      material interface and away from the outer boundary
 *
 * Optional SYCL path:
 *   EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_weighted_gradient.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <iomanip>
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

class HarmonicGradientBoxShape : public ComplexShape
{
  public:
    HarmonicGradientBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Vec3c makeVec3c(const Complex &x, const Complex &y, const Complex &z)
{
    Vec3c result = Vec3c::Zero();
    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
}

Real complexVecErrorNorm(const Vec3c &lhs, const Vec3c &rhs)
{
    Real squared_sum = 0.0;
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        squared_sum += std::norm(lhs[axis] - rhs[axis]);
    }
    return std::sqrt(squared_sum);
}

Real complexVecAbsNorm(const Vec3c &value)
{
    Real squared_sum = 0.0;
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        squared_sum += std::norm(value[axis]);
    }
    return std::sqrt(squared_sum);
}

Real distanceToBoundary(const Vecd &position, Real length, Real height, Real width)
{
    return std::min({position[0], length - position[0], position[1], height - position[1], position[2], width - position[2]});
}

Complex constantFieldValue()
{
    return Complex(-0.15, 0.45);
}

Complex linearFieldValue(const Vecd &position, const Vecd &center)
{
    const Complex gx(0.55, 0.12);
    const Complex gy(-0.35, 0.18);
    const Complex gz(0.22, -0.14);
    return Complex(0.25, -0.05) + gx * (position[0] - center[0]) + gy * (position[1] - center[1]) +
           gz * (position[2] - center[2]);
}

Vec3c linearGradientExact()
{
    return makeVec3c(Complex(0.55, 0.12), Complex(-0.35, 0.18), Complex(0.22, -0.14));
}

struct Summary
{
    size_t total_particles = 0;
    size_t region_core_particles = 0;
    Real constant_max_abs = 0.0;
    Real linear_region_core_mean_error = 0.0;
    Real linear_region_core_max_error = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_DP", 0.10);
    const Real body_length = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_LENGTH", 1.0);
    const Real body_height = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_HEIGHT", 1.0);
    const Real body_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_WIDTH", 1.0);
    const Real boundary_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Real core_shell = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_CORE_SHELL", 1.5 * dp_0);
    const Real interface_shell = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_INTERFACE_SHELL", 1.5 * dp_0);
    const Real sigma_left = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_SIGMA_LEFT", 1.0);
    const Real sigma_right = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_SIGMA_RIGHT", 100.0);
    const Real validation_constant_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_MAX_CONSTANT_ABS", 1.0e-8);
    const Real validation_linear_core_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_HARMONIC_MAX_LINEAR_REGION_CORE", 2.0e-1);
    const bool use_sycl_harmonic_gradient = getEnvBoolLocal("EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT", false);
    const bool use_sph_native_harmonic_gradient = getEnvBoolLocal("EM_APHI_USE_SPH_NATIVE_HARMONIC_GRADIENT", false);

    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Real interface_x = body_center[0];
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system,
                   makeShared<HarmonicGradientBoxShape>("MatrixFreeHarmonicGradientBody", body_center, body_halfsize));
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

#if SPHINXSYS_USE_SYCL
    setMatrixFreeHarmonicGradientUseSycl(use_sycl_harmonic_gradient);
    if (use_sycl_harmonic_gradient)
    {
        prepareMatrixFreeAPhiSyclDeviceResources(graph, number_of_particles);
    }
#else
    (void)use_sycl_harmonic_gradient;
#endif

    StdVec<Real> sigma(number_of_particles, sigma_left);
    StdVec<Complex> constant_field(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> linear_field(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        sigma[i] = positions[i][0] < interface_x ? sigma_left : sigma_right;
        constant_field[i] = constantFieldValue();
        linear_field[i] = linearFieldValue(positions[i], body_center);
    }

    std::unique_ptr<SPHComplexHarmonicWeightedGradientOperator> sph_native_harmonic_gradient_operator;
    if (use_sph_native_harmonic_gradient)
    {
        sph_native_harmonic_gradient_operator =
            std::make_unique<SPHComplexHarmonicWeightedGradientOperator>(body, inner_ck);
    }

    const auto compute_harmonic_gradient = [&](const StdVec<Complex> &scalar_field) -> StdVec<Vec3c>
    {
        if (use_sph_native_harmonic_gradient)
        {
            return sph_native_harmonic_gradient_operator->computeFromField(scalar_field, sigma);
        }
#if SPHINXSYS_USE_SYCL
        if (use_sycl_harmonic_gradient)
        {
            setMatrixFreeHarmonicGradientUseSycl(true);
            const StdVec<Vec3c> result = applyMatrixFreeHarmonicWeightedGradient(graph, scalar_field, sigma);
            setMatrixFreeHarmonicGradientUseSycl(false);
            return result;
        }
#endif
        return applyMatrixFreeHarmonicWeightedGradient(graph, scalar_field, sigma);
    };

    const StdVec<Vec3c> harmonic_gradient_constant = compute_harmonic_gradient(constant_field);
    const StdVec<Vec3c> harmonic_gradient_linear = compute_harmonic_gradient(linear_field);

    Summary summary;
    const Vec3c linear_gradient_exact = linearGradientExact();
    Real linear_region_core_error_sum = 0.0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const bool away_from_boundary = distanceToBoundary(positions[i], body_length, body_height, body_width) > core_shell;
        const bool away_from_interface = std::abs(positions[i][0] - interface_x) > interface_shell;
        const bool region_core = away_from_boundary && away_from_interface;
        const Real sigma_i = sigma[i];
        const Real constant_abs = complexVecAbsNorm(harmonic_gradient_constant[i]);
        const Vec3c recovered_gradient = makeVec3c(harmonic_gradient_linear[i][0] / sigma_i,
                                                   harmonic_gradient_linear[i][1] / sigma_i,
                                                   harmonic_gradient_linear[i][2] / sigma_i);
        const Real linear_error = complexVecErrorNorm(recovered_gradient, linear_gradient_exact);

        summary.total_particles++;
        summary.constant_max_abs = SMAX(summary.constant_max_abs, constant_abs);
        if (region_core)
        {
            summary.region_core_particles++;
            linear_region_core_error_sum += linear_error;
            summary.linear_region_core_max_error = SMAX(summary.linear_region_core_max_error, linear_error);
        }
    }

    summary.linear_region_core_mean_error =
        linear_region_core_error_sum / (static_cast<Real>(summary.region_core_particles) + TinyReal);

    bool validation_pass = true;
    validation_pass = validation_pass && summary.region_core_particles > 0;
    validation_pass = validation_pass && summary.constant_max_abs <= validation_constant_max;
    validation_pass = validation_pass && summary.linear_region_core_max_error <= validation_linear_core_max;

    std::cout << std::setprecision(12)
              << "matrix_free_harmonic_gradient_variable_sigma"
              << " points=" << number_of_particles
              << " dp=" << dp_0
              << " sycl_harmonic_gradient=" << (use_sycl_harmonic_gradient ? 1 : 0)
              << " sph_native_harmonic_gradient=" << (use_sph_native_harmonic_gradient ? 1 : 0)
              << " harmonic_gradient_backend="
              << (use_sph_native_harmonic_gradient ? "sph_native_ck" : (use_sycl_harmonic_gradient ? "graph_sycl" : "graph_cpu"))
              << " sigma_left=" << sigma_left
              << " sigma_right=" << sigma_right
              << " region_core_particles=" << summary.region_core_particles
              << " constant_max_abs=" << summary.constant_max_abs
              << " linear_region_core_mean_error=" << summary.linear_region_core_mean_error
              << " linear_region_core_max_error=" << summary.linear_region_core_max_error
              << " validation_max_constant_abs=" << validation_constant_max
              << " validation_max_linear_region_core=" << validation_linear_core_max
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << std::endl;

    return validation_pass ? 0 : 1;
}
