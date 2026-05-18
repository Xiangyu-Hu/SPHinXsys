/**
 * @file test_3d_em_aphi_matrix_free_gradient_manufactured.cpp
 * @brief Manufactured verification for matrix-free complex gradient operator.
 *
 * Validates:
 *   1) grad(constant) = 0
 *   2) grad(linear) = constant exact complex vector
 *   3) grad(sinusoidal) matches analytical gradient in the core region
 *
 * Optional SYCL path:
 *   EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_gradient.h"

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

class GradientBoxShape : public ComplexShape
{
  public:
    GradientBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
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
    return Complex(0.7, -0.25);
}

Complex linearFieldValue(const Vecd &position, const Vecd &center)
{
    const Complex gx(0.60, -0.10);
    const Complex gy(-0.20, 0.40);
    const Complex gz(0.15, 0.25);
    return Complex(0.35, -0.45) + gx * (position[0] - center[0]) + gy * (position[1] - center[1]) +
           gz * (position[2] - center[2]);
}

Vec3c linearFieldGradientExact()
{
    return makeVec3c(Complex(0.60, -0.10), Complex(-0.20, 0.40), Complex(0.15, 0.25));
}

Complex sinusoidalFieldValue(const Vecd &position, Real length, Real height, Real width)
{
    const Complex ax(1.10, 0.20);
    const Complex ay(-0.45, 0.35);
    const Complex az(0.30, -0.55);
    const Real kx = Pi / length;
    const Real ky = 2.0 * Pi / height;
    const Real kz = Pi / width;
    return ax * std::sin(kx * position[0]) + ay * std::cos(ky * position[1]) + az * std::sin(kz * position[2]);
}

Vec3c sinusoidalFieldGradientExact(const Vecd &position, Real length, Real height, Real width)
{
    const Complex ax(1.10, 0.20);
    const Complex ay(-0.45, 0.35);
    const Complex az(0.30, -0.55);
    const Real kx = Pi / length;
    const Real ky = 2.0 * Pi / height;
    const Real kz = Pi / width;
    return makeVec3c(ax * kx * std::cos(kx * position[0]), -ay * ky * std::sin(ky * position[1]),
                     az * kz * std::cos(kz * position[2]));
}

struct Summary
{
    size_t total_particles = 0;
    size_t core_particles = 0;
    Real constant_max_grad_abs = 0.0;
    Real linear_mean_error = 0.0;
    Real linear_max_error = 0.0;
    Real linear_core_mean_error = 0.0;
    Real linear_core_max_error = 0.0;
    Real sinusoid_mean_error = 0.0;
    Real sinusoid_max_error = 0.0;
    Real sinusoid_core_mean_error = 0.0;
    Real sinusoid_core_max_error = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_DP", 0.10);
    const Real body_length = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_LENGTH", 1.0);
    const Real body_height = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_HEIGHT", 1.0);
    const Real body_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_WIDTH", 1.0);
    const Real boundary_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Real core_shell = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_CORE_SHELL", 2.5 * dp_0);
    const Real validation_constant_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_MAX_CONSTANT_ABS", 1.0e-8);
    const Real validation_linear_core_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_MAX_LINEAR_CORE", 2.0e-1);
    const Real validation_sinusoid_core_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_GRADIENT_MAX_SINUSOID_CORE", 3.5e-1);
    const bool use_sycl_gradient = getEnvBoolLocal("EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT", false);
    const bool use_sph_native_gradient = getEnvBoolLocal("EM_APHI_USE_SPH_NATIVE_GRADIENT", false);

    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system, makeShared<GradientBoxShape>("MatrixFreeGradientBody", body_center, body_halfsize));
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
    MatrixFreeAPhiDiscreteView discrete_view{number_of_particles,
                                             body.getSPHAdaptation().ReferenceSmoothingLength(),
                                             positions,
                                             volumetric_measure,
                                             nullptr,
                                             &inner_relation.inner_configuration_,
                                             nullptr};
    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);

#if SPHINXSYS_USE_SYCL
    setMatrixFreeGradientUseSycl(use_sycl_gradient);
    if (use_sycl_gradient)
    {
        prepareMatrixFreeAPhiSyclDeviceResources(graph, number_of_particles);
    }
#else
    (void)use_sycl_gradient;
#endif

    StdVec<Complex> constant_field(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> linear_field(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> sinusoid_field(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        constant_field[i] = constantFieldValue();
        linear_field[i] = linearFieldValue(positions[i], body_center);
        sinusoid_field[i] = sinusoidalFieldValue(positions[i], body_length, body_height, body_width);
    }

    std::unique_ptr<SPHComplexScalarGradientOperator> sph_native_gradient_operator;
    if (use_sph_native_gradient)
    {
        sph_native_gradient_operator = std::make_unique<SPHComplexScalarGradientOperator>(body, inner_ck);
    }

    const auto compute_gradient = [&](const StdVec<Complex> &scalar_field) -> StdVec<Vec3c>
    {
        if (use_sph_native_gradient)
        {
            return sph_native_gradient_operator->computeFromField(scalar_field);
        }
        return applyMatrixFreeGradient(graph, scalar_field);
    };

    const StdVec<Vec3c> constant_gradient = compute_gradient(constant_field);
    const StdVec<Vec3c> linear_gradient = compute_gradient(linear_field);
    const StdVec<Vec3c> sinusoid_gradient = compute_gradient(sinusoid_field);

    Summary summary;
    const Vec3c linear_gradient_exact = linearFieldGradientExact();
    Real linear_error_sum = 0.0;
    Real linear_core_error_sum = 0.0;
    Real sinusoid_error_sum = 0.0;
    Real sinusoid_core_error_sum = 0.0;

    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const bool is_core = distanceToBoundary(positions[i], body_length, body_height, body_width) > core_shell;
        const Real constant_grad_abs = complexVecAbsNorm(constant_gradient[i]);
        const Real linear_error = complexVecErrorNorm(linear_gradient[i], linear_gradient_exact);
        const Real sinusoid_error =
            complexVecErrorNorm(sinusoid_gradient[i], sinusoidalFieldGradientExact(positions[i], body_length, body_height, body_width));

        summary.total_particles++;
        summary.constant_max_grad_abs = SMAX(summary.constant_max_grad_abs, constant_grad_abs);
        linear_error_sum += linear_error;
        sinusoid_error_sum += sinusoid_error;
        summary.linear_max_error = SMAX(summary.linear_max_error, linear_error);
        summary.sinusoid_max_error = SMAX(summary.sinusoid_max_error, sinusoid_error);

        if (is_core)
        {
            summary.core_particles++;
            linear_core_error_sum += linear_error;
            sinusoid_core_error_sum += sinusoid_error;
            summary.linear_core_max_error = SMAX(summary.linear_core_max_error, linear_error);
            summary.sinusoid_core_max_error = SMAX(summary.sinusoid_core_max_error, sinusoid_error);
        }
    }

    summary.linear_mean_error = linear_error_sum / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.sinusoid_mean_error = sinusoid_error_sum / (static_cast<Real>(summary.total_particles) + TinyReal);
    summary.linear_core_mean_error = linear_core_error_sum / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.sinusoid_core_mean_error = sinusoid_core_error_sum / (static_cast<Real>(summary.core_particles) + TinyReal);

    bool validation_pass = true;
    validation_pass = validation_pass && summary.constant_max_grad_abs <= validation_constant_max;
    validation_pass = validation_pass && summary.linear_core_max_error <= validation_linear_core_max;
    validation_pass = validation_pass && summary.sinusoid_core_max_error <= validation_sinusoid_core_max;

    std::cout << std::setprecision(12)
              << "matrix_free_gradient_manufactured"
              << " points=" << number_of_particles
              << " dp=" << dp_0
              << " sycl_gradient=" << (use_sycl_gradient ? 1 : 0)
              << " sph_native_gradient=" << (use_sph_native_gradient ? 1 : 0)
              << " gradient_backend="
              << (use_sph_native_gradient ? "sph_native_b_corrected" : (use_sycl_gradient ? "graph_sycl" : "graph_cpu"))
              << " core_particles=" << summary.core_particles
              << " constant_max_grad_abs=" << summary.constant_max_grad_abs
              << " linear_mean_error=" << summary.linear_mean_error
              << " linear_max_error=" << summary.linear_max_error
              << " linear_core_mean_error=" << summary.linear_core_mean_error
              << " linear_core_max_error=" << summary.linear_core_max_error
              << " sinusoid_mean_error=" << summary.sinusoid_mean_error
              << " sinusoid_max_error=" << summary.sinusoid_max_error
              << " sinusoid_core_mean_error=" << summary.sinusoid_core_mean_error
              << " sinusoid_core_max_error=" << summary.sinusoid_core_max_error
              << " validation_max_constant_abs=" << validation_constant_max
              << " validation_max_linear_core=" << validation_linear_core_max
              << " validation_max_sinusoid_core=" << validation_sinusoid_core_max
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << std::endl;

    return validation_pass ? 0 : 1;
}
