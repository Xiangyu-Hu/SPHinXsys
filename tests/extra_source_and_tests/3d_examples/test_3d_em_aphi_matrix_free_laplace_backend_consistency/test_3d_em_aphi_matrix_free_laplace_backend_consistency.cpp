/**
 * @file test_3d_em_aphi_matrix_free_laplace_backend_consistency.cpp
 * @brief Manufactured and backend-consistency verification for matrix-free negative Laplace operator.
 *
 * Validates:
 *   1) -Laplace(constant) = 0
 *   2) -Laplace(linear) = 0 in the core
 *   3) -Laplace(sinusoidal) matches analytical value in the core
 *   4) current backend (CPU or SYCL) matches CPU reference
 *
 * Optional SYCL path:
 *   EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_laplace.h"

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

class LaplaceBoxShape : public ComplexShape
{
  public:
    LaplaceBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Real complexAbsError(const Complex &lhs, const Complex &rhs)
{
    return std::abs(lhs - rhs);
}

Real distanceToBoundary(const Vecd &position, Real length, Real height, Real width)
{
    return std::min({position[0], length - position[0], position[1], height - position[1], position[2], width - position[2]});
}

Complex constantFieldValue()
{
    return Complex(-0.4, 0.6);
}

Complex linearFieldValue(const Vecd &position, const Vecd &center)
{
    const Complex gx(0.52, -0.08);
    const Complex gy(-0.27, 0.16);
    const Complex gz(0.11, 0.19);
    return Complex(0.15, -0.35) + gx * (position[0] - center[0]) + gy * (position[1] - center[1]) +
           gz * (position[2] - center[2]);
}

Complex sinusoidalFieldValue(const Vecd &position, Real length, Real height, Real width)
{
    const Complex ax(0.90, 0.20);
    const Complex ay(-0.55, 0.30);
    const Complex az(0.25, -0.45);
    const Real kx = 0.5 * Pi / length;
    const Real ky = 0.5 * Pi / height;
    const Real kz = 0.5 * Pi / width;
    return ax * std::sin(kx * position[0]) + ay * std::cos(ky * position[1]) + az * std::sin(kz * position[2]);
}

Complex sinusoidalNegativeLaplaceExact(const Vecd &position, Real length, Real height, Real width)
{
    const Complex ax(0.90, 0.20);
    const Complex ay(-0.55, 0.30);
    const Complex az(0.25, -0.45);
    const Real kx = 0.5 * Pi / length;
    const Real ky = 0.5 * Pi / height;
    const Real kz = 0.5 * Pi / width;
    return ax * (kx * kx) * std::sin(kx * position[0]) + ay * (ky * ky) * std::cos(ky * position[1]) +
           az * (kz * kz) * std::sin(kz * position[2]);
}

struct Summary
{
    size_t total_particles = 0;
    size_t core_particles = 0;
    Real constant_max_abs = 0.0;
    Real linear_core_mean_error = 0.0;
    Real linear_core_max_error = 0.0;
    Real sinusoid_core_mean_error = 0.0;
    Real sinusoid_core_max_error = 0.0;
    Real backend_constant_max_difference = 0.0;
    Real backend_linear_max_difference = 0.0;
    Real backend_sinusoid_max_difference = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_DP", 0.10);
    const Real body_length = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_LENGTH", 1.0);
    const Real body_height = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_HEIGHT", 1.0);
    const Real body_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_WIDTH", 1.0);
    const Real boundary_width = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Real core_shell = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_CORE_SHELL", 2.5 * dp_0);
    const Real validation_constant_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_MAX_CONSTANT_ABS", 1.0e-8);
    const Real validation_linear_core_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_MAX_LINEAR_CORE", 2.5e-1);
    const Real validation_sinusoid_core_max = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_MAX_SINUSOID_CORE", 6.0e-1);
    const Real validation_backend_max_difference = getEnvRealLocal("EM_APHI_MATRIX_FREE_LAPLACE_MAX_BACKEND_DIFF", 1.0e-5);
    const bool use_sycl_laplace = getEnvBoolLocal("EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL", false);
    const bool use_sph_native_laplace = getEnvBoolLocal("EM_APHI_USE_SPH_NATIVE_LAPLACE", false);

    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system, makeShared<LaplaceBoxShape>("MatrixFreeLaplaceBody", body_center, body_halfsize));
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
    const StdVec<Real> diffusion(number_of_particles, 1.0);

    StdVec<Complex> constant_field(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> linear_field(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> sinusoid_field(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        constant_field[i] = constantFieldValue();
        linear_field[i] = linearFieldValue(positions[i], body_center);
        sinusoid_field[i] = sinusoidalFieldValue(positions[i], body_length, body_height, body_width);
    }

#if SPHINXSYS_USE_SYCL
    setMatrixFreeLaplaceResidualUseSycl(false);
    if (use_sycl_laplace)
    {
        prepareMatrixFreeAPhiSyclDeviceResources(graph, number_of_particles);
    }
#else
    (void)use_sycl_laplace;
#endif

    std::unique_ptr<SPHComplexScalarNegativeLaplaceOperator> sph_native_laplace_operator;
    if (use_sph_native_laplace)
    {
        SPHScalarNegativeLaplaceParameters sph_laplace_parameters;
        sph_laplace_parameters.pair_weight_regularization = parameters.pair_weight_regularization;
        sph_native_laplace_operator =
            std::make_unique<SPHComplexScalarNegativeLaplaceOperator>(body, inner_ck, sph_laplace_parameters);
    }

    const auto apply_negative_laplace = [&](const StdVec<Complex> &scalar_field) -> StdVec<Complex>
    {
        if (use_sph_native_laplace)
        {
            return sph_native_laplace_operator->apply(scalar_field, diffusion);
        }
#if SPHINXSYS_USE_SYCL
        if (use_sycl_laplace)
        {
            setMatrixFreeLaplaceResidualUseSycl(true);
            const StdVec<Complex> result = applyScalarNegativeLaplaceFromGraph(graph, scalar_field, diffusion);
            setMatrixFreeLaplaceResidualUseSycl(false);
            return result;
        }
#endif
        return applyScalarNegativeLaplaceFromGraph(graph, scalar_field, diffusion);
    };

    const StdVec<Complex> cpu_constant = applyScalarNegativeLaplaceFromGraph(graph, constant_field, diffusion);
    const StdVec<Complex> cpu_linear = applyScalarNegativeLaplaceFromGraph(graph, linear_field, diffusion);
    const StdVec<Complex> cpu_sinusoid = applyScalarNegativeLaplaceFromGraph(graph, sinusoid_field, diffusion);

    StdVec<Complex> backend_constant = apply_negative_laplace(constant_field);
    StdVec<Complex> backend_linear = apply_negative_laplace(linear_field);
    StdVec<Complex> backend_sinusoid = apply_negative_laplace(sinusoid_field);

    Summary summary;
    Real linear_core_error_sum = 0.0;
    Real sinusoid_core_error_sum = 0.0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const bool is_core = distanceToBoundary(positions[i], body_length, body_height, body_width) > core_shell;
        const Real constant_abs = std::abs(backend_constant[i]);
        const Real backend_constant_difference = complexAbsError(backend_constant[i], cpu_constant[i]);
        const Real backend_linear_difference = complexAbsError(backend_linear[i], cpu_linear[i]);
        const Real backend_sinusoid_difference = complexAbsError(backend_sinusoid[i], cpu_sinusoid[i]);
        const Real linear_error = std::abs(backend_linear[i]);
        const Real sinusoid_error =
            complexAbsError(backend_sinusoid[i], sinusoidalNegativeLaplaceExact(positions[i], body_length, body_height, body_width));

        summary.total_particles++;
        summary.constant_max_abs = SMAX(summary.constant_max_abs, constant_abs);
        summary.backend_constant_max_difference = SMAX(summary.backend_constant_max_difference, backend_constant_difference);
        summary.backend_linear_max_difference = SMAX(summary.backend_linear_max_difference, backend_linear_difference);
        summary.backend_sinusoid_max_difference = SMAX(summary.backend_sinusoid_max_difference, backend_sinusoid_difference);
        if (is_core)
        {
            summary.core_particles++;
            linear_core_error_sum += linear_error;
            sinusoid_core_error_sum += sinusoid_error;
            summary.linear_core_max_error = SMAX(summary.linear_core_max_error, linear_error);
            summary.sinusoid_core_max_error = SMAX(summary.sinusoid_core_max_error, sinusoid_error);
        }
    }

    summary.linear_core_mean_error = linear_core_error_sum / (static_cast<Real>(summary.core_particles) + TinyReal);
    summary.sinusoid_core_mean_error = sinusoid_core_error_sum / (static_cast<Real>(summary.core_particles) + TinyReal);

    const Real validation_backend_max_difference_effective =
        use_sph_native_laplace ? getEnvRealLocal("EM_APHI_SPH_NATIVE_LAPLACE_MAX_GRAPH_DIFF", 1.0e-2)
                               : validation_backend_max_difference;

    bool validation_pass = true;
    validation_pass = validation_pass && summary.core_particles > 0;
    validation_pass = validation_pass && summary.constant_max_abs <= validation_constant_max;
    validation_pass = validation_pass && summary.linear_core_max_error <= validation_linear_core_max;
    validation_pass = validation_pass && summary.sinusoid_core_max_error <= validation_sinusoid_core_max;
    validation_pass = validation_pass && summary.backend_constant_max_difference <= validation_backend_max_difference_effective;
    validation_pass = validation_pass && summary.backend_linear_max_difference <= validation_backend_max_difference_effective;
    validation_pass = validation_pass && summary.backend_sinusoid_max_difference <= validation_backend_max_difference_effective;

    std::cout << std::setprecision(12)
              << "matrix_free_laplace_backend_consistency"
              << " points=" << number_of_particles
              << " dp=" << dp_0
              << " sycl_laplace=" << (use_sycl_laplace ? 1 : 0)
              << " sph_native_laplace=" << (use_sph_native_laplace ? 1 : 0)
              << " laplace_backend="
              << (use_sph_native_laplace ? "sph_native_ck" : (use_sycl_laplace ? "graph_sycl" : "graph_cpu"))
              << " core_particles=" << summary.core_particles
              << " constant_max_abs=" << summary.constant_max_abs
              << " linear_core_mean_error=" << summary.linear_core_mean_error
              << " linear_core_max_error=" << summary.linear_core_max_error
              << " sinusoid_core_mean_error=" << summary.sinusoid_core_mean_error
              << " sinusoid_core_max_error=" << summary.sinusoid_core_max_error
              << " backend_constant_max_difference=" << summary.backend_constant_max_difference
              << " backend_linear_max_difference=" << summary.backend_linear_max_difference
              << " backend_sinusoid_max_difference=" << summary.backend_sinusoid_max_difference
              << " validation_max_constant_abs=" << validation_constant_max
              << " validation_max_linear_core=" << validation_linear_core_max
              << " validation_max_sinusoid_core=" << validation_sinusoid_core_max
              << " validation_max_backend_difference=" << validation_backend_max_difference_effective
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << std::endl;

    return validation_pass ? 0 : 1;
}
