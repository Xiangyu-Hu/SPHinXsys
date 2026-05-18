/**
 * @file test_3d_em_aphi_matrix_free_complex_helmholtz.cpp
 * @brief First relation-based matrix-free verification case for a complex scalar Helmholtz solve.
 * The current version uses a simple box body with a real SPH inner relation and a matrix-free pairwise graph.
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_solver.h"
#include "aphi_sphinxsys/electromagnetic_aphi_sph_native_context.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;

namespace
{

Real get_env_real_local(const std::string &name, Real default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    char *end_ptr = nullptr;
    Real parsed = static_cast<Real>(std::strtod(value, &end_ptr));
    if (end_ptr == value || !std::isfinite(parsed))
    {
        return default_value;
    }
    return parsed;
}

std::string get_env_string_local(const std::string &name, const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

bool get_env_bool_local(const std::string &name, bool default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    return std::atoi(value) != 0;
}

enum class HelmholtzRhsMode
{
    discrete_manufactured,
    continuum_manufactured
};

HelmholtzRhsMode parse_rhs_mode(const std::string &token)
{
    std::string normalized = token;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (normalized == "continuum" || normalized == "continuum_manufactured")
    {
        return HelmholtzRhsMode::continuum_manufactured;
    }
    return HelmholtzRhsMode::discrete_manufactured;
}

class HelmholtzBoxShape : public ComplexShape
{
  public:
    HelmholtzBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Complex exact_solution(Real x, Real length)
{
    return Complex(1.0, 0.3) * std::sin(Pi * x / length);
}

Complex continuum_rhs(Real x, Real length, const Complex &alpha)
{
    const Real wave_number = Pi / length;
    return (wave_number * wave_number) * exact_solution(x, length) + alpha * exact_solution(x, length);
}

void write_summary_csv(const std::string &filename, size_t number_of_points, Real length, Real dp,
                       const Complex &alpha, const ScalarComplexHelmholtzSolverState &state,
                       Real l2_error, Real max_error, const std::string &rhs_mode)
{
    std::ofstream output(filename);
    output << "number_of_points,length,dp,alpha_real,alpha_imag,iterations,converged,initial_residual_l2,"
              "final_residual_l2,final_mean_abs,final_max_abs,l2_error,max_error,rhs_mode\n";
    output << number_of_points << "," << length << "," << dp << "," << alpha.real() << "," << alpha.imag() << ","
           << state.iterations_ << "," << (state.converged_ ? 1 : 0) << "," << state.initial_residual_l2_ << ","
           << state.current_residual_l2_ << "," << state.current_mean_abs_ << "," << state.current_max_abs_ << ","
           << l2_error << "," << max_error << "," << rhs_mode << "\n";
}

void write_profile_csv(const std::string &filename, const StdVec<Real> &x, const StdVec<Complex> &exact_field,
                       const StdVec<Complex> &numerical_field)
{
    std::ofstream output(filename);
    output << "x,exact_real,exact_imag,solved_real,solved_imag,abs_error\n";
    for (size_t i = 0; i != x.size(); ++i)
    {
        output << x[i] << "," << exact_field[i].real() << "," << exact_field[i].imag() << ","
               << numerical_field[i].real() << "," << numerical_field[i].imag() << ","
               << std::abs(numerical_field[i] - exact_field[i]) << "\n";
    }
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_DP", 0.05);
    const Real body_length = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_LENGTH", 1.0);
    const Real body_height = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_HEIGHT", dp_0);
    const Real body_width = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_WIDTH", dp_0);
    const Real boundary_width = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Complex alpha(get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_ALPHA_REAL", 1.0),
                        get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_ALPHA_IMAG", 2.0));
    const HelmholtzRhsMode rhs_mode =
        parse_rhs_mode(get_env_string_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_RHS_MODE", "discrete_manufactured"));
    const std::string output_tag =
        get_env_string_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_OUTPUT_TAG", "default");

    const Vecd body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                              Vecd(body_length + boundary_width, body_height + boundary_width,
                                   body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    SolidBody body(sph_system, makeShared<HelmholtzBoxShape>("HelmholtzBody", body_center, body_halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);
    Inner<> inner_ck(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t number_of_points = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    parameters.max_change_rate = MaxReal;

    MatrixFreeAPhiDiscreteView discrete_view{
        number_of_points,
        body.getSPHAdaptation().ReferenceSmoothingLength(),
        positions,
        volumetric_measure,
        nullptr,
        &inner_relation.inner_configuration_,
        nullptr};

    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);

    const bool use_sph_native_helmholtz = get_env_bool_local("EM_APHI_USE_SPH_NATIVE_HELMHOLTZ", false) ||
                                          get_env_bool_local("EM_APHI_USE_SPH_NATIVE_RESIDUALS", false);
    std::unique_ptr<MatrixFreeAPhiSphNativeContext> sph_native_context_storage;
    if (use_sph_native_helmholtz)
    {
        sph_native_context_storage = std::make_unique<MatrixFreeAPhiSphNativeContext>(body, inner_ck);
    }
    MatrixFreeAPhiSphNativeContext *sph_native_context = sph_native_context_storage.get();

    StdVec<Real> x(number_of_points, 0.0);
    StdVec<Complex> exact_field(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> field(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> rhs(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> reaction_coefficient(number_of_points, alpha);
    StdVec<Real> diffusion_coefficient(number_of_points, 1.0);

    for (size_t i = 0; i != number_of_points; ++i)
    {
        x[i] = positions[i][0];
        exact_field[i] = exact_solution(x[i], body_length);
    }

    if (rhs_mode == HelmholtzRhsMode::discrete_manufactured)
    {
        ScalarComplexHelmholtzResiduals exact_residuals;
        exact_residuals.resize(number_of_points);
        exact_residuals.clear();
        if (use_sph_native_helmholtz)
        {
            sph_native_context->operatorAssembly().accumulateScalarLaplaceHelmholtzResiduals(
                exact_field, diffusion_coefficient, exact_residuals);
        }
        else
        {
            accumulateScalarLaplaceResidualsFromGraph(graph, exact_field, diffusion_coefficient, exact_residuals);
        }
        for (size_t i = 0; i != number_of_points; ++i)
        {
            rhs[i] = exact_residuals.laplace_term_[i] + alpha * exact_field[i];
        }
    }
    else
    {
        for (size_t i = 0; i != number_of_points; ++i)
        {
            rhs[i] = continuum_rhs(x[i], body_length, alpha);
        }
    }

    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(number_of_points);

    ScalarComplexHelmholtzSolverParameters solver_parameters;
    solver_parameters.max_iterations_ = static_cast<size_t>(get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_MAX_ITERATIONS", 5000));
    solver_parameters.relaxation_factor_ = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_RELAXATION", 1.0);
    solver_parameters.absolute_tolerance_ = get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_ABS_TOL", 1.0e-5);
    solver_parameters.diagonal_regularization_ =
        get_env_real_local("EM_APHI_MATRIX_FREE_HELMHOLTZ_DIAGONAL_REGULARIZATION", 1.0e-12);

    const ScalarComplexHelmholtzSolverState state = solveScalarComplexHelmholtz(
        field, rhs, reaction_coefficient, residuals, solver_parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            if (use_sph_native_helmholtz)
            {
                sph_native_context->operatorAssembly().accumulateScalarLaplaceHelmholtzResiduals(
                    current_field, diffusion_coefficient, current_residuals);
            }
            else
            {
                accumulateScalarLaplaceResidualsFromGraph(graph, current_field, diffusion_coefficient, current_residuals);
            }
        });

    Real squared_error_sum = 0.0;
    Real max_error = 0.0;
    for (size_t i = 0; i != number_of_points; ++i)
    {
        const Real abs_error = std::abs(field[i] - exact_field[i]);
        squared_error_sum += abs_error * abs_error;
        max_error = std::max(max_error, abs_error);
    }
    const Real l2_error = number_of_points == 0 ? 0.0 : std::sqrt(squared_error_sum / static_cast<Real>(number_of_points));

    const std::string rhs_mode_name = rhs_mode == HelmholtzRhsMode::discrete_manufactured
                                          ? "discrete_manufactured"
                                          : "continuum_manufactured";

    write_summary_csv("em_aphi_matrix_free_complex_helmholtz_summary_" + output_tag + ".csv", number_of_points,
                      body_length, dp_0, alpha, state, l2_error, max_error, rhs_mode_name);
    write_profile_csv("em_aphi_matrix_free_complex_helmholtz_profile_" + output_tag + ".csv", x, exact_field, field);

    std::cout << std::setprecision(12)
              << "matrix_free_helmholtz_summary"
              << " points=" << number_of_points
              << " dp=" << dp_0
              << " alpha_real=" << alpha.real()
              << " alpha_imag=" << alpha.imag()
              << " rhs_mode=" << rhs_mode_name
              << " sph_native_helmholtz=" << (use_sph_native_helmholtz ? 1 : 0)
              << " helmholtz_backend=" << (use_sph_native_helmholtz ? "sph_native_ck" : "graph_cpu")
              << " iterations=" << state.iterations_
              << " converged=" << (state.converged_ ? 1 : 0)
              << " initial_residual_l2=" << state.initial_residual_l2_
              << " final_residual_l2=" << state.current_residual_l2_
              << " l2_error=" << l2_error
              << " max_error=" << max_error
              << std::endl;

    return state.iterations_ == 0 ? 1 : 0;
}
