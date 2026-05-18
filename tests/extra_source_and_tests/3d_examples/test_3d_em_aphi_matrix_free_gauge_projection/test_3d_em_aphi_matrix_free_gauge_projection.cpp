/**
 * @file test_3d_em_aphi_matrix_free_gauge_projection.cpp
 * @brief First relation-based matrix-free gauge-projection prototype.
 * The current version uses a thin body with a real SPH inner relation and a matrix-free pairwise graph.
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_gauge_projection.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_pairwise_graph.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
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

size_t get_env_size_t_local(const std::string &name, size_t default_value)
{
    const Real parsed = get_env_real_local(name, static_cast<Real>(default_value));
    return parsed < 1.0 ? default_value : static_cast<size_t>(parsed);
}

std::string get_env_string_local(const std::string &name, const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

class GaugeBoxShape : public ComplexShape
{
  public:
    GaugeBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Complex exact_chi_value(Real x, Real length, const Complex &scale)
{
    return scale * std::sin(Pi * x / length);
}

Complex background_ax_value(const Complex &scale)
{
    return scale;
}

void write_summary_csv(const std::string &filename, size_t number_of_points, Real length, Real dp,
                       Real angular_frequency, const ScalarComplexHelmholtzSolverState &solver_state,
                       const LineGaugeProjectionDiagnostics &diagnostics)
{
    std::ofstream output(filename);
    output << "number_of_points,length,dp,angular_frequency,iterations,converged,initial_residual_l2,"
              "final_residual_l2,chi_l2_error,chi_max_error,projection_before_l2,projection_after_l2,"
              "gradient_div_a_before_l2,gradient_div_a_after_l2,electric_field_before_l2,electric_field_after_l2,"
              "electric_field_change_l2\n";
    output << number_of_points << "," << length << "," << dp << "," << angular_frequency << ","
           << solver_state.iterations_ << "," << (solver_state.converged_ ? 1 : 0) << ","
           << solver_state.initial_residual_l2_ << "," << solver_state.current_residual_l2_ << ","
           << diagnostics.chi_l2_error_ << "," << diagnostics.chi_max_error_ << ","
           << diagnostics.projection_before_l2_ << "," << diagnostics.projection_after_l2_ << ","
           << diagnostics.gradient_div_a_before_l2_ << "," << diagnostics.gradient_div_a_after_l2_ << ","
           << diagnostics.electric_field_before_l2_ << "," << diagnostics.electric_field_after_l2_ << ","
           << diagnostics.electric_field_change_l2_ << "\n";
}

void write_profile_csv(const std::string &filename, const StdVec<Real> &x, const StdVec<Complex> &exact_chi,
                       const StdVec<Complex> &solved_chi, const StdVec<Complex> &projection_before,
                       const StdVec<Complex> &projection_after, const StdVec<Complex> &divergence_before,
                       const StdVec<Complex> &divergence_after, const StdVec<Complex> &ax_before,
                       const StdVec<Complex> &ax_after, const StdVec<Complex> &electric_field_before,
                       const StdVec<Complex> &electric_field_after)
{
    std::ofstream output(filename);
    output << "x,exact_chi_real,exact_chi_imag,solved_chi_real,solved_chi_imag,projection_before_real,"
              "projection_before_imag,projection_after_real,projection_after_imag,div_before_real,div_before_imag,"
              "div_after_real,div_after_imag,ax_before_real,ax_before_imag,ax_after_real,ax_after_imag,"
              "e_before_real,e_before_imag,e_after_real,e_after_imag\n";
    for (size_t i = 0; i != x.size(); ++i)
    {
        output << x[i] << "," << exact_chi[i].real() << "," << exact_chi[i].imag() << ","
               << solved_chi[i].real() << "," << solved_chi[i].imag() << "," << projection_before[i].real() << ","
               << projection_before[i].imag() << "," << projection_after[i].real() << ","
               << projection_after[i].imag() << "," << divergence_before[i].real() << ","
               << divergence_before[i].imag() << "," << divergence_after[i].real() << ","
               << divergence_after[i].imag() << "," << ax_before[i].real() << "," << ax_before[i].imag() << ","
               << ax_after[i].real() << "," << ax_after[i].imag() << "," << electric_field_before[i].real() << ","
               << electric_field_before[i].imag() << "," << electric_field_after[i].real() << ","
               << electric_field_after[i].imag() << "\n";
    }
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_DP", 0.05);
    const Real body_length = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_LENGTH", 1.0);
    const Real body_height = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_HEIGHT", dp_0);
    const Real body_width = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_WIDTH", dp_0);
    const Real boundary_width = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_BOUNDARY_WIDTH", 3.0 * dp_0);
    const Real angular_frequency = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_OMEGA", 5.0);
    const Complex chi_scale(get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_CHI_REAL", 1.0),
                            get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_CHI_IMAG", 0.2));
    const Complex background_ax(get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_AX0_REAL", 0.15),
                                get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_AX0_IMAG", -0.05));
    const std::string output_tag = get_env_string_local("EM_APHI_MATRIX_FREE_GAUGE_OUTPUT_TAG", "default");

    ScalarComplexHelmholtzSolverParameters solver_parameters;
    solver_parameters.max_iterations_ = get_env_size_t_local("EM_APHI_MATRIX_FREE_GAUGE_MAX_ITERATIONS", 5000);
    solver_parameters.relaxation_factor_ = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_RELAXATION", 1.0);
    solver_parameters.absolute_tolerance_ = get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_ABS_TOL", 1.0e-5);
    solver_parameters.diagonal_regularization_ =
        get_env_real_local("EM_APHI_MATRIX_FREE_GAUGE_DIAGONAL_REGULARIZATION", 1.0e-12);

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

    SolidBody body(sph_system, makeShared<GaugeBoxShape>("GaugeBody", body_center, body_halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t number_of_points = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    MatrixFreeAPhiDiscreteView discrete_view{
        number_of_points,
        body.getSPHAdaptation().ReferenceSmoothingLength(),
        positions,
        volumetric_measure,
        nullptr,
        &inner_relation.inner_configuration_,
        nullptr};

    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);
    StdVec<Real> unit_diffusion(number_of_points, 1.0);

    StdVec<Real> x(number_of_points, 0.0);
    StdVec<Complex> exact_chi(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> ax_before(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> phi_before(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> solved_chi(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> rhs(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> reaction_coefficient(number_of_points, Complex(0.0, 0.0));

    for (size_t i = 0; i != number_of_points; ++i)
    {
        x[i] = positions[i][0];
        exact_chi[i] = exact_chi_value(x[i], body_length, chi_scale);
    }

    const StdVec<Vec3c> exact_grad_chi = applyMatrixFreeGradient(graph, exact_chi);
    for (size_t i = 0; i != number_of_points; ++i)
    {
        ax_before[i] = background_ax_value(background_ax) + exact_grad_chi[i][0];
    }

    const StdVec<Complex> projection_before = applyScalarNegativeLaplaceFromGraph(graph, exact_chi, unit_diffusion);
    rhs = projection_before;

    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(number_of_points);

    const ScalarComplexHelmholtzSolverState solver_state = solveScalarComplexHelmholtz(
        solved_chi, rhs, reaction_coefficient, residuals, solver_parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromGraph(graph, current_field, unit_diffusion, current_residuals);
        });

    Complex mean_offset(0.0, 0.0);
    if (number_of_points != 0)
    {
        for (size_t i = 0; i != number_of_points; ++i)
        {
            mean_offset += solved_chi[i] - exact_chi[i];
        }
        mean_offset /= static_cast<Real>(number_of_points);
        for (size_t i = 0; i != number_of_points; ++i)
        {
            solved_chi[i] -= mean_offset;
        }
    }

    StdVec<Complex> ax_after = ax_before;
    StdVec<Complex> phi_after = phi_before;
    const StdVec<Vec3c> solved_grad_chi = applyMatrixFreeGradient(graph, solved_chi);
    const Complex imaginary_unit(0.0, 1.0);
    for (size_t i = 0; i != number_of_points; ++i)
    {
        ax_after[i] -= solved_grad_chi[i][0];
        phi_after[i] += imaginary_unit * angular_frequency * solved_chi[i];
    }

    const StdVec<Complex> solved_projection = applyScalarNegativeLaplaceFromGraph(graph, solved_chi, unit_diffusion);
    StdVec<Complex> projection_after(number_of_points, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_points; ++i)
    {
        projection_after[i] = rhs[i] - solved_projection[i];
    }

    const StdVec<Vec3c> grad_ax_before = applyMatrixFreeGradient(graph, ax_before);
    const StdVec<Vec3c> grad_ax_after = applyMatrixFreeGradient(graph, ax_after);
    StdVec<Complex> divergence_before(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> divergence_after(number_of_points, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_points; ++i)
    {
        divergence_before[i] = grad_ax_before[i][0];
        divergence_after[i] = grad_ax_after[i][0];
    }

    const StdVec<Vec3c> grad_phi_before = applyMatrixFreeGradient(graph, phi_before);
    const StdVec<Vec3c> grad_phi_after = applyMatrixFreeGradient(graph, phi_after);
    StdVec<Complex> electric_field_before(number_of_points, Complex(0.0, 0.0));
    StdVec<Complex> electric_field_after(number_of_points, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_points; ++i)
    {
        electric_field_before[i] = -imaginary_unit * angular_frequency * ax_before[i] - grad_phi_before[i][0];
        electric_field_after[i] = -imaginary_unit * angular_frequency * ax_after[i] - grad_phi_after[i][0];
    }

    const LineGaugeProjectionDiagnostics diagnostics = evaluateLineGaugeProjectionDiagnostics(
        exact_chi, solved_chi, projection_before, projection_after, divergence_before, divergence_after,
        electric_field_before, electric_field_after);

    const std::string summary_filename = "em_aphi_matrix_free_gauge_projection_summary_" + output_tag + ".csv";
    const std::string profile_filename = "em_aphi_matrix_free_gauge_projection_profile_" + output_tag + ".csv";

    write_summary_csv(summary_filename, number_of_points, body_length, dp_0, angular_frequency, solver_state,
                      diagnostics);
    write_profile_csv(profile_filename, x, exact_chi, solved_chi, projection_before, projection_after,
                      divergence_before, divergence_after, ax_before, ax_after, electric_field_before,
                      electric_field_after);

    std::cout << std::setprecision(12)
              << "matrix_free_gauge_projection_summary"
              << " points=" << number_of_points
              << " dp=" << dp_0
              << " omega=" << angular_frequency
              << " iterations=" << solver_state.iterations_
              << " converged=" << (solver_state.converged_ ? 1 : 0)
              << " initial_residual_l2=" << solver_state.initial_residual_l2_
              << " final_residual_l2=" << solver_state.current_residual_l2_
              << " chi_l2_error=" << diagnostics.chi_l2_error_
              << " projection_before_l2=" << diagnostics.projection_before_l2_
              << " projection_after_l2=" << diagnostics.projection_after_l2_
              << " divergence_before_l2=" << diagnostics.gradient_div_a_before_l2_
              << " divergence_after_l2=" << diagnostics.gradient_div_a_after_l2_
              << " electric_field_change_l2=" << diagnostics.electric_field_change_l2_
              << std::endl;

    return solver_state.iterations_ == 0 ? 1 : 0;
}
