/**
 * @file test_3d_em_aphi_matrix_free_staggered_smoke.cpp
 * @brief Simple physical-response baseline for the matrix-free staggered A-phi solver.
 */

#include "sphinxsys.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.h"

#include <cstdlib>
#include <iomanip>
#include <string>
#include <Eigen/Dense>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::matrix_free;

namespace
{
class SmokeBoxShape : public ComplexShape
{
  public:
    SmokeBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

Real computeL2Error(const StdVec<Complex> &solution, const StdVec<Complex> &reference)
{
    if (solution.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        squared_sum += std::norm(solution[i] - reference[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

Real computeComplexFieldL2NormLocal(const StdVec<Complex> &field)
{
    if (field.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (const Complex &value : field)
    {
        squared_sum += std::norm(value);
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

StdVec<Complex> computeDivergenceOfVectorFieldLocal(const MatrixFreePairwiseGraph &graph,
                                                    const MatrixFreeAPhiFields &fields)
{
    const size_t number_of_particles = fields.ax.size();
    const StdVec<Vec3c> grad_ax = applyMatrixFreeGradient(graph, fields.ax);
    const StdVec<Vec3c> grad_ay = applyMatrixFreeGradient(graph, fields.ay);
    const StdVec<Vec3c> grad_az = applyMatrixFreeGradient(graph, fields.az);
    StdVec<Complex> divergence(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        divergence[i] = grad_ax[i][0] + grad_ay[i][1] + grad_az[i][2];
    }
    return divergence;
}

Real computeMaxError(const StdVec<Complex> &solution, const StdVec<Complex> &reference)
{
    Real max_error = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        max_error = SMAX(max_error, std::abs(solution[i] - reference[i]));
    }
    return max_error;
}


Real computeVectorFieldL2Error(const StdVec<Vec3c> &solution, const StdVec<Vec3c> &reference)
{
    if (solution.size() != reference.size() || solution.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            squared_sum += std::norm(solution[i][axis] - reference[i][axis]);
        }
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

Real computeRealFieldL2Error(const StdVec<Real> &solution, const StdVec<Real> &reference)
{
    if (solution.size() != reference.size() || solution.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (size_t i = 0; i != solution.size(); ++i)
    {
        const Real difference = solution[i] - reference[i];
        squared_sum += difference * difference;
    }
    return std::sqrt(squared_sum / static_cast<Real>(solution.size()));
}

StdVec<Real> computeJouleHeatingDensity(const StdVec<Vec3c> &electric_field,
                                        const StdVec<Vec3c> &current_density)
{
    const size_t number_of_particles = electric_field.size();
    StdVec<Real> joule_density(number_of_particles, 0.0);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        Complex power_density(0.0, 0.0);
        for (int axis = 0; axis != Dimensions; ++axis)
        {
            power_density += electric_field[i][axis] * std::conj(current_density[i][axis]);
        }
        joule_density[i] = 0.5 * power_density.real();
    }
    return joule_density;
}


Real computeRealFieldL2Norm(const StdVec<Real> &field)
{
    if (field.empty())
    {
        return 0.0;
    }
    Real squared_sum = 0.0;
    for (const Real value : field)
    {
        squared_sum += value * value;
    }
    return std::sqrt(squared_sum / static_cast<Real>(field.size()));
}

int findMirrorIndex(const Vecd *positions, size_t number_of_particles, size_t target_index,
                    Real center_y, Real tolerance)
{
    const Real mirror_y = 2.0 * center_y - positions[target_index][1];
    for (size_t candidate = 0; candidate != number_of_particles; ++candidate)
    {
        if (std::abs(positions[candidate][0] - positions[target_index][0]) <= tolerance &&
            std::abs(positions[candidate][1] - mirror_y) <= tolerance &&
            std::abs(positions[candidate][2] - positions[target_index][2]) <= tolerance)
        {
            return static_cast<int>(candidate);
        }
    }
    return -1;
}

Real computeComplexFieldMirrorSymmetryError(const StdVec<Complex> &field, const Vecd *positions,
                                            size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        squared_sum += std::norm(field[i] - field[static_cast<size_t>(mirror_index)]);
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

Real computeComplexFieldMirrorAntisymmetryError(const StdVec<Complex> &field, const Vecd *positions,
                                                size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        squared_sum += std::norm(field[i] + field[static_cast<size_t>(mirror_index)]);
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

StdVec<Complex> makeMeanCenteredCopy(const StdVec<Complex> &field)
{
    StdVec<Complex> centered = field;
    removeMeanOffset(centered);
    return centered;
}

Real computeRealFieldMirrorSymmetryError(const StdVec<Real> &field, const Vecd *positions,
                                         size_t number_of_particles, Real center_y, Real tolerance)
{
    Real squared_sum = 0.0;
    size_t pair_count = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        if (positions[i][1] > center_y + tolerance)
        {
            continue;
        }
        const int mirror_index = findMirrorIndex(positions, number_of_particles, i, center_y, tolerance);
        if (mirror_index < 0)
        {
            continue;
        }
        const Real difference = field[i] - field[static_cast<size_t>(mirror_index)];
        squared_sum += difference * difference;
        ++pair_count;
    }
    return pair_count == 0 ? 0.0 : std::sqrt(squared_sum / static_cast<Real>(pair_count));
}

size_t readEnvSizeT(const char *name, size_t default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : static_cast<size_t>(std::strtoull(value, nullptr, 10));
}

Real readEnvReal(const char *name, Real default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : static_cast<Real>(std::strtod(value, nullptr));
}

bool readEnvBool(const char *name, bool default_value)
{
    const char *value = std::getenv(name);
    if (value == nullptr)
    {
        return default_value;
    }
    return std::atoi(value) != 0;
}

std::string readEnvString(const char *name, const std::string &default_value)
{
    const char *value = std::getenv(name);
    return value == nullptr ? default_value : std::string(value);
}

void enforceReference(StdVec<Complex> &field, size_t reference_index, Complex reference_value)
{
    if (field.empty() || reference_index >= field.size())
    {
        return;
    }
    const Complex offset = field[reference_index] - reference_value;
    for (Complex &value : field)
    {
        value -= offset;
    }
}

ScalarComplexHelmholtzSolverState solveDiscretePhiForSourceFreeCase(const MatrixFreePairwiseGraph &graph,
                                                                    const StdVec<Real> &sigma,
                                                                    const StdVec<Complex> &rhs,
                                                                    StdVec<Complex> &phi)
{
    ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(phi.size());
    ScalarComplexHelmholtzSolverParameters parameters;
    parameters.max_iterations_ = 8000;
    parameters.relaxation_factor_ = 1.0;
    parameters.absolute_tolerance_ = 1.0e-6;
    parameters.diagonal_regularization_ = 1.0e-12;

    const StdVec<Complex> zero_reaction(phi.size(), Complex(0.0, 0.0));
    return solveScalarComplexHelmholtz(
        phi, rhs, zero_reaction, residuals, parameters,
        [&](const StdVec<Complex> &current_field, ScalarComplexHelmholtzResiduals &current_residuals)
        {
            accumulateScalarLaplaceResidualsFromGraph(graph, current_field, sigma, current_residuals);
        });
}

Real computeResidualL2FromDifference(const StdVec<Complex> &lhs, const StdVec<Complex> &rhs)
{
    if (lhs.empty())
    {
        return 0.0;
    }

    Real squared_sum = 0.0;
    for (size_t i = 0; i != lhs.size(); ++i)
    {
        squared_sum += std::norm(lhs[i] - rhs[i]);
    }
    return std::sqrt(squared_sum / static_cast<Real>(lhs.size()));
}

StdVec<Complex> solveReferenceGaugeChi(const MatrixFreePairwiseGraph &graph, const StdVec<Complex> &rhs,
                                       Real &reference_residual_l2)
{
    using DenseComplexMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
    using DenseComplexVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

    const size_t number_of_particles = rhs.size();
    DenseComplexMatrix matrix(number_of_particles, number_of_particles);
    DenseComplexVector rhs_vector(number_of_particles);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        rhs_vector(static_cast<Eigen::Index>(i)) = rhs[i];
    }

    for (size_t j = 0; j != number_of_particles; ++j)
    {
        StdVec<Complex> basis(number_of_particles, Complex(0.0, 0.0));
        basis[j] = Complex(1.0, 0.0);
        const StdVec<Complex> column = applyScalarDivergenceOfGradientFromGraph(graph, basis);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            matrix(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = column[i];
        }
    }

    DenseComplexVector solution = matrix.completeOrthogonalDecomposition().solve(rhs_vector);
    StdVec<Complex> chi(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        chi[i] = solution(static_cast<Eigen::Index>(i));
    }
    removeMeanOffset(chi);

    const StdVec<Complex> operator_response = applyScalarDivergenceOfGradientFromGraph(graph, chi);
    reference_residual_l2 = computeResidualL2FromDifference(operator_response, rhs);
    return chi;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = readEnvReal("EM_APHI_STAGGERED_SMOKE_DP", 0.05);
    const Real body_length = 1.0;
    const std::string case_mode = readEnvString("EM_APHI_STAGGERED_SMOKE_CASE", "coulomb_variable_sigma_forced_response");

    const bool use_two_dimensional_band =
        case_mode == "coulomb_variable_sigma_source_free" || case_mode == "coulomb_variable_sigma_driven" ||
        case_mode == "coulomb_variable_sigma_forced_response";
    const Real body_height = use_two_dimensional_band ? 0.30 : dp_0;
    const Vecd halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * dp_0);
    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * dp_0);
    BoundingBoxd system_bounds(Vecd(-3.0 * dp_0, -3.0 * dp_0, -3.0 * dp_0),
                               Vecd(body_length + 3.0 * dp_0, body_height + 3.0 * dp_0, 4.0 * dp_0));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system, makeShared<SmokeBoxShape>("SmokeBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = body.getBaseParticles();
    const size_t number_of_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");

    MatrixFreeAPhiParameters parameters;
    parameters.angular_frequency = readEnvReal("EM_APHI_STAGGERED_SMOKE_OMEGA", 5.0);
    parameters.phi_reference_index = 0;
    parameters.phi_reference_value = Complex(0.0, 0.0);

    MatrixFreeAPhiDiscreteView discrete_view{
        number_of_particles,
        body.getSPHAdaptation().ReferenceSmoothingLength(),
        positions,
        volumetric_measure,
        nullptr,
        &inner_relation.inner_configuration_,
        nullptr};

    const MatrixFreePairwiseGraph graph = buildMatrixFreePairwiseGraph(discrete_view, parameters);

    StdVec<Real> sigma(number_of_particles, 1.0);
    StdVec<Real> nu(number_of_particles, 1.0);
    if (case_mode == "coulomb_variable_sigma_source_free" || case_mode == "coulomb_variable_sigma_driven" ||
        case_mode == "coulomb_variable_sigma_forced_response")
    {
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Real y = positions[i][1];
            const Real phase_y = Pi * y / body_height;
            sigma[i] = 1.0 + 0.35 * std::sin(phase_y);
        }
    }

    MatrixFreeAPhiFields exact_fields;
    exact_fields.ax.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.ay.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.az.resize(number_of_particles, Complex(0.0, 0.0));
    exact_fields.phi.resize(number_of_particles, Complex(0.0, 0.0));

    const Complex transverse_amplitude(1.0, 0.1);
    const Complex coupled_ax_amplitude(0.25, -0.08);
    const Complex driven_phi_amplitude(0.12, -0.03);
    const bool has_discrete_manufactured_reference = case_mode != "coulomb_variable_sigma_forced_response";
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const Real x = positions[i][0];
        const Real phase_x = Pi * x / body_length;
        if (case_mode == "coupled_source_free")
        {
            exact_fields.ax[i] = coupled_ax_amplitude * std::sin(phase_x);
        }
        else if (case_mode != "coulomb_variable_sigma_forced_response")
        {
            exact_fields.ay[i] = transverse_amplitude * std::sin(phase_x);
            if (case_mode == "coulomb_variable_sigma_driven")
            {
                const Real y = positions[i][1];
                const Real phase_y = Pi * y / body_height;
                exact_fields.phi[i] = driven_phi_amplitude * std::cos(phase_x) * std::sin(phase_y);
            }
        }
    }
    enforceReference(exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);

    Real reference_phi_build_residual = 0.0;
    if (case_mode == "coupled_source_free" || case_mode == "coulomb_variable_sigma_source_free")
    {
        StdVec<Complex> sigma_ax(number_of_particles, Complex(0.0, 0.0));
        StdVec<Complex> sigma_ay(number_of_particles, Complex(0.0, 0.0));
        StdVec<Complex> sigma_az(number_of_particles, Complex(0.0, 0.0));
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            sigma_ax[i] = sigma[i] * exact_fields.ax[i];
            sigma_ay[i] = sigma[i] * exact_fields.ay[i];
            sigma_az[i] = sigma[i] * exact_fields.az[i];
        }
        const StdVec<Vec3c> grad_sigma_ax = applyMatrixFreeGradient(graph, sigma_ax);
        const StdVec<Vec3c> grad_sigma_ay = applyMatrixFreeGradient(graph, sigma_ay);
        const StdVec<Vec3c> grad_sigma_az = applyMatrixFreeGradient(graph, sigma_az);
        StdVec<Complex> rhs_phi(number_of_particles, Complex(0.0, 0.0));
        const Complex imaginary_unit(0.0, 1.0);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Complex divergence_sigma_a = grad_sigma_ax[i][0] + grad_sigma_ay[i][1] + grad_sigma_az[i][2];
            rhs_phi[i] = imaginary_unit * parameters.angular_frequency * divergence_sigma_a;
        }
        ScalarComplexHelmholtzSolverState exact_phi_state =
            solveDiscretePhiForSourceFreeCase(graph, sigma, rhs_phi, exact_fields.phi);
        enforceReference(exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);
        reference_phi_build_residual = exact_phi_state.current_residual_l2_;
        const Real acceptable_reference_residual = 1.0e-5;
        if (!exact_phi_state.converged_ && exact_phi_state.current_residual_l2_ > acceptable_reference_residual)
        {
            std::cerr << "failed_to_build_coupled_source_free_phi_reference"
                      << " iterations=" << exact_phi_state.iterations_
                      << " residual_l2=" << exact_phi_state.current_residual_l2_ << std::endl;
            return 1;
        }
    }

    const StdVec<Complex> laplace_ax = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.ax, nu);
    const StdVec<Complex> laplace_ay = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.ay, nu);
    const StdVec<Complex> laplace_az = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.az, nu);
    const StdVec<Complex> laplace_phi_exact = applyScalarNegativeLaplaceFromGraph(graph, exact_fields.phi, sigma);
    const StdVec<Vec3c> grad_phi = applyMatrixFreeGradient(graph, exact_fields.phi);

    MatrixFreeAPhiSources sources;
    sources.source_ax.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_ay.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_az.resize(number_of_particles, Complex(0.0, 0.0));
    sources.source_phi.resize(number_of_particles, Complex(0.0, 0.0));

    const Complex imaginary_unit(0.0, 1.0);
    StdVec<Complex> sigma_ax(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> sigma_ay(number_of_particles, Complex(0.0, 0.0));
    StdVec<Complex> sigma_az(number_of_particles, Complex(0.0, 0.0));
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        sigma_ax[i] = sigma[i] * exact_fields.ax[i];
        sigma_ay[i] = sigma[i] * exact_fields.ay[i];
        sigma_az[i] = sigma[i] * exact_fields.az[i];
    }
    const StdVec<Vec3c> grad_sigma_ax_for_source = applyMatrixFreeGradient(graph, sigma_ax);
    const StdVec<Vec3c> grad_sigma_ay_for_source = applyMatrixFreeGradient(graph, sigma_ay);
    const StdVec<Vec3c> grad_sigma_az_for_source = applyMatrixFreeGradient(graph, sigma_az);
    const Complex forced_source_amplitude(0.8, 0.15);
    const std::string forced_profile = readEnvString("EM_APHI_FORCED_PROFILE", "sin");
    const Real forced_gaussian_center = readEnvReal("EM_APHI_FORCED_CENTER_X", 0.5 * body_length);
    const Real forced_gaussian_width = readEnvReal("EM_APHI_FORCED_WIDTH", 0.15 * body_length);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        sources.source_ax[i] = laplace_ax[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.ax[i] +
                               sigma[i] * grad_phi[i][0];
        sources.source_ay[i] = laplace_ay[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.ay[i] +
                               sigma[i] * grad_phi[i][1];
        sources.source_az[i] = laplace_az[i] + imaginary_unit * parameters.angular_frequency * sigma[i] * exact_fields.az[i] +
                               sigma[i] * grad_phi[i][2];
        if (case_mode == "coulomb_variable_sigma_driven")
        {
            const Complex divergence_sigma_a = grad_sigma_ax_for_source[i][0] + grad_sigma_ay_for_source[i][1] +
                                               grad_sigma_az_for_source[i][2];
            sources.source_phi[i] = laplace_phi_exact[i] - imaginary_unit * parameters.angular_frequency * divergence_sigma_a;
        }
        else
        {
            sources.source_phi[i] = Complex(0.0, 0.0);
        }
    }

    if (case_mode == "coulomb_variable_sigma_forced_response")
    {
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Real x = positions[i][0];
            const Real phase_x = Pi * x / body_length;
            Real profile_weight = std::sin(phase_x);
            if (forced_profile == "gaussian")
            {
                const Real width = SMAX(forced_gaussian_width, static_cast<Real>(0.1) * dp_0);
                const Real normalized_distance = (x - forced_gaussian_center) / width;
                profile_weight = std::exp(-normalized_distance * normalized_distance);
            }
            sources.source_ax[i] = Complex(0.0, 0.0);
            sources.source_ay[i] = forced_source_amplitude * profile_weight;
            sources.source_az[i] = Complex(0.0, 0.0);
            sources.source_phi[i] = Complex(0.0, 0.0);
        }
    }

    MatrixFreeAPhiSolverParameters solver_parameters;
    const size_t default_outer_iterations =
        case_mode == "coupled_source_free" || case_mode == "coulomb_variable_sigma_source_free" ||
                case_mode == "coulomb_variable_sigma_driven" ||
                case_mode == "coulomb_variable_sigma_forced_response"
            ? 160
            : 80;
    solver_parameters.max_outer_iterations =
        readEnvSizeT("EM_APHI_STAGGERED_SMOKE_OUTER_ITERS", default_outer_iterations);
    solver_parameters.a_component_solver.max_iterations_ = readEnvSizeT("EM_APHI_STAGGERED_SMOKE_A_ITERS", 800);
    solver_parameters.a_component_solver.absolute_tolerance_ =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_A_TOL", 1.0e-6);
    solver_parameters.phi_solver.max_iterations_ = readEnvSizeT("EM_APHI_STAGGERED_SMOKE_PHI_ITERS", 800);
    solver_parameters.phi_solver.absolute_tolerance_ = readEnvReal("EM_APHI_STAGGERED_SMOKE_PHI_TOL", 1.0e-6);
    solver_parameters.enable_gauge_projection = readEnvBool("EM_APHI_STAGGERED_SMOKE_ENABLE_GAUGE", false);
    const std::string gauge_operator_mode = readEnvString("EM_APHI_STAGGERED_SMOKE_GAUGE_OPERATOR", "consistent");
    solver_parameters.use_operator_consistent_gauge_projection = gauge_operator_mode != "laplace";
    if (case_mode == "coupled_source_free" && solver_parameters.enable_gauge_projection)
    {
        std::cerr << "coupled_source_free_does_not_match_coulomb_gauge_in_quasi_1d" << std::endl;
        return 2;
    }
    solver_parameters.residual_tolerance = readEnvReal("EM_APHI_STAGGERED_SMOKE_RESIDUAL_TOL", 1.0e-5);
    solver_parameters.divergence_tolerance = readEnvReal("EM_APHI_STAGGERED_SMOKE_DIVERGENCE_TOL", 1.0e-5);
    solver_parameters.enable_gauge_penalty = readEnvBool("EM_APHI_ENABLE_GAUGE_PENALTY", false);
    solver_parameters.gauge_penalty_coefficient =
        readEnvReal("EM_APHI_GAUGE_PENALTY_COEFF", solver_parameters.enable_gauge_penalty ? 1.0 : 0.0);
    solver_parameters.gauge_penalty_ramp_iterations =
        readEnvSizeT("EM_APHI_GAUGE_PENALTY_RAMP_ITERS", 0);
    solver_parameters.gauge_penalty_initial_ratio =
        readEnvReal("EM_APHI_GAUGE_PENALTY_INITIAL_RATIO", 1.0);
    solver_parameters.outer_relaxation_factor =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_OUTER_RELAX", solver_parameters.enable_gauge_penalty ? 0.8 : 1.0);
    solver_parameters.update_tolerance =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_UPDATE_TOL", 1.0e-7);

    if (solver_parameters.enable_gauge_penalty)
    {
        MatrixFreeAPhiFields exact_penalty_fields = exact_fields;
        const StdVec<Complex> exact_divergence_a = computeDivergenceOfVectorFieldLocal(graph, exact_penalty_fields);
        const StdVec<Vec3c> exact_gauge_penalty_gradient = applyMatrixFreeGradient(graph, exact_divergence_a);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            sources.source_ax[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][0];
            sources.source_ay[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][1];
            sources.source_az[i] += solver_parameters.gauge_penalty_coefficient * exact_gauge_penalty_gradient[i][2];
        }
    }

    if (solver_parameters.enable_gauge_projection)
    {
        solver_parameters.gauge_solver.max_iterations_ =
            readEnvSizeT("EM_APHI_STAGGERED_SMOKE_GAUGE_ITERS",
                         solver_parameters.use_operator_consistent_gauge_projection ? 12000 : 500);
        solver_parameters.gauge_solver.relaxation_factor_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_RELAX",
                        solver_parameters.use_operator_consistent_gauge_projection ? 5.0e-4 : 1.0);
        solver_parameters.gauge_solver.diagonal_regularization_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_DIAG_REG",
                        solver_parameters.use_operator_consistent_gauge_projection ? 1.0e-6 : 1.0e-12);
        solver_parameters.gauge_solver.absolute_tolerance_ =
            readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_TOL", 1.0e-6);
    }

    const Real source_ax_l2 = computeComplexFieldL2NormLocal(sources.source_ax);
    const Real source_ay_l2 = computeComplexFieldL2NormLocal(sources.source_ay);
    const Real source_phi_l2 = computeComplexFieldL2NormLocal(sources.source_phi);

    const bool gauge_only_analysis = readEnvBool("EM_APHI_STAGGERED_SMOKE_GAUGE_ONLY", false);
    if (case_mode != "coulomb_variable_sigma_forced_response")
    {
        std::cerr << "physical_response_target_only_supports_coulomb_variable_sigma_forced_response" << std::endl;
        return 3;
    }

    if (gauge_only_analysis && !solver_parameters.enable_gauge_projection)
    {
        std::cerr << "gauge_only_requires_enable_gauge" << std::endl;
        return 4;
    }

    const StdVec<Complex> exact_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, exact_fields);
    const Real exact_divergence_a_l2 =
        has_discrete_manufactured_reference ? computeComplexFieldL2NormLocal(exact_divergence_a_field) : Real(-1.0);

    Real gauge_seed_div_delta_l2 = 0.0;
    Real gauge_seed_operator_response_l2 = 0.0;
    Real gauge_seed_match_l2 = 0.0;
    Real gauge_seed_match_neg_l2 = 0.0;

    MatrixFreeAPhiFields fields;
    fields.ax.resize(number_of_particles, Complex(0.0, 0.0));
    fields.ay.resize(number_of_particles, Complex(0.0, 0.0));
    fields.az.resize(number_of_particles, Complex(0.0, 0.0));
    fields.phi.resize(number_of_particles, Complex(0.0, 0.0));

    if (solver_parameters.enable_gauge_projection)
    {
        const Real seed_amplitude = readEnvReal("EM_APHI_STAGGERED_SMOKE_GAUGE_SEED", 0.2);
        StdVec<Complex> chi_seed(number_of_particles, Complex(0.0, 0.0));
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Real x = positions[i][0];
            chi_seed[i] = seed_amplitude * std::cos(Pi * x / body_length);
        }
        const StdVec<Vec3c> grad_chi_seed = applyMatrixFreeGradient(graph, chi_seed);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            fields.ax[i] = exact_fields.ax[i] + grad_chi_seed[i][0];
            fields.ay[i] = exact_fields.ay[i] + grad_chi_seed[i][1];
            fields.az[i] = exact_fields.az[i] + grad_chi_seed[i][2];
            fields.phi[i] = exact_fields.phi[i] - imaginary_unit * parameters.angular_frequency * chi_seed[i];
        }

        const StdVec<Complex> seeded_divergence_a = computeDivergenceOfVectorFieldLocal(graph, fields);
        const StdVec<Complex> operator_response = applyScalarDivergenceOfGradientFromGraph(graph, chi_seed);
        StdVec<Complex> seeded_divergence_delta(number_of_particles, Complex(0.0, 0.0));
        StdVec<Complex> operator_response_negative(number_of_particles, Complex(0.0, 0.0));
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            seeded_divergence_delta[i] = seeded_divergence_a[i] - exact_divergence_a_field[i];
            operator_response_negative[i] = -operator_response[i];
        }
        gauge_seed_div_delta_l2 = computeComplexFieldL2NormLocal(seeded_divergence_delta);
        gauge_seed_operator_response_l2 = computeComplexFieldL2NormLocal(operator_response);
        gauge_seed_match_l2 = computeL2Error(seeded_divergence_delta, operator_response);
        gauge_seed_match_neg_l2 = computeL2Error(seeded_divergence_delta, operator_response_negative);
    }

    if (gauge_only_analysis)
    {
        MatrixFreeAPhiFields pregauge_fields = fields;
        MatrixFreeAPhiSolverParameters pregauge_solver = solver_parameters;
        pregauge_solver.enable_gauge_projection = false;
        pregauge_solver.max_outer_iterations =
            readEnvSizeT("EM_APHI_STAGGERED_SMOKE_PREGAUGE_OUTER_ITERS", default_outer_iterations);

        const MatrixFreeAPhiSolverState pregauge_state =
            solveMatrixFreeAPhiStaggered(graph, pregauge_fields, sigma, nu, parameters, sources, pregauge_solver);
        const StdVec<Complex> pregauge_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        const Real pregauge_divergence_a_error_l2 = computeL2Error(pregauge_divergence_a_field, exact_divergence_a_field);
        const MatrixFreeAPhiGaugeProjectionResult isolated_gauge =
            applyMatrixFreeAPhiGaugeProjectionStep(graph, pregauge_fields, sigma, parameters,
                                                  solver_parameters.gauge_solver,
                                                  solver_parameters.remove_gauge_mean_offset, true,
                                                  solver_parameters.use_operator_consistent_gauge_projection);
        const MatrixFreeAPhiResiduals post_gauge_residuals =
            evaluateMatrixFreeAPhiResiduals(graph, pregauge_fields, sigma, nu, parameters, sources,
                                          solver_parameters.enable_gauge_penalty,
                                          solver_parameters.gauge_penalty_coefficient);
        const StdVec<Complex> post_gauge_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        const Real post_gauge_divergence_a_error_l2 = computeL2Error(post_gauge_divergence_a_field, exact_divergence_a_field);

        StdVec<Complex> reference_rhs = computeDivergenceOfVectorFieldLocal(graph, pregauge_fields);
        removeMeanOffset(reference_rhs);
        Real reference_gauge_solver_residual_l2 = 0.0;
        const StdVec<Complex> reference_chi = solveReferenceGaugeChi(graph, reference_rhs, reference_gauge_solver_residual_l2);
        MatrixFreeAPhiFields reference_gauge_fields = pregauge_fields;
        const StdVec<Vec3c> reference_grad_chi = applyMatrixFreeGradient(graph, reference_chi);
        for (size_t i = 0; i != reference_gauge_fields.ax.size(); ++i)
        {
            reference_gauge_fields.ax[i] -= reference_grad_chi[i][0];
            reference_gauge_fields.ay[i] -= reference_grad_chi[i][1];
            reference_gauge_fields.az[i] -= reference_grad_chi[i][2];
            reference_gauge_fields.phi[i] += imaginary_unit * parameters.angular_frequency * reference_chi[i];
        }
        enforceReference(reference_gauge_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);
        const MatrixFreeAPhiResiduals reference_post_gauge_residuals =
            evaluateMatrixFreeAPhiResiduals(graph, reference_gauge_fields, sigma, nu, parameters, sources,
                                          solver_parameters.enable_gauge_penalty,
                                          solver_parameters.gauge_penalty_coefficient);
        const StdVec<Complex> reference_divergence_a = computeDivergenceOfVectorFieldLocal(graph, reference_gauge_fields);
        const Real reference_gauge_div_a_after_l2 = computeComplexFieldL2NormLocal(reference_divergence_a);
        const Real reference_gauge_div_a_error_l2 = computeL2Error(reference_divergence_a, exact_divergence_a_field);
        const Real reference_gauge_chi_l2 = computeComplexFieldL2NormLocal(reference_chi);

        const Real ax_l2_error = computeL2Error(pregauge_fields.ax, exact_fields.ax);
        const Real ay_l2_error = computeL2Error(pregauge_fields.ay, exact_fields.ay);
        const Real phi_l2_error = computeL2Error(pregauge_fields.phi, exact_fields.phi);

        std::cout << std::setprecision(12)
                  << "matrix_free_aphi_gauge_isolation_smoke"
                  << " case_mode=" << case_mode
                  << " dp=" << dp_0
                  << " points=" << number_of_particles
                  << " gauge_operator_mode=" << gauge_operator_mode
                  << " gauge_penalty_enabled=" << (solver_parameters.enable_gauge_penalty ? 1 : 0)
                  << " gauge_penalty_coeff=" << solver_parameters.gauge_penalty_coefficient
                  << " effective_gauge_penalty_coeff=" << pregauge_state.effective_gauge_penalty_coefficient
                  << " outer_relaxation_factor=" << solver_parameters.outer_relaxation_factor
                  << " update_tolerance=" << solver_parameters.update_tolerance
                  << " gauge_seed_div_delta_l2=" << gauge_seed_div_delta_l2
                  << " gauge_seed_operator_response_l2=" << gauge_seed_operator_response_l2
                  << " gauge_seed_match_l2=" << gauge_seed_match_l2
                  << " gauge_seed_match_neg_l2=" << gauge_seed_match_neg_l2
                  << " pregauge_outer_iterations=" << pregauge_state.outer_iterations
                  << " pregauge_converged=" << (pregauge_state.converged ? 1 : 0)
                  << " pregauge_residual_ax_l2=" << pregauge_state.residuals.residual_ax_l2
                  << " pregauge_residual_ay_l2=" << pregauge_state.residuals.residual_ay_l2
                  << " pregauge_residual_phi_l2=" << pregauge_state.residuals.residual_phi_l2
                  << " pregauge_divergence_a_l2=" << pregauge_state.residuals.divergence_a_l2
                  << " pregauge_divergence_j_l2=" << pregauge_state.residuals.divergence_j_l2
                  << " exact_divergence_a_l2=" << exact_divergence_a_l2
                  << " pregauge_divergence_a_error_l2=" << pregauge_divergence_a_error_l2
                  << " gauge_solver_iterations=" << isolated_gauge.solver_state.iterations_
                  << " gauge_solver_converged=" << (isolated_gauge.solver_state.converged_ ? 1 : 0)
                  << " gauge_solver_residual_l2=" << isolated_gauge.solver_state.current_residual_l2_
                  << " gauge_solver_min_diagonal_abs=" << isolated_gauge.solver_state.current_min_diagonal_abs_
                  << " gauge_solver_max_diagonal_abs=" << isolated_gauge.solver_state.current_max_diagonal_abs_
                  << " gauge_solver_nonfinite_diagonal_count=" << isolated_gauge.solver_state.current_nonfinite_diagonal_count_
                  << " gauge_rhs_mean_abs_before_projection=" << isolated_gauge.diagnostics.gauge_rhs_mean_abs_before_projection
                  << " gauge_rhs_mean_abs_after_projection=" << isolated_gauge.diagnostics.gauge_rhs_mean_abs_after_projection
                  << " gauge_chi_mean_abs_after_solve=" << isolated_gauge.diagnostics.chi_mean_abs_after_solve
                  << " gauge_chi_l2=" << isolated_gauge.diagnostics.chi_l2
                  << " gauge_grad_chi_l2=" << isolated_gauge.diagnostics.grad_chi_l2
                  << " gauge_chi_max_abs=" << isolated_gauge.diagnostics.chi_max_abs
                  << " gauge_phi_ref_after_phi_abs=" << pregauge_state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs
                  << " gauge_phi_ref_after_update_abs=" << isolated_gauge.diagnostics.phi_reference_offset_after_gauge_update_abs
                  << " gauge_phi_ref_after_final_abs=" << isolated_gauge.diagnostics.phi_reference_offset_after_final_reference_abs
                  << " gauge_div_a_before_l2=" << isolated_gauge.diagnostics.divergence_a_before_l2
                  << " gauge_div_a_after_raw_l2=" << isolated_gauge.diagnostics.divergence_a_after_raw_l2
                  << " gauge_div_a_after_final_l2=" << isolated_gauge.diagnostics.divergence_a_after_final_l2
                  << " gauge_div_j_before_l2=" << isolated_gauge.diagnostics.divergence_j_before_l2
                  << " gauge_div_j_after_raw_l2=" << isolated_gauge.diagnostics.divergence_j_after_raw_l2
                  << " gauge_div_j_after_final_l2=" << isolated_gauge.diagnostics.divergence_j_after_final_l2
                  << " gauge_e_before_l2=" << isolated_gauge.diagnostics.electric_field_before_l2
                  << " gauge_e_after_raw_l2=" << isolated_gauge.diagnostics.electric_field_after_raw_l2
                  << " gauge_e_after_final_l2=" << isolated_gauge.diagnostics.electric_field_after_final_l2
                  << " gauge_e_change_raw_l2=" << isolated_gauge.diagnostics.electric_field_change_raw_l2
                  << " gauge_e_change_final_l2=" << isolated_gauge.diagnostics.electric_field_change_final_l2
                  << " gauge_j_before_l2=" << isolated_gauge.diagnostics.current_density_before_l2
                  << " gauge_j_after_raw_l2=" << isolated_gauge.diagnostics.current_density_after_raw_l2
                  << " gauge_j_after_final_l2=" << isolated_gauge.diagnostics.current_density_after_final_l2
                  << " gauge_j_change_raw_l2=" << isolated_gauge.diagnostics.current_density_change_raw_l2
                  << " gauge_j_change_final_l2=" << isolated_gauge.diagnostics.current_density_change_final_l2
                  << " post_gauge_residual_ax_l2=" << post_gauge_residuals.residual_ax_l2
                  << " post_gauge_residual_ay_l2=" << post_gauge_residuals.residual_ay_l2
                  << " post_gauge_residual_phi_l2=" << post_gauge_residuals.residual_phi_l2
                  << " post_gauge_divergence_a_l2=" << post_gauge_residuals.divergence_a_l2
                  << " post_gauge_divergence_j_l2=" << post_gauge_residuals.divergence_j_l2
                  << " post_gauge_divergence_a_error_l2=" << post_gauge_divergence_a_error_l2
                  << " reference_gauge_solver_residual_l2=" << reference_gauge_solver_residual_l2
                  << " reference_gauge_chi_l2=" << reference_gauge_chi_l2
                  << " reference_gauge_div_a_after_l2=" << reference_gauge_div_a_after_l2
                  << " reference_gauge_div_a_error_l2=" << reference_gauge_div_a_error_l2
                  << " reference_post_gauge_residual_ax_l2=" << reference_post_gauge_residuals.residual_ax_l2
                  << " reference_post_gauge_residual_phi_l2=" << reference_post_gauge_residuals.residual_phi_l2
                  << " ax_l2_error=" << ax_l2_error
                  << " ay_l2_error=" << ay_l2_error
                  << " phi_l2_error=" << phi_l2_error
                  << std::endl;

        return pregauge_state.outer_iterations == 0 ? 1 : 0;
    }

    const MatrixFreeAPhiSolverState state =
        solveMatrixFreeAPhiStaggered(graph, fields, sigma, nu, parameters, sources, solver_parameters);

    const StdVec<Complex> final_divergence_a_field = computeDivergenceOfVectorFieldLocal(graph, fields);
    const Real final_divergence_a_error_l2 =
        has_discrete_manufactured_reference ? computeL2Error(final_divergence_a_field, exact_divergence_a_field) : Real(-1.0);
    const Real ax_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.ax, exact_fields.ax) : Real(-1.0);
    const Real ay_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.ay, exact_fields.ay) : Real(-1.0);
    const Real phi_l2_error =
        has_discrete_manufactured_reference ? computeL2Error(fields.phi, exact_fields.phi) : Real(-1.0);
    const Real ax_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.ax, exact_fields.ax) : Real(-1.0);
    const Real ay_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.ay, exact_fields.ay) : Real(-1.0);
    const Real phi_max_error =
        has_discrete_manufactured_reference ? computeMaxError(fields.phi, exact_fields.phi) : Real(-1.0);

    const StdVec<Vec3c> exact_electric_field =
        has_discrete_manufactured_reference ? computeElectricField(graph, exact_fields, parameters) : StdVec<Vec3c>();
    const StdVec<Vec3c> exact_current_density =
        has_discrete_manufactured_reference ? computeCurrentDensity(graph, exact_fields, sigma, parameters) : StdVec<Vec3c>();
    const StdVec<Real> exact_joule_density =
        has_discrete_manufactured_reference ? computeJouleHeatingDensity(exact_electric_field, exact_current_density) : StdVec<Real>();
    const StdVec<Vec3c> final_electric_field = computeElectricField(graph, fields, parameters);
    const StdVec<Vec3c> final_current_density = computeCurrentDensity(graph, fields, sigma, parameters);
    const StdVec<Real> final_joule_density = computeJouleHeatingDensity(final_electric_field, final_current_density);
    const Real electric_l2_error =
        has_discrete_manufactured_reference ? computeVectorFieldL2Error(final_electric_field, exact_electric_field) : Real(-1.0);
    const Real current_l2_error =
        has_discrete_manufactured_reference ? computeVectorFieldL2Error(final_current_density, exact_current_density) : Real(-1.0);
    const Real joule_l2_error =
        has_discrete_manufactured_reference ? computeRealFieldL2Error(final_joule_density, exact_joule_density) : Real(-1.0);
    const Real electric_field_l2 = computeVectorComplexFieldL2Norm(final_electric_field);
    const Real current_density_l2 = computeVectorComplexFieldL2Norm(final_current_density);
    const Real joule_density_l2 = computeRealFieldL2Norm(final_joule_density);
    const Real joule_density_max = final_joule_density.empty() ? 0.0 : *std::max_element(final_joule_density.begin(), final_joule_density.end());
    const Real symmetry_tolerance = 0.25 * dp_0;
    const Real ay_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(fields.ay, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(fields.phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_mirror_antisymmetry_l2 =
        computeComplexFieldMirrorAntisymmetryError(fields.phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const StdVec<Complex> centered_phi = makeMeanCenteredCopy(fields.phi);
    const Real phi_centered_mirror_symmetry_l2 =
        computeComplexFieldMirrorSymmetryError(centered_phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real phi_centered_mirror_antisymmetry_l2 =
        computeComplexFieldMirrorAntisymmetryError(centered_phi, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);
    const Real joule_mirror_symmetry_l2 =
        computeRealFieldMirrorSymmetryError(final_joule_density, positions, number_of_particles, 0.5 * body_height, symmetry_tolerance);

    const bool require_converged = readEnvBool("EM_APHI_STAGGERED_SMOKE_REQUIRE_CONVERGED", false);
    const Real validation_max_residual_ay = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_RESIDUAL_AY", -1.0);
    const Real validation_max_divergence_a_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_DIVA_ERROR", -1.0);
    const Real validation_max_ay_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_AY_ERROR", -1.0);
    const Real validation_max_phi_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_ERROR", -1.0);
    const Real validation_max_e_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_E_ERROR", -1.0);
    const Real validation_max_j_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_J_ERROR", -1.0);
    const Real validation_max_joule_error = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_JOULE_ERROR", -1.0);
    const Real validation_max_ay_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_AY_MIRROR", -1.0);
    const Real validation_max_phi_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_MIRROR", -1.0);
    const Real validation_max_joule_mirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_JOULE_MIRROR", -1.0);
    const Real validation_max_phi_antimirror = readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_ANTIMIRROR", -1.0);
    const Real validation_max_phi_centered_mirror =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_MIRROR", -1.0);
    const Real validation_max_phi_centered_antimirror =
        readEnvReal("EM_APHI_STAGGERED_SMOKE_MAX_PHI_CENTERED_ANTIMIRROR", -1.0);

    bool validation_pass = true;
    if (require_converged && !state.converged)
    {
        validation_pass = false;
    }
    if (validation_max_residual_ay >= Real(0.0) && state.residuals.residual_ay_l2 > validation_max_residual_ay)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_divergence_a_error >= Real(0.0) &&
        final_divergence_a_error_l2 > validation_max_divergence_a_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_ay_error >= Real(0.0) &&
        ay_l2_error > validation_max_ay_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_phi_error >= Real(0.0) &&
        phi_l2_error > validation_max_phi_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_e_error >= Real(0.0) &&
        electric_l2_error > validation_max_e_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_j_error >= Real(0.0) &&
        current_l2_error > validation_max_j_error)
    {
        validation_pass = false;
    }
    if (has_discrete_manufactured_reference && validation_max_joule_error >= Real(0.0) &&
        joule_l2_error > validation_max_joule_error)
    {
        validation_pass = false;
    }
    if (validation_max_ay_mirror >= Real(0.0) && ay_mirror_symmetry_l2 > validation_max_ay_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_mirror >= Real(0.0) && phi_mirror_symmetry_l2 > validation_max_phi_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_joule_mirror >= Real(0.0) && joule_mirror_symmetry_l2 > validation_max_joule_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_antimirror >= Real(0.0) &&
        phi_mirror_antisymmetry_l2 > validation_max_phi_antimirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_centered_mirror >= Real(0.0) &&
        phi_centered_mirror_symmetry_l2 > validation_max_phi_centered_mirror)
    {
        validation_pass = false;
    }
    if (validation_max_phi_centered_antimirror >= Real(0.0) &&
        phi_centered_mirror_antisymmetry_l2 > validation_max_phi_centered_antimirror)
    {
        validation_pass = false;
    }

    std::cout << std::setprecision(12)
              << "matrix_free_aphi_physical_response"
              << " case_mode=" << case_mode
              << " dp=" << dp_0
              << " points=" << number_of_particles
              << " has_discrete_reference=" << (has_discrete_manufactured_reference ? 1 : 0)
              << " gauge_operator_mode=" << gauge_operator_mode
              << " gauge_enabled=" << (solver_parameters.enable_gauge_projection ? 1 : 0)
              << " gauge_penalty_enabled=" << (solver_parameters.enable_gauge_penalty ? 1 : 0)
              << " gauge_penalty_coeff=" << solver_parameters.gauge_penalty_coefficient
              << " effective_gauge_penalty_coeff=" << state.effective_gauge_penalty_coefficient
              << " outer_relaxation_factor=" << solver_parameters.outer_relaxation_factor
              << " update_tolerance=" << solver_parameters.update_tolerance
              << " field_update_l2=" << state.field_update_l2
              << " relative_field_update_l2=" << state.relative_field_update_l2
              << " outer_iterations=" << state.outer_iterations
              << " converged=" << (state.converged ? 1 : 0)
              << " residual_ax_l2=" << state.residuals.residual_ax_l2
              << " residual_ay_l2=" << state.residuals.residual_ay_l2
              << " residual_phi_l2=" << state.residuals.residual_phi_l2
              << " divergence_a_l2=" << state.residuals.divergence_a_l2
              << " divergence_j_l2=" << state.residuals.divergence_j_l2
              << " exact_divergence_a_l2=" << exact_divergence_a_l2
              << " divergence_a_error_l2=" << final_divergence_a_error_l2
              << " ax_l2_error=" << ax_l2_error
              << " ay_l2_error=" << ay_l2_error
              << " phi_l2_error=" << phi_l2_error
              << " electric_l2_error=" << electric_l2_error
              << " current_l2_error=" << current_l2_error
              << " joule_l2_error=" << joule_l2_error
              << " electric_field_l2=" << electric_field_l2
              << " current_density_l2=" << current_density_l2
              << " joule_density_l2=" << joule_density_l2
              << " joule_density_max=" << joule_density_max
              << " ay_mirror_symmetry_l2=" << ay_mirror_symmetry_l2
              << " phi_mirror_symmetry_l2=" << phi_mirror_symmetry_l2
              << " phi_mirror_antisymmetry_l2=" << phi_mirror_antisymmetry_l2
              << " phi_centered_mirror_symmetry_l2=" << phi_centered_mirror_symmetry_l2
              << " phi_centered_mirror_antisymmetry_l2=" << phi_centered_mirror_antisymmetry_l2
              << " joule_mirror_symmetry_l2=" << joule_mirror_symmetry_l2
              << " ax_max_error=" << ax_max_error
              << " ay_max_error=" << ay_max_error
              << " phi_max_error=" << phi_max_error
              << " source_ax_l2=" << source_ax_l2
              << " source_ay_l2=" << source_ay_l2
              << " source_phi_l2=" << source_phi_l2
              << " forced_profile=" << forced_profile
              << " forced_center_x=" << forced_gaussian_center
              << " forced_width=" << forced_gaussian_width
              << " reference_phi_build_residual=" << reference_phi_build_residual
              << " gauge_applied=" << (state.gauge_diagnostics.applied ? 1 : 0)
              << " gauge_solver_min_diagonal_abs=" << state.gauge_state.current_min_diagonal_abs_
              << " gauge_solver_max_diagonal_abs=" << state.gauge_state.current_max_diagonal_abs_
              << " gauge_solver_nonfinite_diagonal_count=" << state.gauge_state.current_nonfinite_diagonal_count_
              << " gauge_rhs_mean_abs_before_projection=" << state.gauge_diagnostics.gauge_rhs_mean_abs_before_projection
              << " gauge_rhs_mean_abs_after_projection=" << state.gauge_diagnostics.gauge_rhs_mean_abs_after_projection
              << " gauge_chi_mean_abs_after_solve=" << state.gauge_diagnostics.chi_mean_abs_after_solve
              << " gauge_chi_l2=" << state.gauge_diagnostics.chi_l2
              << " gauge_grad_chi_l2=" << state.gauge_diagnostics.grad_chi_l2
              << " gauge_chi_max_abs=" << state.gauge_diagnostics.chi_max_abs
              << " gauge_phi_ref_after_phi_abs=" << state.gauge_diagnostics.phi_reference_offset_after_phi_solve_abs
              << " gauge_phi_ref_after_update_abs=" << state.gauge_diagnostics.phi_reference_offset_after_gauge_update_abs
              << " gauge_phi_ref_after_final_abs=" << state.gauge_diagnostics.phi_reference_offset_after_final_reference_abs
              << " gauge_div_a_before_l2=" << state.gauge_diagnostics.divergence_a_before_l2
              << " gauge_div_a_after_raw_l2=" << state.gauge_diagnostics.divergence_a_after_raw_l2
              << " gauge_div_a_after_final_l2=" << state.gauge_diagnostics.divergence_a_after_final_l2
              << " gauge_div_j_before_l2=" << state.gauge_diagnostics.divergence_j_before_l2
              << " gauge_div_j_after_raw_l2=" << state.gauge_diagnostics.divergence_j_after_raw_l2
              << " gauge_div_j_after_final_l2=" << state.gauge_diagnostics.divergence_j_after_final_l2
              << " gauge_e_before_l2=" << state.gauge_diagnostics.electric_field_before_l2
              << " gauge_e_after_raw_l2=" << state.gauge_diagnostics.electric_field_after_raw_l2
              << " gauge_e_after_final_l2=" << state.gauge_diagnostics.electric_field_after_final_l2
              << " gauge_e_change_raw_l2=" << state.gauge_diagnostics.electric_field_change_raw_l2
              << " gauge_e_change_final_l2=" << state.gauge_diagnostics.electric_field_change_final_l2
              << " gauge_j_before_l2=" << state.gauge_diagnostics.current_density_before_l2
              << " gauge_j_after_raw_l2=" << state.gauge_diagnostics.current_density_after_raw_l2
              << " gauge_j_after_final_l2=" << state.gauge_diagnostics.current_density_after_final_l2
              << " gauge_j_change_raw_l2=" << state.gauge_diagnostics.current_density_change_raw_l2
              << " gauge_j_change_final_l2=" << state.gauge_diagnostics.current_density_change_final_l2
              << " validation_require_converged=" << (require_converged ? 1 : 0)
              << " validation_max_residual_ay=" << validation_max_residual_ay
              << " validation_max_diva_error=" << validation_max_divergence_a_error
              << " validation_max_ay_error=" << validation_max_ay_error
              << " validation_max_phi_error=" << validation_max_phi_error
              << " validation_max_e_error=" << validation_max_e_error
              << " validation_max_j_error=" << validation_max_j_error
              << " validation_max_joule_error=" << validation_max_joule_error
              << " validation_max_ay_mirror=" << validation_max_ay_mirror
              << " validation_max_phi_mirror=" << validation_max_phi_mirror
              << " validation_max_joule_mirror=" << validation_max_joule_mirror
              << " validation_max_phi_antimirror=" << validation_max_phi_antimirror
              << " validation_max_phi_centered_mirror=" << validation_max_phi_centered_mirror
              << " validation_max_phi_centered_antimirror=" << validation_max_phi_centered_antimirror
              << " validation_pass=" << (validation_pass ? 1 : 0)
              << std::endl;

    if (state.outer_iterations == 0)
    {
        return 1;
    }
    return validation_pass ? 0 : 4;
}
