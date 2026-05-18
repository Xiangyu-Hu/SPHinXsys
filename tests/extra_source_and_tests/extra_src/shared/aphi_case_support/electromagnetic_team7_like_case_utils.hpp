#ifndef ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_HPP
#define ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_HPP

#include "aphi_case_support/electromagnetic_team7_like_case_utils.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_aphi_solver.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>

namespace SPH
{
namespace electromagnetics
{
namespace team7_like
{
namespace
{
Real readEnvReal(const char *name, Real default_value)
{
    const char *value = std::getenv(name);
    if (value == nullptr)
    {
        return default_value;
    }
    char *end_ptr = nullptr;
    const Real parsed = static_cast<Real>(std::strtod(value, &end_ptr));
    return end_ptr == value || !std::isfinite(parsed) ? default_value : parsed;
}
} // namespace

inline bool insideBoxRegion(const Vecd &position, Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax)
{
    return position[0] >= xmin && position[0] <= xmax && position[1] >= ymin && position[1] <= ymax &&
           position[2] >= zmin && position[2] <= zmax;
}

inline Real coordinateFromFraction(Real fraction, Real extent)
{
    return fraction * extent;
}

inline Real smoothBoxProfile(const Vecd &position, Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax,
                             Real amplitude)
{
    if (!insideBoxRegion(position, xmin, xmax, ymin, ymax, zmin, zmax))
    {
        return 0.0;
    }
    const Real xi = (position[0] - xmin) / (xmax - xmin + TinyReal);
    const Real eta = (position[1] - ymin) / (ymax - ymin + TinyReal);
    const Real zeta = (position[2] - zmin) / (zmax - zmin + TinyReal);
    return amplitude * std::sin(Pi * xi) * std::sin(Pi * eta) * std::sin(Pi * zeta);
}

inline bool isExteriorBoundaryParticle(const Vecd &position, const Team7LikeGeometry &geometry, Real shell_thickness)
{
    return position[0] < shell_thickness || position[0] > geometry.body_length - shell_thickness ||
           position[1] < shell_thickness || position[1] > geometry.body_height - shell_thickness ||
           position[2] < shell_thickness || position[2] > geometry.body_width - shell_thickness;
}

inline Team7LikeGeometry readTeam7LikeGeometryFromEnvironment(Real dp_0)
{
    Team7LikeGeometry geometry;
    geometry.body_length = readEnvReal("EM_APHI_TEAM7LIKE_LENGTH", 1.20);
    geometry.body_height = readEnvReal("EM_APHI_TEAM7LIKE_HEIGHT", 1.00);
    geometry.body_width = readEnvReal("EM_APHI_TEAM7LIKE_WIDTH", 0.30);
    geometry.boundary_width = readEnvReal("EM_APHI_TEAM7LIKE_BOUNDARY_WIDTH", 3.0 * dp_0);
    return geometry;
}

inline Team7LikeRegionBoxes readTeam7LikeRegionBoxesFromEnvironment(const Team7LikeGeometry &geometry)
{
    Team7LikeRegionBoxes boxes;
    const Real conductor_xmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_XMIN_FRACTION", 0.52);
    const Real conductor_xmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_XMAX_FRACTION", 0.68);
    const Real conductor_ymin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_YMIN_FRACTION", 0.30);
    const Real conductor_ymax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_YMAX_FRACTION", 0.70);
    const Real conductor_zmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_ZMIN_FRACTION", 0.20);
    const Real conductor_zmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_CONDUCTOR_ZMAX_FRACTION", 0.80);
    const Real coil_xmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_XMIN_FRACTION", 0.24);
    const Real coil_xmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_XMAX_FRACTION", 0.38);
    const Real coil_ymin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_YMIN_FRACTION", 0.15);
    const Real coil_ymax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_YMAX_FRACTION", 0.85);
    const Real coil_zmin_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_ZMIN_FRACTION", 0.20);
    const Real coil_zmax_fraction = readEnvReal("EM_APHI_TEAM7LIKE_COIL_ZMAX_FRACTION", 0.80);

    boxes.conductor_xmin = coordinateFromFraction(conductor_xmin_fraction, geometry.body_length);
    boxes.conductor_xmax = coordinateFromFraction(conductor_xmax_fraction, geometry.body_length);
    boxes.conductor_ymin = coordinateFromFraction(conductor_ymin_fraction, geometry.body_height);
    boxes.conductor_ymax = coordinateFromFraction(conductor_ymax_fraction, geometry.body_height);
    boxes.conductor_zmin = coordinateFromFraction(conductor_zmin_fraction, geometry.body_width);
    boxes.conductor_zmax = coordinateFromFraction(conductor_zmax_fraction, geometry.body_width);
    boxes.coil_xmin = coordinateFromFraction(coil_xmin_fraction, geometry.body_length);
    boxes.coil_xmax = coordinateFromFraction(coil_xmax_fraction, geometry.body_length);
    boxes.coil_ymin = coordinateFromFraction(coil_ymin_fraction, geometry.body_height);
    boxes.coil_ymax = coordinateFromFraction(coil_ymax_fraction, geometry.body_height);
    boxes.coil_zmin = coordinateFromFraction(coil_zmin_fraction, geometry.body_width);
    boxes.coil_zmax = coordinateFromFraction(coil_zmax_fraction, geometry.body_width);
    return boxes;
}

inline Team7LikeMaterialProperties readTeam7LikeMaterialPropertiesFromEnvironment()
{
    Team7LikeMaterialProperties materials;
    materials.sigma_air = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_AIR", 1.0e-4);
    materials.sigma_conductor = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_CONDUCTOR", 1.0);
    materials.sigma_coil = readEnvReal("EM_APHI_TEAM7LIKE_SIGMA_COIL", 1.0e-4);
    materials.nu_air = readEnvReal("EM_APHI_TEAM7LIKE_NU_AIR", 1.0);
    materials.nu_conductor = readEnvReal("EM_APHI_TEAM7LIKE_NU_CONDUCTOR", 1.0);
    materials.nu_coil = readEnvReal("EM_APHI_TEAM7LIKE_NU_COIL", 1.0);
    return materials;
}

inline std::string normalizeManufacturedCaseMode(const std::string &case_mode_token)
{
    if (case_mode_token == "source_free")
    {
        return "coulomb_variable_sigma_source_free";
    }
    if (case_mode_token == "driven")
    {
        return "coulomb_variable_sigma_driven";
    }
    return case_mode_token;
}

inline Team7LikeMaterialAssignment buildTeam7LikeMaterialAssignment(const Vecd *positions, size_t number_of_particles,
                                                                    const Team7LikeGeometry &geometry,
                                                                    const Team7LikeRegionBoxes &region_boxes,
                                                                    const Team7LikeMaterialProperties &materials)
{
    Team7LikeMaterialAssignment assignment;
    assignment.sigma.assign(number_of_particles, materials.sigma_air);
    assignment.nu.assign(number_of_particles, materials.nu_air);
    assignment.is_conductor.assign(number_of_particles, false);
    assignment.is_coil.assign(number_of_particles, false);
    assignment.is_source.assign(number_of_particles, false);

    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const bool conductor_region =
            insideBoxRegion(positions[i], region_boxes.conductor_xmin, region_boxes.conductor_xmax,
                            region_boxes.conductor_ymin, region_boxes.conductor_ymax, region_boxes.conductor_zmin,
                            region_boxes.conductor_zmax);
        const bool coil_region = insideBoxRegion(positions[i], region_boxes.coil_xmin, region_boxes.coil_xmax,
                                                 region_boxes.coil_ymin, region_boxes.coil_ymax, region_boxes.coil_zmin,
                                                 region_boxes.coil_zmax);
        assignment.is_conductor[i] = conductor_region;
        assignment.is_coil[i] = coil_region;
        if (conductor_region)
        {
            assignment.sigma[i] = materials.sigma_conductor;
            assignment.nu[i] = materials.nu_conductor;
            ++assignment.conductor_particles;
        }
        else if (coil_region)
        {
            assignment.sigma[i] = materials.sigma_coil;
            assignment.nu[i] = materials.nu_coil;
            ++assignment.coil_particles;
        }
        else
        {
            ++assignment.air_particles;
        }
    }
    return assignment;
}

inline void enforcePhiReference(StdVec<Complex> &phi, size_t reference_index, Complex reference_value)
{
    if (phi.empty() || reference_index >= phi.size())
    {
        return;
    }
    const Complex offset = phi[reference_index] - reference_value;
    for (Complex &value : phi)
    {
        value -= offset;
    }
#if SPHINXSYS_USE_SYCL
    matrix_free::invalidateMatrixFreeAPhiSyclValueFieldCache();
#endif
}

inline matrix_free::ScalarComplexHelmholtzSolverState
solveDiscretePhiForSourceFreeCase(const matrix_free::MatrixFreePairwiseGraph &graph, const StdVec<Real> &sigma,
                                  const StdVec<Complex> &rhs, StdVec<Complex> &phi)
{
    matrix_free::ScalarComplexHelmholtzResiduals residuals;
    residuals.resize(phi.size());
    matrix_free::ScalarComplexHelmholtzSolverParameters parameters;
    parameters.max_iterations_ = 8000;
    parameters.relaxation_factor_ = 1.0;
    parameters.absolute_tolerance_ = 1.0e-6;
    parameters.diagonal_regularization_ = 1.0e-12;

    const StdVec<Complex> zero_reaction(phi.size(), Complex(0.0, 0.0));
    return matrix_free::solveScalarComplexHelmholtz(
        phi, rhs, zero_reaction, residuals, parameters,
        [&](const StdVec<Complex> &current_field, matrix_free::ScalarComplexHelmholtzResiduals &current_residuals)
        {
            matrix_free::accumulateScalarLaplaceResidualsFromClearedGraph(graph, current_field, sigma,
                                                                          current_residuals);
        });
}

inline bool buildTeam7LikeManufacturedSetup(const std::string &case_mode,
                                            const matrix_free::MatrixFreePairwiseGraph &graph, const Vecd *positions,
                                            size_t number_of_particles, const Team7LikeGeometry &geometry,
                                            const Team7LikeRegionBoxes &region_boxes,
                                            const Team7LikeMaterialAssignment &material_assignment,
                                            const MatrixFreeAPhiParameters &parameters,
                                            const StdVec<Real> &sigma, const StdVec<Real> &nu,
                                            Team7LikeManufacturedSetup &setup, std::string &error_message)
{
    setup = Team7LikeManufacturedSetup{};
    setup.exact_fields.ax.assign(number_of_particles, Complex(0.0, 0.0));
    setup.exact_fields.ay.assign(number_of_particles, Complex(0.0, 0.0));
    setup.exact_fields.az.assign(number_of_particles, Complex(0.0, 0.0));
    setup.exact_fields.phi.assign(number_of_particles, Complex(0.0, 0.0));
    setup.sources.source_ax.assign(number_of_particles, Complex(0.0, 0.0));
    setup.sources.source_ay.assign(number_of_particles, Complex(0.0, 0.0));
    setup.sources.source_az.assign(number_of_particles, Complex(0.0, 0.0));
    setup.sources.source_phi.assign(number_of_particles, Complex(0.0, 0.0));

    const Complex transverse_amplitude(1.0, 0.1);
    const Complex coupled_ax_amplitude(0.25, -0.08);
    const Complex driven_phi_amplitude(0.12, -0.03);
    setup.has_discrete_manufactured_reference = case_mode != "coulomb_variable_sigma_forced_response";

    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const Real x = positions[i][0];
        const Real phase_x = Pi * x / geometry.body_length;
        if (case_mode == "coupled_source_free")
        {
            setup.exact_fields.ax[i] = coupled_ax_amplitude * std::sin(phase_x);
        }
        else if (case_mode != "coulomb_variable_sigma_forced_response")
        {
            setup.exact_fields.ay[i] = transverse_amplitude * std::sin(phase_x);
            if (case_mode == "coulomb_variable_sigma_driven")
            {
                const Real y = positions[i][1];
                const Real phase_y = Pi * y / geometry.body_height;
                setup.exact_fields.phi[i] = driven_phi_amplitude * std::cos(phase_x) * std::sin(phase_y);
            }
        }
    }
    enforcePhiReference(setup.exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);

    if (case_mode == "coupled_source_free" || case_mode == "coulomb_variable_sigma_source_free")
    {
        const StdVec<Vec3c> sigma_grad_ax =
            matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.ax, sigma);
        const StdVec<Vec3c> sigma_grad_ay =
            matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.ay, sigma);
        const StdVec<Vec3c> sigma_grad_az =
            matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.az, sigma);

        StdVec<Complex> rhs_phi(number_of_particles, Complex(0.0, 0.0));
        const Complex imaginary_unit(0.0, 1.0);
        for (size_t i = 0; i != number_of_particles; ++i)
        {
            const Complex divergence_sigma_a =
                sigma_grad_ax[i][0] + sigma_grad_ay[i][1] + sigma_grad_az[i][2];
            rhs_phi[i] = imaginary_unit * parameters.angular_frequency * divergence_sigma_a;
        }
        const matrix_free::ScalarComplexHelmholtzSolverState exact_phi_state =
            solveDiscretePhiForSourceFreeCase(graph, sigma, rhs_phi, setup.exact_fields.phi);
        enforcePhiReference(setup.exact_fields.phi, parameters.phi_reference_index, parameters.phi_reference_value);
        setup.reference_phi_build_residual_l2 = exact_phi_state.current_residual_l2_;
        const Real acceptable_reference_residual = 1.0e-5;
        if (!exact_phi_state.converged_ && exact_phi_state.current_residual_l2_ > acceptable_reference_residual)
        {
            error_message = "failed_to_build_source_free_phi_reference";
            return false;
        }
    }

    const StdVec<Complex> laplace_ax =
        matrix_free::applyScalarNegativeLaplaceFromGraph(graph, setup.exact_fields.ax, nu);
    const StdVec<Complex> laplace_ay =
        matrix_free::applyScalarNegativeLaplaceFromGraph(graph, setup.exact_fields.ay, nu);
    const StdVec<Complex> laplace_az =
        matrix_free::applyScalarNegativeLaplaceFromGraph(graph, setup.exact_fields.az, nu);
    const StdVec<Complex> laplace_phi_exact =
        matrix_free::applyScalarNegativeLaplaceFromGraph(graph, setup.exact_fields.phi, sigma);
    const StdVec<Vec3c> grad_phi = matrix_free::applyMatrixFreeGradient(graph, setup.exact_fields.phi);
    const StdVec<Vec3c> sigma_grad_phi =
        matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.phi, sigma);
    const StdVec<Vec3c> sigma_grad_ax_for_source =
        matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.ax, sigma);
    const StdVec<Vec3c> sigma_grad_ay_for_source =
        matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.ay, sigma);
    const StdVec<Vec3c> sigma_grad_az_for_source =
        matrix_free::applyMatrixFreeHarmonicWeightedGradient(graph, setup.exact_fields.az, sigma);

    const Complex imaginary_unit(0.0, 1.0);
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        const Complex divergence_sigma_a = sigma_grad_ax_for_source[i][0] + sigma_grad_ay_for_source[i][1] +
                                           sigma_grad_az_for_source[i][2];
        setup.sources.source_ax[i] = laplace_ax[i] + imaginary_unit * parameters.angular_frequency * sigma[i] *
                                                         setup.exact_fields.ax[i] + sigma_grad_phi[i][0];
        setup.sources.source_ay[i] = laplace_ay[i] + imaginary_unit * parameters.angular_frequency * sigma[i] *
                                                         setup.exact_fields.ay[i] + sigma_grad_phi[i][1];
        setup.sources.source_az[i] = laplace_az[i] + imaginary_unit * parameters.angular_frequency * sigma[i] *
                                                         setup.exact_fields.az[i] + sigma_grad_phi[i][2];
        setup.sources.source_phi[i] =
            laplace_phi_exact[i] - imaginary_unit * parameters.angular_frequency * divergence_sigma_a;
    }

    if (case_mode == "coulomb_variable_sigma_forced_response")
    {
        error_message = "forced_response_not_supported_in_manufactured_setup_helper";
        return false;
    }

    return true;
}

inline LaplaceStructuredAPhiBoundaryCondition buildSparseBaselineBoundaryCondition(
    const Vecd *positions, size_t number_of_particles, const Team7LikeGeometry &geometry, Real boundary_shell_thickness,
    const MatrixFreeAPhiFields &exact_fields, size_t &phi_reference_index_out)
{
    LaplaceStructuredAPhiBoundaryCondition boundary_condition;
    boundary_condition.is_dirichlet.assign(number_of_particles, false);
    boundary_condition.ax.assign(number_of_particles, Complex(0.0, 0.0));
    boundary_condition.ay.assign(number_of_particles, Complex(0.0, 0.0));
    boundary_condition.az.assign(number_of_particles, Complex(0.0, 0.0));
    boundary_condition.phi.assign(number_of_particles, Complex(0.0, 0.0));

    phi_reference_index_out = 0;
    for (size_t i = 0; i != number_of_particles; ++i)
    {
        boundary_condition.is_dirichlet[i] =
            isExteriorBoundaryParticle(positions[i], geometry, boundary_shell_thickness);
        boundary_condition.ax[i] = exact_fields.ax[i];
        boundary_condition.ay[i] = exact_fields.ay[i];
        boundary_condition.az[i] = exact_fields.az[i];
        boundary_condition.phi[i] = exact_fields.phi[i];
        if (!boundary_condition.is_dirichlet[i])
        {
            phi_reference_index_out = i;
        }
    }
    return boundary_condition;
}

} // namespace team7_like
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_TEAM7_LIKE_CASE_UTILS_HPP
