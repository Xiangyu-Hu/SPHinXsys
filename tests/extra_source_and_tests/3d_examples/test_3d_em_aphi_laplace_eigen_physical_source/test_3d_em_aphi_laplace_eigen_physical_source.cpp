/**
 * @file test_3d_em_aphi_laplace_eigen_physical_source.cpp
 * @brief Minimal physical Laplace-structured A-phi case with air, a conductive block,
 * and a smooth source-current region that can be either a strip, a full shell,
 * or a single-sided shell hugging one face of the conductor.
 */

#include "sphinxsys.h"
#include "electromagnetic_aphi_laplace_eigen.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>

using namespace SPH;

namespace
{
using Complex = electromagnetics::Complex;

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

bool get_env_bool_local(const std::string &name, bool default_value)
{
    const char *value = std::getenv(name.c_str());
    if (value == nullptr)
    {
        return default_value;
    }
    std::string token(value);
    std::transform(token.begin(), token.end(), token.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (token == "1" || token == "true" || token == "on" || token == "yes")
    {
        return true;
    }
    if (token == "0" || token == "false" || token == "off" || token == "no")
    {
        return false;
    }
    return default_value;
}

std::string get_env_string_local(const std::string &name,
                                 const std::string &default_value = "")
{
    const char *value = std::getenv(name.c_str());
    return value == nullptr ? default_value : std::string(value);
}

const Real dp_0 = get_env_real_local("EM_APHI_LAPLACE_DP", 1.0);
const Real body_length = get_env_real_local("EM_APHI_LAPLACE_LENGTH", 10.0);
const Real body_height = get_env_real_local("EM_APHI_LAPLACE_HEIGHT", 8.0);
const Real body_width = get_env_real_local("EM_APHI_LAPLACE_WIDTH", 6.0);
const Real boundary_width = get_env_real_local("EM_APHI_LAPLACE_BOUNDARY_WIDTH", 3.0 * dp_0);
const Real boundary_shell_thickness =
    get_env_real_local("EM_APHI_LAPLACE_BOUNDARY_SHELL_THICKNESS", 1.5 * dp_0);

const Real sigma_air =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SIGMA_AIR", 1.0e-8);
const Real sigma_conductor =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SIGMA_CONDUCTOR", 100.0);
const Real nu_air =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_NU_AIR", 1.0);
const Real nu_conductor =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_NU_CONDUCTOR", 1.0);
const Real frequency_hz =
    get_env_real_local("EM_APHI_LAPLACE_FREQUENCY_HZ", 10.0);
const Real source_current_amplitude =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_J0", 1000.0);

const Real conductor_x_min_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_XMIN_FRACTION", 0.40);
const Real conductor_x_max_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_XMAX_FRACTION", 0.60);
const Real conductor_y_min_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_YMIN_FRACTION", 0.25);
const Real conductor_y_max_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_YMAX_FRACTION", 0.75);
const Real conductor_z_min_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_ZMIN_FRACTION", 0.25);
const Real conductor_z_max_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_CONDUCTOR_ZMAX_FRACTION", 0.75);

const Real source_x_min_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_XMIN_FRACTION", 0.28);
const Real source_x_max_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_XMAX_FRACTION", 0.38);
const Real source_y_min_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_YMIN_FRACTION", 0.20);
const Real source_y_max_fraction =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_YMAX_FRACTION", 0.80);
const Real source_shell_thickness =
    get_env_real_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_SHELL_THICKNESS", 1.0 * dp_0);
const std::string source_mode =
    get_env_string_local("EM_APHI_LAPLACE_PHYSICAL_SOURCE_MODE", "shell");

const bool enable_gauge_penalty =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_GAUGE_PENALTY", false);
const Real gauge_penalty =
    get_env_real_local("EM_APHI_LAPLACE_GAUGE_PENALTY", 0.0);
const bool enable_block_scaling =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_BLOCK_SCALING", false);
const Real phi_block_scale =
    get_env_real_local("EM_APHI_LAPLACE_PHI_BLOCK_SCALE", 0.0);
const bool enable_diagonal_equilibration =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_DIAGONAL_EQUILIBRATION", true);
const int diagonal_equilibration_iterations = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_EQUILIBRATION_ITERATIONS", 2.0));
const std::string solver_backend =
    get_env_string_local("EM_APHI_LAPLACE_SOLVER_BACKEND", "bicgstab");
const int iterative_max_iterations = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_ITERATIVE_MAX_ITERATIONS", 5000.0));
const Real iterative_tolerance =
    get_env_real_local("EM_APHI_LAPLACE_ITERATIVE_TOLERANCE", 1.0e-8);
const bool write_particle_csv =
    get_env_bool_local("EM_APHI_LAPLACE_WRITE_PARTICLE_CSV", true);
const bool write_vtp =
    get_env_bool_local("EM_APHI_LAPLACE_PHYSICAL_WRITE_VTP", true);

const Vec3d body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_lower_bound = body_center - body_halfsize;
const Vec3d body_upper_bound = body_center + body_halfsize;
BoundingBoxd system_domain_bounds(
    Vec3d(-boundary_width, -boundary_width, -boundary_width),
    Vec3d(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

class PhysicalBoxShape : public ComplexShape
{
  public:
    explicit PhysicalBoxShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(body_center), body_halfsize);
    }
};

bool is_boundary_particle(const Vec3d &position)
{
    return (position[0] - body_lower_bound[0] < boundary_shell_thickness) ||
           (body_upper_bound[0] - position[0] < boundary_shell_thickness) ||
           (position[1] - body_lower_bound[1] < boundary_shell_thickness) ||
           (body_upper_bound[1] - position[1] < boundary_shell_thickness) ||
           (position[2] - body_lower_bound[2] < boundary_shell_thickness) ||
           (body_upper_bound[2] - position[2] < boundary_shell_thickness);
}

Real coordinate_from_fraction(Real fraction, int axis)
{
    Real lengths[3] = {body_length, body_height, body_width};
    return body_lower_bound[axis] + fraction * lengths[axis];
}

bool inside_box_region(const Vec3d &position,
                       Real xmin, Real xmax,
                       Real ymin, Real ymax,
                       Real zmin, Real zmax)
{
    return position[0] >= xmin && position[0] <= xmax &&
           position[1] >= ymin && position[1] <= ymax &&
           position[2] >= zmin && position[2] <= zmax;
}

struct PhysicalRegionSummary
{
    size_t particles = 0;
    Real mean_abs_a = 0.0;
    Real max_abs_a = 0.0;
    Real mean_abs_phi = 0.0;
    Real max_abs_phi = 0.0;
    Real mean_abs_e = 0.0;
    Real max_abs_e = 0.0;
    Real mean_abs_j = 0.0;
    Real max_abs_j = 0.0;
    Real mean_joule = 0.0;
    Real max_joule = 0.0;
};

void finalize_region_summary(PhysicalRegionSummary &summary)
{
    if (summary.particles == 0)
    {
        return;
    }
    Real inv = 1.0 / static_cast<Real>(summary.particles);
    summary.mean_abs_a *= inv;
    summary.mean_abs_phi *= inv;
    summary.mean_abs_e *= inv;
    summary.mean_abs_j *= inv;
    summary.mean_joule *= inv;
}

Real source_profile_z_current(const Vec3d &position,
                              Real xmin, Real xmax,
                              Real ymin, Real ymax,
                              Real amplitude)
{
    if (position[0] < xmin || position[0] > xmax || position[1] < ymin || position[1] > ymax)
    {
        return 0.0;
    }

    Real xi = (position[0] - xmin) / (xmax - xmin + TinyReal);
    Real eta = (position[1] - ymin) / (ymax - ymin + TinyReal);
    return amplitude * std::sin(Pi * xi) * std::sin(Pi * eta);
}

Real source_profile_z_current_shell(const Vec3d &position,
                                    Real conductor_xmin, Real conductor_xmax,
                                    Real conductor_ymin, Real conductor_ymax,
                                    Real conductor_zmin, Real conductor_zmax,
                                    Real shell_thickness,
                                    Real amplitude)
{
    if (position[2] < conductor_zmin || position[2] > conductor_zmax)
    {
        return 0.0;
    }

    bool inside_conductor_xy =
        position[0] >= conductor_xmin && position[0] <= conductor_xmax &&
        position[1] >= conductor_ymin && position[1] <= conductor_ymax;
    if (inside_conductor_xy)
    {
        return 0.0;
    }

    Real dx = 0.0;
    if (position[0] < conductor_xmin)
    {
        dx = conductor_xmin - position[0];
    }
    else if (position[0] > conductor_xmax)
    {
        dx = position[0] - conductor_xmax;
    }

    Real dy = 0.0;
    if (position[1] < conductor_ymin)
    {
        dy = conductor_ymin - position[1];
    }
    else if (position[1] > conductor_ymax)
    {
        dy = position[1] - conductor_ymax;
    }

    Real outside_distance = std::sqrt(dx * dx + dy * dy);
    if (outside_distance > shell_thickness)
    {
        return 0.0;
    }

    Real normalized = 1.0 - outside_distance / (shell_thickness + TinyReal);
    Real taper = std::sin(0.5 * Pi * normalized);
    return amplitude * taper;
}

Real source_profile_z_current_shell_left(const Vec3d &position,
                                         Real conductor_xmin,
                                         Real conductor_ymin, Real conductor_ymax,
                                         Real conductor_zmin, Real conductor_zmax,
                                         Real shell_thickness,
                                         Real amplitude)
{
    if (position[2] < conductor_zmin || position[2] > conductor_zmax)
    {
        return 0.0;
    }

    if (position[0] > conductor_xmin || position[0] < conductor_xmin - shell_thickness)
    {
        return 0.0;
    }

    Real y_lower = conductor_ymin - shell_thickness;
    Real y_upper = conductor_ymax + shell_thickness;
    if (position[1] < y_lower || position[1] > y_upper)
    {
        return 0.0;
    }

    Real dx = conductor_xmin - position[0];
    Real x_weight = std::sin(0.5 * Pi * (1.0 - dx / (shell_thickness + TinyReal)));

    Real y_mid = 0.5 * (conductor_ymin + conductor_ymax);
    Real y_half_span = 0.5 * (y_upper - y_lower);
    Real eta = (position[1] - y_mid) / (y_half_span + TinyReal);
    Real y_weight = std::sqrt(SMAX(0.0, 1.0 - eta * eta));

    return amplitude * x_weight * y_weight;
}

Real compute_source_jz(const Vec3d &position,
                       const std::string &mode_token,
                       Real source_xmin, Real source_xmax,
                       Real source_ymin, Real source_ymax,
                       Real conductor_xmin, Real conductor_xmax,
                       Real conductor_ymin, Real conductor_ymax,
                       Real conductor_zmin, Real conductor_zmax,
                       Real shell_thickness,
                       Real amplitude)
{
    std::string mode = mode_token;
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });

    if (mode == "strip")
    {
        return source_profile_z_current(position, source_xmin, source_xmax, source_ymin, source_ymax, amplitude);
    }

    if (mode == "shell_left" || mode == "single_sided_shell")
    {
        return source_profile_z_current_shell_left(position,
                                                   conductor_xmin,
                                                   conductor_ymin, conductor_ymax,
                                                   conductor_zmin, conductor_zmax,
                                                   shell_thickness, amplitude);
    }

    return source_profile_z_current_shell(position,
                                          conductor_xmin, conductor_xmax,
                                          conductor_ymin, conductor_ymax,
                                          conductor_zmin, conductor_zmax,
                                          shell_thickness, amplitude);
}

} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_APHI_LAPLACE_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-aphi-laplace-physical-config] dp=" << dp_0
              << ", sigma_air=" << sigma_air
              << ", sigma_conductor=" << sigma_conductor
              << ", nu_air=" << nu_air
              << ", nu_conductor=" << nu_conductor
              << ", frequency_hz=" << frequency_hz
              << ", source_current_amplitude=" << source_current_amplitude
              << ", source_mode=" << source_mode
              << ", source_shell_thickness=" << source_shell_thickness
              << ", enable_block_scaling=" << enable_block_scaling
              << ", phi_block_scale=" << phi_block_scale
              << ", enable_diagonal_equilibration=" << enable_diagonal_equilibration
              << ", diagonal_equilibration_iterations=" << diagonal_equilibration_iterations
              << ", solver_backend=" << solver_backend
              << ", iterative_max_iterations=" << iterative_max_iterations
              << ", iterative_tolerance=" << iterative_tolerance
              << std::endl;

    SolidBody physical_body(sph_system, makeShared<PhysicalBoxShape>("PhysicalBox"));
    physical_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    physical_body.defineMaterial<Solid>();
    physical_body.defineBodyLevelSetShape();
    physical_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(physical_body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = physical_body.getBaseParticles();
    size_t total_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Real conductor_xmin = coordinate_from_fraction(conductor_x_min_fraction, 0);
    Real conductor_xmax = coordinate_from_fraction(conductor_x_max_fraction, 0);
    Real conductor_ymin = coordinate_from_fraction(conductor_y_min_fraction, 1);
    Real conductor_ymax = coordinate_from_fraction(conductor_y_max_fraction, 1);
    Real conductor_zmin = coordinate_from_fraction(conductor_z_min_fraction, 2);
    Real conductor_zmax = coordinate_from_fraction(conductor_z_max_fraction, 2);

    Real source_xmin = coordinate_from_fraction(source_x_min_fraction, 0);
    Real source_xmax = coordinate_from_fraction(source_x_max_fraction, 0);
    Real source_ymin = coordinate_from_fraction(source_y_min_fraction, 1);
    Real source_ymax = coordinate_from_fraction(source_y_max_fraction, 1);

    StdVec<Real> sigma(total_particles, sigma_air);
    StdVec<Real> nu(total_particles, nu_air);
    StdVec<bool> is_conductor(total_particles, false);
    StdVec<bool> is_source(total_particles, false);
    StdVec<Real> source_jz(total_particles, 0.0);

    size_t conductor_particles = 0;
    size_t air_particles = 0;
    size_t source_particles = 0;

    for (size_t i = 0; i != total_particles; ++i)
    {
        bool conductor_region = inside_box_region(
            positions[i], conductor_xmin, conductor_xmax, conductor_ymin, conductor_ymax, conductor_zmin, conductor_zmax);
        is_conductor[i] = conductor_region;
        if (conductor_region)
        {
            sigma[i] = sigma_conductor;
            nu[i] = nu_conductor;
            conductor_particles++;
        }
        else
        {
            air_particles++;
        }

        Real jz = compute_source_jz(
            positions[i], source_mode,
            source_xmin, source_xmax, source_ymin, source_ymax,
            conductor_xmin, conductor_xmax, conductor_ymin, conductor_ymax, conductor_zmin, conductor_zmax,
            source_shell_thickness, source_current_amplitude);
        source_jz[i] = jz;
        is_source[i] = std::abs(jz) > TinyReal;
        if (is_source[i])
        {
            source_particles++;
        }
    }

    electromagnetics::LaplaceStructuredAPhiFields zero_fields;
    zero_fields.ax.resize(total_particles, Complex(0.0, 0.0));
    zero_fields.ay.resize(total_particles, Complex(0.0, 0.0));
    zero_fields.az.resize(total_particles, Complex(0.0, 0.0));
    zero_fields.phi.resize(total_particles, Complex(0.0, 0.0));

    electromagnetics::LaplaceStructuredAPhiBoundaryCondition boundary_condition;
    boundary_condition.is_dirichlet.resize(total_particles, false);
    boundary_condition.ax.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.ay.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.az.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.phi.resize(total_particles, Complex(0.0, 0.0));

    size_t reference_phi_index = 0;
    for (size_t i = 0; i != total_particles; ++i)
    {
        boundary_condition.is_dirichlet[i] = is_boundary_particle(positions[i]);
        if (!boundary_condition.is_dirichlet[i])
        {
            reference_phi_index = i;
        }
    }

    electromagnetics::LaplaceStructuredAPhiParameters parameters;
    parameters.angular_frequency = 2.0 * Pi * frequency_hz;
    parameters.gauge_penalty = gauge_penalty;
    parameters.enable_gauge_penalty = enable_gauge_penalty;
    parameters.enable_block_scaling = enable_block_scaling;
    parameters.phi_block_scale = phi_block_scale;
    parameters.enable_diagonal_equilibration = enable_diagonal_equilibration;
    parameters.diagonal_equilibration_iterations = diagonal_equilibration_iterations;
    parameters.solver_backend = solver_backend;
    parameters.iterative_max_iterations = iterative_max_iterations;
    parameters.iterative_tolerance = iterative_tolerance;
    parameters.fix_phi_reference = true;
    parameters.phi_reference_index = reference_phi_index;
    parameters.phi_reference_value = 0.0;

    electromagnetics::LaplaceStructuredAPhiEigenSolver solver(inner_relation, parameters);
    electromagnetics::LaplaceStructuredAPhiEigenSolver *active_solver = &solver;
    std::unique_ptr<electromagnetics::LaplaceStructuredAPhiEigenSolver> fallback_solver;
    electromagnetics::VecC rhs(static_cast<int>(4 * total_particles));
    rhs.setZero();
    for (size_t i = 0; i != total_particles; ++i)
    {
        rhs[solver.idxAz(i)] = Complex(source_jz[i], 0.0);
    }

    solver.assembleSystem(sigma, nu, rhs, boundary_condition);

    electromagnetics::LaplaceStructuredAPhiFields solved_fields;
    std::string solver_message;
    bool solve_success = solver.solve(solved_fields, solver_message);
    if (!solve_success)
    {
        std::cerr << "[em-aphi-laplace-physical] primary solve failed"
                  << ", backend=" << solver.LastSolverBackend()
                  << ", iterations=" << solver.LastSolverIterations()
                  << ", estimated_error=" << solver.LastSolverEstimatedError()
                  << ", message=" << solver_message
                  << std::endl;

        std::string requested_backend = solver_backend;
        std::transform(requested_backend.begin(), requested_backend.end(), requested_backend.begin(),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        if (requested_backend == "bicgstab")
        {
            electromagnetics::LaplaceStructuredAPhiParameters fallback_parameters = parameters;
            fallback_parameters.solver_backend = "sparse_lu";
            fallback_solver = std::make_unique<electromagnetics::LaplaceStructuredAPhiEigenSolver>(
                inner_relation, fallback_parameters);
            fallback_solver->assembleSystem(sigma, nu, rhs, boundary_condition);

            solve_success = fallback_solver->solve(solved_fields, solver_message);
            if (!solve_success)
            {
                std::cerr << "[em-aphi-laplace-physical] fallback solve failed"
                          << ", backend=" << fallback_solver->LastSolverBackend()
                          << ", iterations=" << fallback_solver->LastSolverIterations()
                          << ", estimated_error=" << fallback_solver->LastSolverEstimatedError()
                          << ", message=" << solver_message
                          << std::endl;
                return 1;
            }

            std::cout << "[em-aphi-laplace-physical] fallback solve succeeded"
                      << ", backend=" << fallback_solver->LastSolverBackend()
                      << ", message=" << solver_message
                      << std::endl;
            active_solver = fallback_solver.get();
        }
        else
        {
            return 1;
        }
    }

    electromagnetics::LaplaceStructuredAPhiDiagnostics solved_diagnostics =
        active_solver->computeDiagnostics(sigma, nu, solved_fields);

    PhysicalRegionSummary conductor_summary;
    PhysicalRegionSummary air_summary;
    PhysicalRegionSummary source_summary;

    for (size_t i = 0; i != total_particles; ++i)
    {
        Real abs_a = std::sqrt(std::norm(solved_fields.ax[i]) +
                               std::norm(solved_fields.ay[i]) +
                               std::norm(solved_fields.az[i]));
        Real abs_phi = std::abs(solved_fields.phi[i]);
        Real abs_e = solved_fields.electric_field.empty() ? 0.0 :
            static_cast<Real>(solved_fields.electric_field[i].norm());
        Real abs_j = solved_fields.current_density.empty() ? 0.0 :
            static_cast<Real>(solved_fields.current_density[i].norm());
        Real joule = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];

        auto accumulate_summary = [&](PhysicalRegionSummary &summary)
        {
            summary.particles++;
            summary.mean_abs_a += abs_a;
            summary.max_abs_a = SMAX(summary.max_abs_a, abs_a);
            summary.mean_abs_phi += abs_phi;
            summary.max_abs_phi = SMAX(summary.max_abs_phi, abs_phi);
            summary.mean_abs_e += abs_e;
            summary.max_abs_e = SMAX(summary.max_abs_e, abs_e);
            summary.mean_abs_j += abs_j;
            summary.max_abs_j = SMAX(summary.max_abs_j, abs_j);
            summary.mean_joule += joule;
            summary.max_joule = SMAX(summary.max_joule, joule);
        };

        if (is_conductor[i])
        {
            accumulate_summary(conductor_summary);
        }
        else
        {
            accumulate_summary(air_summary);
        }
        if (is_source[i])
        {
            accumulate_summary(source_summary);
        }
    }

    finalize_region_summary(conductor_summary);
    finalize_region_summary(air_summary);
    finalize_region_summary(source_summary);

    if (write_vtp)
    {
        physical_body.setNewlyUpdated();
        write_states.writeToFile(0.0);
    }

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_aphi_laplace_physical_source_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "solve_success,particles,conductor_particles,air_particles,source_particles,"
                 << "sigma_air,sigma_conductor,nu_air,nu_conductor,frequency_hz,source_current_amplitude,"
                 << "source_mode,source_shell_thickness,"
                 << "solver_backend,solver_iterations,solver_estimated_error,"
                 << "used_phi_block_scale,used_diagonal_equilibration_iterations,"
                 << "solved_div_a_l2,solved_current_continuity_l2,solved_phi_equation_residual_l2,"
                 << "solved_physical_current_continuity_l2,solved_relative_current_continuity_l2,"
                 << "solved_linear_residual_norm,solved_total_joule_power,solved_gauge_penalty_energy,solved_magnetic_energy,"
                 << "solved_max_abs_a,solved_max_abs_e,solved_max_abs_j,"
                 << "conductor_mean_abs_a,conductor_max_abs_a,conductor_mean_abs_phi,conductor_max_abs_phi,"
                 << "conductor_mean_abs_e,conductor_max_abs_e,conductor_mean_abs_j,conductor_max_abs_j,"
                 << "conductor_mean_joule,conductor_max_joule,"
                 << "air_mean_abs_a,air_max_abs_a,air_mean_abs_phi,air_max_abs_phi,"
                 << "air_mean_abs_e,air_max_abs_e,air_mean_abs_j,air_max_abs_j,"
                 << "air_mean_joule,air_max_joule,"
                 << "source_mean_abs_a,source_max_abs_a,source_mean_abs_phi,source_max_abs_phi,"
                 << "source_mean_abs_e,source_max_abs_e,source_mean_abs_j,source_max_abs_j,"
                 << "source_mean_joule,source_max_joule\n";

    summary_file << (solve_success ? 1 : 0) << ","
                 << total_particles << ","
                 << conductor_particles << ","
                 << air_particles << ","
                 << source_particles << ","
                 << sigma_air << ","
                 << sigma_conductor << ","
                 << nu_air << ","
                 << nu_conductor << ","
                 << frequency_hz << ","
                 << source_current_amplitude << ","
                 << source_mode << ","
                 << source_shell_thickness << ","
                 << active_solver->LastSolverBackend() << ","
                 << active_solver->LastSolverIterations() << ","
                 << active_solver->LastSolverEstimatedError() << ","
                 << active_solver->UsedPhiBlockScale() << ","
                 << active_solver->UsedDiagonalEquilibrationIterations() << ","
                 << solved_diagnostics.divergence_a_l2 << ","
                 << solved_diagnostics.current_continuity_l2 << ","
                 << solved_diagnostics.phi_equation_residual_l2 << ","
                 << solved_diagnostics.physical_current_continuity_l2 << ","
                 << solved_diagnostics.relative_current_continuity_l2 << ","
                 << solved_diagnostics.linear_residual_norm << ","
                 << solved_diagnostics.total_joule_power << ","
                 << solved_diagnostics.gauge_penalty_energy << ","
                 << solved_diagnostics.magnetic_energy << ","
                 << solved_diagnostics.max_abs_a << ","
                 << solved_diagnostics.max_abs_e << ","
                 << solved_diagnostics.max_abs_j << ","
                 << conductor_summary.mean_abs_a << ","
                 << conductor_summary.max_abs_a << ","
                 << conductor_summary.mean_abs_phi << ","
                 << conductor_summary.max_abs_phi << ","
                 << conductor_summary.mean_abs_e << ","
                 << conductor_summary.max_abs_e << ","
                 << conductor_summary.mean_abs_j << ","
                 << conductor_summary.max_abs_j << ","
                 << conductor_summary.mean_joule << ","
                 << conductor_summary.max_joule << ","
                 << air_summary.mean_abs_a << ","
                 << air_summary.max_abs_a << ","
                 << air_summary.mean_abs_phi << ","
                 << air_summary.max_abs_phi << ","
                 << air_summary.mean_abs_e << ","
                 << air_summary.max_abs_e << ","
                 << air_summary.mean_abs_j << ","
                 << air_summary.max_abs_j << ","
                 << air_summary.mean_joule << ","
                 << air_summary.max_joule << ","
                 << source_summary.mean_abs_a << ","
                 << source_summary.max_abs_a << ","
                 << source_summary.mean_abs_phi << ","
                 << source_summary.max_abs_phi << ","
                 << source_summary.mean_abs_e << ","
                 << source_summary.max_abs_e << ","
                 << source_summary.mean_abs_j << ","
                 << source_summary.max_abs_j << ","
                 << source_summary.mean_joule << ","
                 << source_summary.max_joule << "\n";

    if (write_particle_csv)
    {
        std::ofstream particle_file(
            io_environment.OutputFolder() + "/em_aphi_laplace_physical_source_particles.csv",
            std::ios::out | std::ios::trunc);
        particle_file << std::setprecision(12);
        particle_file << "particle_id,x,y,z,is_boundary,is_conductor,is_source,sigma,nu,source_jz,"
                      << "ax_real,ax_imag,ay_real,ay_imag,az_real,az_imag,phi_real,phi_imag,"
                      << "abs_a,abs_phi,abs_e,abs_j,joule_heat\n";
        for (size_t i = 0; i != total_particles; ++i)
        {
            Real abs_a = std::sqrt(std::norm(solved_fields.ax[i]) +
                                   std::norm(solved_fields.ay[i]) +
                                   std::norm(solved_fields.az[i]));
            Real abs_phi = std::abs(solved_fields.phi[i]);
            Real abs_e = solved_fields.electric_field.empty() ? 0.0 :
                static_cast<Real>(solved_fields.electric_field[i].norm());
            Real abs_j = solved_fields.current_density.empty() ? 0.0 :
                static_cast<Real>(solved_fields.current_density[i].norm());
            Real joule = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];

            particle_file << i << ","
                          << positions[i][0] << ","
                          << positions[i][1] << ","
                          << positions[i][2] << ","
                          << (boundary_condition.is_dirichlet[i] ? 1 : 0) << ","
                          << (is_conductor[i] ? 1 : 0) << ","
                          << (is_source[i] ? 1 : 0) << ","
                          << sigma[i] << ","
                          << nu[i] << ","
                          << source_jz[i] << ","
                          << solved_fields.ax[i].real() << ","
                          << solved_fields.ax[i].imag() << ","
                          << solved_fields.ay[i].real() << ","
                          << solved_fields.ay[i].imag() << ","
                          << solved_fields.az[i].real() << ","
                          << solved_fields.az[i].imag() << ","
                          << solved_fields.phi[i].real() << ","
                          << solved_fields.phi[i].imag() << ","
                          << abs_a << ","
                          << abs_phi << ","
                          << abs_e << ","
                          << abs_j << ","
                          << joule << "\n";
        }
    }

    std::cout << "[em-aphi-laplace-physical] " << solver_message
              << ", solver_backend=" << active_solver->LastSolverBackend()
              << ", solver_iterations=" << active_solver->LastSolverIterations()
              << ", solver_estimated_error=" << active_solver->LastSolverEstimatedError()
              << ", conductor_particles=" << conductor_particles
              << ", source_particles=" << source_particles
              << ", solved_total_joule_power=" << solved_diagnostics.total_joule_power
              << ", solved_max_abs_j=" << solved_diagnostics.max_abs_j
              << ", conductor_mean_abs_j=" << conductor_summary.mean_abs_j
              << ", conductor_mean_joule=" << conductor_summary.mean_joule
              << ", solved_linear_residual_norm=" << solved_diagnostics.linear_residual_norm
              << std::endl;

    return 0;
}
