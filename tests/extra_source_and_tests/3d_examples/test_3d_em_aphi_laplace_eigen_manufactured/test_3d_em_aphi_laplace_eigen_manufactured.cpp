/**
 * @file test_3d_em_aphi_laplace_eigen_manufactured.cpp
 * @brief Minimal manufactured-solution verification for the Laplace-structured
 * frequency-domain A-phi prototype assembled as one global complex linear system.
 */

#include "sphinxsys.h"
#include "electromagnetic_aphi_laplace_eigen.hpp"
#include "electromagnetic_aphi_operator_comparison_eigen.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>

using namespace SPH;

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
const Real conductivity = get_env_real_local("EM_APHI_LAPLACE_SIGMA", 1.0);
const Real magnetic_reluctivity = get_env_real_local("EM_APHI_LAPLACE_NU", 1.0);
const Real frequency_hz = get_env_real_local("EM_APHI_LAPLACE_FREQUENCY_HZ", 1.0);
const Real gauge_penalty = get_env_real_local("EM_APHI_LAPLACE_GAUGE_PENALTY", 0.0);
const bool enable_gauge_penalty =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_GAUGE_PENALTY", false);
const bool enable_block_scaling =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_BLOCK_SCALING", false);
const Real phi_block_scale =
    get_env_real_local("EM_APHI_LAPLACE_PHI_BLOCK_SCALE", 0.0);
const bool enable_diagonal_equilibration =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_DIAGONAL_EQUILIBRATION", false);
const int diagonal_equilibration_iterations = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_EQUILIBRATION_ITERATIONS", 2.0));
const std::string solver_backend =
    get_env_string_local("EM_APHI_LAPLACE_SOLVER_BACKEND", "sparse_lu");
const int iterative_max_iterations = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_ITERATIVE_MAX_ITERATIONS", 5000.0));
const Real iterative_tolerance =
    get_env_real_local("EM_APHI_LAPLACE_ITERATIVE_TOLERANCE", 1.0e-8);
const bool enable_layered_material =
    get_env_bool_local("EM_APHI_LAPLACE_ENABLE_LAYERED_MATERIAL", false);
const int layered_material_axis = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_LAYER_AXIS", 0.0));
const Real layered_material_split_fraction =
    get_env_real_local("EM_APHI_LAPLACE_LAYER_SPLIT_FRACTION", 0.5);
const Real sigma_upper_ratio =
    get_env_real_local("EM_APHI_LAPLACE_SIGMA_UPPER_RATIO", 10.0);
const Real nu_upper_ratio =
    get_env_real_local("EM_APHI_LAPLACE_NU_UPPER_RATIO", 1.0);
const Real interface_shell_thickness =
    get_env_real_local("EM_APHI_LAPLACE_INTERFACE_SHELL_THICKNESS", 2.0 * dp_0);
const bool write_particle_csv =
    get_env_bool_local("EM_APHI_LAPLACE_WRITE_PARTICLE_CSV", true);
const std::string rhs_mode_token =
    get_env_string_local("EM_APHI_LAPLACE_RHS_MODE", "discrete_manufactured");
const std::string case_mode_token =
    get_env_string_local("EM_APHI_LAPLACE_CASE_MODE", "general_complex");

const Vec3d body_halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
const Vec3d body_lower_bound = body_center - body_halfsize;
const Vec3d body_upper_bound = body_center + body_halfsize;
BoundingBoxd system_domain_bounds(
    Vec3d(-boundary_width, -boundary_width, -boundary_width),
    Vec3d(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

class ConductorShape : public ComplexShape
{
  public:
    explicit ConductorShape(const std::string &shape_name) : ComplexShape(shape_name)
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

enum class ManufacturedRhsMode
{
    discrete_manufactured,
    continuum_manufactured
};

enum class ManufacturedCaseMode
{
    general_complex,
    source_free_continuity
};

ManufacturedRhsMode parse_rhs_mode(const std::string &token)
{
    if (token == "continuum" || token == "continuum_manufactured")
    {
        return ManufacturedRhsMode::continuum_manufactured;
    }
    return ManufacturedRhsMode::discrete_manufactured;
}

ManufacturedCaseMode parse_case_mode(const std::string &token)
{
    if (token == "source_free" || token == "source_free_continuity")
    {
        return ManufacturedCaseMode::source_free_continuity;
    }
    return ManufacturedCaseMode::general_complex;
}

struct ManufacturedPointData
{
    electromagnetics::Complex ax = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex ay = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex az = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex phi = electromagnetics::Complex(0.0, 0.0);

    electromagnetics::Complex grad_phi_x = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex grad_phi_y = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex grad_phi_z = electromagnetics::Complex(0.0, 0.0);

    electromagnetics::Complex rhs_ax = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex rhs_ay = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex rhs_az = electromagnetics::Complex(0.0, 0.0);
    electromagnetics::Complex rhs_phi = electromagnetics::Complex(0.0, 0.0);
};

ManufacturedPointData evaluate_manufactured_point(const Vec3d &position,
                                                  Real electrical_conductivity,
                                                  Real magnetic_reluctivity,
                                                  Real angular_frequency,
                                                  ManufacturedCaseMode case_mode)
{
    ManufacturedPointData point_data;

    Real x = position[0] - body_lower_bound[0];
    Real y = position[1] - body_lower_bound[1];
    Real z = position[2] - body_lower_bound[2];

    Real kx = Pi / body_length;
    Real ky = Pi / body_height;
    Real kz = Pi / body_width;

    Real sx = std::sin(kx * x);
    Real cx = std::cos(kx * x);
    Real sy = std::sin(ky * y);
    Real cy = std::cos(ky * y);
    Real sz = std::sin(kz * z);
    Real cz = std::cos(kz * z);

    Real s2x = std::sin(2.0 * kx * x);
    Real c2x = std::cos(2.0 * kx * x);
    Real s2y = std::sin(2.0 * ky * y);
    Real c2y = std::cos(2.0 * ky * y);
    Real s2z = std::sin(2.0 * kz * z);
    Real c2z = std::cos(2.0 * kz * z);

    electromagnetics::Complex psi_scale(1.0, 0.2);
    Real base_stream = sx * sy * sz;
    point_data.ax = -psi_scale * (ky * sx * cy * sz);
    point_data.ay = psi_scale * (kx * cx * sy * sz);
    point_data.az = electromagnetics::Complex(0.0, 0.0);

    if (case_mode == ManufacturedCaseMode::source_free_continuity)
    {
        point_data.phi = electromagnetics::Complex(0.0, 0.0);
        point_data.grad_phi_x = electromagnetics::Complex(0.0, 0.0);
        point_data.grad_phi_y = electromagnetics::Complex(0.0, 0.0);
        point_data.grad_phi_z = electromagnetics::Complex(0.0, 0.0);
    }
    else
    {
        Real phi_real = 0.15 * sx * s2y * sz;
        Real phi_imag = 0.08 * s2x * sy * s2z;
        point_data.phi = electromagnetics::Complex(phi_real, phi_imag);
        point_data.grad_phi_x = electromagnetics::Complex(0.15 * kx * cx * s2y * sz,
                                                          0.16 * kx * c2x * sy * s2z);
        point_data.grad_phi_y = electromagnetics::Complex(0.30 * ky * sx * c2y * sz,
                                                          0.08 * ky * s2x * cy * s2z);
        point_data.grad_phi_z = electromagnetics::Complex(0.15 * kz * sx * s2y * cz,
                                                          0.16 * kz * s2x * sy * c2z);
    }

    Real lambda_stream = kx * kx + ky * ky + kz * kz;
    electromagnetics::Complex minus_nu_laplace_ax =
        magnetic_reluctivity * lambda_stream * point_data.ax;
    electromagnetics::Complex minus_nu_laplace_ay =
        magnetic_reluctivity * lambda_stream * point_data.ay;
    electromagnetics::Complex minus_nu_laplace_az(0.0, 0.0);

    electromagnetics::Complex minus_sigma_laplace_phi(0.0, 0.0);
    if (case_mode != ManufacturedCaseMode::source_free_continuity)
    {
        Real phi_real = point_data.phi.real();
        Real phi_imag = point_data.phi.imag();
        Real lambda_phi_real = kx * kx + 4.0 * ky * ky + kz * kz;
        Real lambda_phi_imag = 4.0 * kx * kx + ky * ky + 4.0 * kz * kz;
        minus_sigma_laplace_phi =
            electromagnetics::Complex(electrical_conductivity * lambda_phi_real * phi_real,
                                      electrical_conductivity * lambda_phi_imag * phi_imag);
    }

    electromagnetics::Complex jw_sigma(0.0, angular_frequency * electrical_conductivity);
    point_data.rhs_ax = minus_nu_laplace_ax + jw_sigma * point_data.ax +
                        electrical_conductivity * point_data.grad_phi_x;
    point_data.rhs_ay = minus_nu_laplace_ay + jw_sigma * point_data.ay +
                        electrical_conductivity * point_data.grad_phi_y;
    point_data.rhs_az = minus_nu_laplace_az + jw_sigma * point_data.az +
                        electrical_conductivity * point_data.grad_phi_z;
    point_data.rhs_phi = minus_sigma_laplace_phi;

    return point_data;
}

int parity_from_position(const Vec3d &position)
{
    int ix = static_cast<int>(std::llround((position[0] - 0.5 * dp_0) / dp_0));
    int iy = static_cast<int>(std::llround((position[1] - 0.5 * dp_0) / dp_0));
    int iz = static_cast<int>(std::llround((position[2] - 0.5 * dp_0) / dp_0));
    return ((ix + iy + iz) % 2 == 0) ? 1 : -1;
}

struct ErrorSummary
{
    Real mean_a_error = 0.0;
    Real max_a_error = 0.0;
    Real l2_a_error = 0.0;
    Real mean_phi_error = 0.0;
    Real max_phi_error = 0.0;
    Real l2_phi_error = 0.0;
    Real mean_abs_j_error = 0.0;
    Real max_abs_j_error = 0.0;
    Real l2_abs_j_error = 0.0;
    Real mean_rel_abs_j_error = 0.0;
    Real mean_joule_error = 0.0;
    Real max_joule_error = 0.0;
    Real l2_joule_error = 0.0;
    Real mean_rel_joule_error = 0.0;
};

ErrorSummary evaluate_error_summary(BaseParticles &particles,
                                    const electromagnetics::LaplaceStructuredAPhiFields &exact_fields,
                                    const electromagnetics::LaplaceStructuredAPhiFields &solved_fields)
{
    ErrorSummary summary;
    size_t total_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_particles; ++i)
    {
        Real a_error = std::sqrt(std::norm(exact_fields.ax[i] - solved_fields.ax[i]) +
                                 std::norm(exact_fields.ay[i] - solved_fields.ay[i]) +
                                 std::norm(exact_fields.az[i] - solved_fields.az[i]));
        Real phi_error = std::abs(exact_fields.phi[i] - solved_fields.phi[i]);
        Real exact_abs_j = exact_fields.current_density.empty() ? 0.0 :
            static_cast<Real>(exact_fields.current_density[i].norm());
        Real solved_abs_j = solved_fields.current_density.empty() ? 0.0 :
            static_cast<Real>(solved_fields.current_density[i].norm());
        Real abs_j_error = std::abs(solved_abs_j - exact_abs_j);
        Real rel_abs_j_error = abs_j_error / (exact_abs_j + TinyReal);
        Real exact_joule = exact_fields.joule_heat.empty() ? 0.0 : exact_fields.joule_heat[i];
        Real solved_joule = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];
        Real joule_error = std::abs(solved_joule - exact_joule);
        Real rel_joule_error = joule_error / (std::abs(exact_joule) + TinyReal);

        summary.mean_a_error += a_error;
        summary.mean_phi_error += phi_error;
        summary.l2_a_error += a_error * a_error;
        summary.l2_phi_error += phi_error * phi_error;
        summary.max_a_error = SMAX(summary.max_a_error, a_error);
        summary.max_phi_error = SMAX(summary.max_phi_error, phi_error);
        summary.mean_abs_j_error += abs_j_error;
        summary.max_abs_j_error = SMAX(summary.max_abs_j_error, abs_j_error);
        summary.l2_abs_j_error += abs_j_error * abs_j_error;
        summary.mean_rel_abs_j_error += rel_abs_j_error;
        summary.mean_joule_error += joule_error;
        summary.max_joule_error = SMAX(summary.max_joule_error, joule_error);
        summary.l2_joule_error += joule_error * joule_error;
        summary.mean_rel_joule_error += rel_joule_error;
    }
    summary.mean_a_error /= static_cast<Real>(total_particles) + TinyReal;
    summary.mean_phi_error /= static_cast<Real>(total_particles) + TinyReal;
    summary.l2_a_error = std::sqrt(summary.l2_a_error / (static_cast<Real>(total_particles) + TinyReal));
    summary.l2_phi_error = std::sqrt(summary.l2_phi_error / (static_cast<Real>(total_particles) + TinyReal));
    summary.mean_abs_j_error /= static_cast<Real>(total_particles) + TinyReal;
    summary.l2_abs_j_error = std::sqrt(summary.l2_abs_j_error / (static_cast<Real>(total_particles) + TinyReal));
    summary.mean_rel_abs_j_error /= static_cast<Real>(total_particles) + TinyReal;
    summary.mean_joule_error /= static_cast<Real>(total_particles) + TinyReal;
    summary.l2_joule_error = std::sqrt(summary.l2_joule_error / (static_cast<Real>(total_particles) + TinyReal));
    summary.mean_rel_joule_error /= static_cast<Real>(total_particles) + TinyReal;
    return summary;
}

int clamped_material_axis(int requested_axis)
{
    return requested_axis < 0 ? 0 : (requested_axis > 2 ? 2 : requested_axis);
}

Real material_split_coordinate(int axis, Real split_fraction)
{
    Real safe_fraction = SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), split_fraction));
    return body_lower_bound[axis] + safe_fraction * (body_upper_bound[axis] - body_lower_bound[axis]);
}

bool is_upper_material(const Vec3d &position, int axis, Real split_coordinate)
{
    return position[axis] >= split_coordinate;
}

struct MaterialSummary
{
    Real lower_sigma = 0.0;
    Real upper_sigma = 0.0;
    Real lower_nu = 0.0;
    Real upper_nu = 0.0;
    size_t lower_particles = 0;
    size_t upper_particles = 0;
};

struct InterfaceRegionSummary
{
    size_t particles = 0;
    Real mean_a_error = 0.0;
    Real max_a_error = 0.0;
    Real mean_phi_error = 0.0;
    Real max_phi_error = 0.0;
    Real mean_abs_j = 0.0;
    Real max_abs_j = 0.0;
    Real mean_exact_abs_j = 0.0;
    Real mean_abs_j_error = 0.0;
    Real mean_rel_abs_j_error = 0.0;
    Real mean_joule = 0.0;
    Real max_joule = 0.0;
    Real mean_exact_joule = 0.0;
    Real mean_joule_error = 0.0;
    Real mean_rel_joule_error = 0.0;
};

Real interface_distance(const Vec3d &position, int axis, Real split_coordinate)
{
    return std::abs(position[axis] - split_coordinate);
}

bool is_interface_shell_particle(const Vec3d &position, int axis, Real split_coordinate, Real shell_thickness)
{
    return interface_distance(position, axis, split_coordinate) <= shell_thickness;
}

void finalize_interface_region_summary(InterfaceRegionSummary &summary)
{
    if (summary.particles == 0)
    {
        return;
    }
    Real inv_count = 1.0 / (static_cast<Real>(summary.particles) + TinyReal);
    summary.mean_a_error *= inv_count;
    summary.mean_phi_error *= inv_count;
    summary.mean_abs_j *= inv_count;
    summary.mean_exact_abs_j *= inv_count;
    summary.mean_abs_j_error *= inv_count;
    summary.mean_rel_abs_j_error *= inv_count;
    summary.mean_joule *= inv_count;
    summary.mean_exact_joule *= inv_count;
    summary.mean_joule_error *= inv_count;
    summary.mean_rel_joule_error *= inv_count;
}
} // namespace

int main(int ac, char *av[])
{
    ManufacturedRhsMode rhs_mode = parse_rhs_mode(rhs_mode_token);
    ManufacturedCaseMode case_mode = parse_case_mode(case_mode_token);
    ManufacturedRhsMode effective_rhs_mode = rhs_mode;
    std::string effective_rhs_mode_token = rhs_mode_token;
    if (case_mode == ManufacturedCaseMode::source_free_continuity)
    {
        effective_rhs_mode = ManufacturedRhsMode::discrete_manufactured;
        effective_rhs_mode_token = "discrete_manufactured";
    }

    if (effective_rhs_mode == ManufacturedRhsMode::continuum_manufactured && enable_layered_material)
    {
        std::cerr << "[em-aphi-laplace] continuum_manufactured mode currently supports constant coefficients only."
                  << std::endl;
        return 2;
    }

    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_APHI_LAPLACE_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-aphi-laplace-config] dp=" << dp_0
              << ", sigma=" << conductivity
              << ", nu=" << magnetic_reluctivity
              << ", frequency_hz=" << frequency_hz
              << ", gauge_penalty=" << gauge_penalty
              << ", enable_gauge_penalty=" << enable_gauge_penalty
              << ", enable_block_scaling=" << enable_block_scaling
              << ", phi_block_scale=" << phi_block_scale
              << ", enable_diagonal_equilibration=" << enable_diagonal_equilibration
              << ", diagonal_equilibration_iterations=" << diagonal_equilibration_iterations
              << ", solver_backend=" << solver_backend
              << ", iterative_max_iterations=" << iterative_max_iterations
              << ", iterative_tolerance=" << iterative_tolerance
              << ", enable_layered_material=" << enable_layered_material
              << ", layer_axis=" << layered_material_axis
              << ", layer_split_fraction=" << layered_material_split_fraction
              << ", sigma_upper_ratio=" << sigma_upper_ratio
              << ", nu_upper_ratio=" << nu_upper_ratio
              << ", rhs_mode=" << rhs_mode_token
              << ", effective_rhs_mode=" << effective_rhs_mode_token
              << ", case_mode=" << case_mode_token
              << ", interface_shell_thickness=" << interface_shell_thickness
              << std::endl;

    SolidBody conductor_body(sph_system, makeShared<ConductorShape>("Conductor"));
    conductor_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    conductor_body.defineMaterial<Solid>();
    conductor_body.defineBodyLevelSetShape();
    conductor_body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(conductor_body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    BaseParticles &particles = conductor_body.getBaseParticles();
    size_t total_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    StdVec<Real> sigma(total_particles, conductivity);
    StdVec<Real> nu(total_particles, magnetic_reluctivity);
    MaterialSummary material_summary;
    material_summary.lower_sigma = conductivity;
    material_summary.upper_sigma = conductivity * sigma_upper_ratio;
    material_summary.lower_nu = magnetic_reluctivity;
    material_summary.upper_nu = magnetic_reluctivity * nu_upper_ratio;
    int material_axis = clamped_material_axis(layered_material_axis);
    Real split_coordinate = material_split_coordinate(material_axis, layered_material_split_fraction);

    electromagnetics::LaplaceStructuredAPhiFields exact_fields;
    exact_fields.ax.resize(total_particles);
    exact_fields.ay.resize(total_particles);
    exact_fields.az.resize(total_particles);
    exact_fields.phi.resize(total_particles);

    electromagnetics::LaplaceStructuredAPhiBoundaryCondition boundary_condition;
    boundary_condition.is_dirichlet.resize(total_particles, false);
    boundary_condition.ax.resize(total_particles, electromagnetics::Complex(0.0, 0.0));
    boundary_condition.ay.resize(total_particles, electromagnetics::Complex(0.0, 0.0));
    boundary_condition.az.resize(total_particles, electromagnetics::Complex(0.0, 0.0));
    boundary_condition.phi.resize(total_particles, electromagnetics::Complex(0.0, 0.0));

    size_t reference_phi_index = 0;
    for (size_t i = 0; i != total_particles; ++i)
    {
        if (enable_layered_material)
        {
            bool upper_material = is_upper_material(positions[i], material_axis, split_coordinate);
            sigma[i] = upper_material ? material_summary.upper_sigma : material_summary.lower_sigma;
            nu[i] = upper_material ? material_summary.upper_nu : material_summary.lower_nu;
            if (upper_material)
            {
                material_summary.upper_particles++;
            }
            else
            {
                material_summary.lower_particles++;
            }
        }
        else
        {
            material_summary.lower_particles++;
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

    electromagnetics::VecC rhs_exact(static_cast<int>(4 * total_particles));
    rhs_exact.setZero();
    for (size_t i = 0; i != total_particles; ++i)
    {
        ManufacturedPointData point_data =
            evaluate_manufactured_point(positions[i], sigma[i], nu[i], parameters.angular_frequency, case_mode);

        exact_fields.ax[i] = point_data.ax;
        exact_fields.ay[i] = point_data.ay;
        exact_fields.az[i] = point_data.az;
        exact_fields.phi[i] = point_data.phi;

        boundary_condition.is_dirichlet[i] = is_boundary_particle(positions[i]);
        boundary_condition.ax[i] = exact_fields.ax[i];
        boundary_condition.ay[i] = exact_fields.ay[i];
        boundary_condition.az[i] = exact_fields.az[i];
        boundary_condition.phi[i] = exact_fields.phi[i];

        if (!boundary_condition.is_dirichlet[i])
        {
            reference_phi_index = i;
        }

        if (effective_rhs_mode == ManufacturedRhsMode::continuum_manufactured)
        {
            rhs_exact[static_cast<int>(i)] = point_data.rhs_ax;
            rhs_exact[static_cast<int>(total_particles + i)] = point_data.rhs_ay;
            rhs_exact[static_cast<int>(2 * total_particles + i)] = point_data.rhs_az;
            rhs_exact[static_cast<int>(3 * total_particles + i)] = point_data.rhs_phi;
        }
    }
    parameters.phi_reference_index = reference_phi_index;
    electromagnetics::LaplaceStructuredAPhiEigenSolver solver(inner_relation, parameters);
    if (case_mode == ManufacturedCaseMode::source_free_continuity)
    {
        std::string compatible_phi_message;
        bool compatible_phi_success =
            solver.solveDiscreteCompatiblePhi(sigma, nu, boundary_condition, exact_fields, compatible_phi_message);
        if (!compatible_phi_success)
        {
            std::cerr << "[em-aphi-laplace] " << compatible_phi_message << std::endl;
            return 1;
        }
        for (size_t i = 0; i != total_particles; ++i)
        {
            boundary_condition.phi[i] = exact_fields.phi[i];
        }
    }

    if (effective_rhs_mode == ManufacturedRhsMode::discrete_manufactured)
    {
        rhs_exact = solver.assembleRightHandSideFromExactSolution(sigma, nu, exact_fields);
    }
    solver.assembleSystem(sigma, nu, rhs_exact, boundary_condition);
    electromagnetics::VecC exact_residual =
        solver.systemMatrix() * solver.gatherScaledUnknownVector(exact_fields) - solver.rightHandSide();
    Real exact_linear_residual_norm = exact_residual.norm();

    electromagnetics::LaplaceStructuredAPhiFields solved_fields;
    std::string solver_message;
    bool solve_success = solver.solve(solved_fields, solver_message);
    if (!solve_success)
    {
        std::cerr << "[em-aphi-laplace] " << solver_message << std::endl;
        return 1;
    }

    electromagnetics::LaplaceStructuredAPhiDiagnostics exact_diagnostics =
        solver.computeDiagnostics(sigma, nu, exact_fields);
    electromagnetics::LaplaceStructuredAPhiDiagnostics solved_diagnostics =
        solver.computeDiagnostics(sigma, nu, solved_fields);
    ErrorSummary error_summary =
        evaluate_error_summary(particles, exact_fields, solved_fields);
    InterfaceRegionSummary interface_shell_summary;
    InterfaceRegionSummary non_interface_summary;
    for (size_t i = 0; i != total_particles; ++i)
    {
        Real a_error = std::sqrt(std::norm(exact_fields.ax[i] - solved_fields.ax[i]) +
                                 std::norm(exact_fields.ay[i] - solved_fields.ay[i]) +
                                 std::norm(exact_fields.az[i] - solved_fields.az[i]));
        Real phi_error = std::abs(exact_fields.phi[i] - solved_fields.phi[i]);
        Real exact_abs_j = exact_fields.current_density.empty() ? 0.0 :
            static_cast<Real>(exact_fields.current_density[i].norm());
        Real abs_j = solved_fields.current_density.empty() ? 0.0 :
            static_cast<Real>(solved_fields.current_density[i].norm());
        Real abs_j_error = std::abs(abs_j - exact_abs_j);
        Real rel_abs_j_error = abs_j_error / (exact_abs_j + TinyReal);
        Real exact_joule = exact_fields.joule_heat.empty() ? 0.0 : exact_fields.joule_heat[i];
        Real joule = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];
        Real joule_error = std::abs(joule - exact_joule);
        Real rel_joule_error = joule_error / (std::abs(exact_joule) + TinyReal);

        bool in_interface_shell = enable_layered_material &&
                                  is_interface_shell_particle(positions[i], material_axis,
                                                              split_coordinate, interface_shell_thickness);
        InterfaceRegionSummary &region_summary =
            in_interface_shell ? interface_shell_summary : non_interface_summary;
        region_summary.particles++;
        region_summary.mean_a_error += a_error;
        region_summary.max_a_error = SMAX(region_summary.max_a_error, a_error);
        region_summary.mean_phi_error += phi_error;
        region_summary.max_phi_error = SMAX(region_summary.max_phi_error, phi_error);
        region_summary.mean_abs_j += abs_j;
        region_summary.max_abs_j = SMAX(region_summary.max_abs_j, abs_j);
        region_summary.mean_exact_abs_j += exact_abs_j;
        region_summary.mean_abs_j_error += abs_j_error;
        region_summary.mean_rel_abs_j_error += rel_abs_j_error;
        region_summary.mean_joule += joule;
        region_summary.max_joule = SMAX(region_summary.max_joule, joule);
        region_summary.mean_exact_joule += exact_joule;
        region_summary.mean_joule_error += joule_error;
        region_summary.mean_rel_joule_error += rel_joule_error;
    }
    finalize_interface_region_summary(interface_shell_summary);
    finalize_interface_region_summary(non_interface_summary);

    StdVec<electromagnetics::Complex> checker_ax_field(total_particles);
    StdVec<electromagnetics::Complex> checker_ay_field(total_particles, electromagnetics::Complex(0.0, 0.0));
    StdVec<electromagnetics::Complex> checker_az_field(total_particles, electromagnetics::Complex(0.0, 0.0));
    for (size_t i = 0; i != total_particles; ++i)
    {
        checker_ax_field[i] = electromagnetics::Complex(
            static_cast<Real>(parity_from_position(positions[i])), 0.0);
    }

    electromagnetics::DiscreteOperatorComparisonParameters comparison_parameters;
    comparison_parameters.gauge_penalty = gauge_penalty;
    comparison_parameters.enable_gauge_penalty = enable_gauge_penalty;

    electromagnetics::DiscreteEMOperatorComparator comparator(inner_relation, comparison_parameters);
    comparator.buildOperators(nu);
    electromagnetics::DiscreteModeEnergySummary smooth_mode_summary =
        comparator.evaluateMode(exact_fields.ax, exact_fields.ay, exact_fields.az);
    electromagnetics::DiscreteModeEnergySummary checker_mode_summary =
        comparator.evaluateMode(checker_ax_field, checker_ay_field, checker_az_field);

    auto normalized_energy = [](const electromagnetics::DiscreteModeEnergySummary &summary,
                                Real energy) -> Real
    {
        return energy / (summary.max_abs_a * summary.max_abs_a + TinyReal);
    };

    Real smooth_weak_to_laplace_ratio =
        normalized_energy(smooth_mode_summary, smooth_mode_summary.weak_curlcurl_energy) /
        (normalized_energy(smooth_mode_summary, smooth_mode_summary.laplace_energy) + TinyReal);
    Real checker_weak_to_laplace_ratio =
        normalized_energy(checker_mode_summary, checker_mode_summary.weak_curlcurl_energy) /
        (normalized_energy(checker_mode_summary, checker_mode_summary.laplace_energy) + TinyReal);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_aphi_laplace_eigen_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "solve_success,"
                 << "rhs_mode,"
                 << "case_mode,"
                 << "particles,"
                 << "enable_block_scaling,"
                 << "used_phi_block_scale,"
                 << "enable_diagonal_equilibration,"
                 << "used_diagonal_equilibration_iterations,"
                 << "solver_backend,"
                 << "solver_iterations,"
                 << "solver_estimated_error,"
                 << "enable_layered_material,"
                 << "layer_axis,"
                 << "layer_split_fraction,"
                 << "lower_sigma,"
                 << "upper_sigma,"
                 << "lower_nu,"
                 << "upper_nu,"
                 << "lower_material_particles,"
                 << "upper_material_particles,"
                 << "exact_linear_residual_norm,"
                 << "mean_a_error,"
                 << "max_a_error,"
                 << "l2_a_error,"
                 << "mean_phi_error,"
                 << "max_phi_error,"
                 << "l2_phi_error,"
                 << "mean_abs_j_error,"
                 << "max_abs_j_error,"
                 << "l2_abs_j_error,"
                 << "mean_rel_abs_j_error,"
                 << "mean_joule_error,"
                 << "max_joule_error,"
                 << "l2_joule_error,"
                 << "mean_rel_joule_error,"
                 << "exact_div_a_l2,"
                 << "solved_div_a_l2,"
                 << "exact_current_continuity_l2,"
                 << "solved_current_continuity_l2,"
                 << "exact_phi_equation_residual_l2,"
                 << "solved_phi_equation_residual_l2,"
                 << "exact_physical_current_continuity_l2,"
                 << "solved_physical_current_continuity_l2,"
                 << "exact_relative_current_continuity_l2,"
                 << "solved_relative_current_continuity_l2,"
                 << "exact_phi_diffusion_term_l2,"
                 << "solved_phi_diffusion_term_l2,"
                 << "exact_sigma_a_coupling_term_l2,"
                 << "solved_sigma_a_coupling_term_l2,"
                 << "exact_rhs_phi_l2,"
                 << "solved_rhs_phi_l2,"
                 << "solved_linear_residual_norm,"
                 << "solved_total_joule_power,"
                 << "solved_gauge_penalty_energy,"
                 << "solved_magnetic_energy,"
                 << "smooth_mode_laplace_energy,"
                 << "smooth_mode_weak_curlcurl_energy,"
                 << "smooth_mode_weak_to_laplace_ratio,"
                 << "checker_mode_laplace_energy,"
                 << "checker_mode_weak_curlcurl_energy,"
                 << "checker_mode_weak_to_laplace_ratio,"
                 << "interface_shell_thickness,"
                 << "interface_shell_particles,"
                 << "non_interface_particles,"
                 << "interface_shell_mean_a_error,"
                 << "interface_shell_max_a_error,"
                 << "interface_shell_mean_phi_error,"
                 << "interface_shell_max_phi_error,"
                 << "interface_shell_mean_exact_abs_j,"
                 << "interface_shell_mean_abs_j,"
                 << "interface_shell_max_abs_j,"
                 << "interface_shell_mean_abs_j_error,"
                 << "interface_shell_mean_rel_abs_j_error,"
                 << "interface_shell_mean_exact_joule,"
                 << "interface_shell_mean_joule,"
                 << "interface_shell_max_joule,"
                 << "interface_shell_mean_joule_error,"
                 << "interface_shell_mean_rel_joule_error,"
                 << "non_interface_mean_a_error,"
                 << "non_interface_max_a_error,"
                 << "non_interface_mean_phi_error,"
                 << "non_interface_max_phi_error,"
                 << "non_interface_mean_exact_abs_j,"
                 << "non_interface_mean_abs_j,"
                 << "non_interface_max_abs_j,"
                 << "non_interface_mean_abs_j_error,"
                 << "non_interface_mean_rel_abs_j_error,"
                 << "non_interface_mean_exact_joule,"
                 << "non_interface_mean_joule,"
                 << "non_interface_max_joule,"
                 << "non_interface_mean_joule_error,"
                 << "non_interface_mean_rel_joule_error\n";
    summary_file << (solve_success ? 1 : 0) << ","
                 << effective_rhs_mode_token << ","
                 << case_mode_token << ","
                 << total_particles << ","
                 << (enable_block_scaling ? 1 : 0) << ","
                 << solver.UsedPhiBlockScale() << ","
                 << (enable_diagonal_equilibration ? 1 : 0) << ","
                 << solver.UsedDiagonalEquilibrationIterations() << ","
                 << solver.LastSolverBackend() << ","
                 << solver.LastSolverIterations() << ","
                 << solver.LastSolverEstimatedError() << ","
                 << (enable_layered_material ? 1 : 0) << ","
                 << material_axis << ","
                 << layered_material_split_fraction << ","
                 << material_summary.lower_sigma << ","
                 << material_summary.upper_sigma << ","
                 << material_summary.lower_nu << ","
                 << material_summary.upper_nu << ","
                 << material_summary.lower_particles << ","
                 << material_summary.upper_particles << ","
                 << exact_linear_residual_norm << ","
                 << error_summary.mean_a_error << ","
                 << error_summary.max_a_error << ","
                 << error_summary.l2_a_error << ","
                 << error_summary.mean_phi_error << ","
                 << error_summary.max_phi_error << ","
                 << error_summary.l2_phi_error << ","
                 << error_summary.mean_abs_j_error << ","
                 << error_summary.max_abs_j_error << ","
                 << error_summary.l2_abs_j_error << ","
                 << error_summary.mean_rel_abs_j_error << ","
                 << error_summary.mean_joule_error << ","
                 << error_summary.max_joule_error << ","
                 << error_summary.l2_joule_error << ","
                 << error_summary.mean_rel_joule_error << ","
                 << exact_diagnostics.divergence_a_l2 << ","
                 << solved_diagnostics.divergence_a_l2 << ","
                 << exact_diagnostics.current_continuity_l2 << ","
                 << solved_diagnostics.current_continuity_l2 << ","
                 << exact_diagnostics.phi_equation_residual_l2 << ","
                 << solved_diagnostics.phi_equation_residual_l2 << ","
                 << exact_diagnostics.physical_current_continuity_l2 << ","
                 << solved_diagnostics.physical_current_continuity_l2 << ","
                 << exact_diagnostics.relative_current_continuity_l2 << ","
                 << solved_diagnostics.relative_current_continuity_l2 << ","
                 << exact_diagnostics.phi_diffusion_term_l2 << ","
                 << solved_diagnostics.phi_diffusion_term_l2 << ","
                 << exact_diagnostics.sigma_a_coupling_term_l2 << ","
                 << solved_diagnostics.sigma_a_coupling_term_l2 << ","
                 << exact_diagnostics.rhs_phi_l2 << ","
                 << solved_diagnostics.rhs_phi_l2 << ","
                 << solved_diagnostics.linear_residual_norm << ","
                 << solved_diagnostics.total_joule_power << ","
                 << solved_diagnostics.gauge_penalty_energy << ","
                 << solved_diagnostics.magnetic_energy << ","
                 << smooth_mode_summary.laplace_energy << ","
                 << smooth_mode_summary.weak_curlcurl_energy << ","
                 << smooth_weak_to_laplace_ratio << ","
                 << checker_mode_summary.laplace_energy << ","
                 << checker_mode_summary.weak_curlcurl_energy << ","
                 << checker_weak_to_laplace_ratio << ","
                 << interface_shell_thickness << ","
                 << interface_shell_summary.particles << ","
                 << non_interface_summary.particles << ","
                 << interface_shell_summary.mean_a_error << ","
                 << interface_shell_summary.max_a_error << ","
                 << interface_shell_summary.mean_phi_error << ","
                 << interface_shell_summary.max_phi_error << ","
                 << interface_shell_summary.mean_exact_abs_j << ","
                 << interface_shell_summary.mean_abs_j << ","
                 << interface_shell_summary.max_abs_j << ","
                 << interface_shell_summary.mean_abs_j_error << ","
                 << interface_shell_summary.mean_rel_abs_j_error << ","
                 << interface_shell_summary.mean_exact_joule << ","
                 << interface_shell_summary.mean_joule << ","
                 << interface_shell_summary.max_joule << ","
                 << interface_shell_summary.mean_joule_error << ","
                 << interface_shell_summary.mean_rel_joule_error << ","
                 << non_interface_summary.mean_a_error << ","
                 << non_interface_summary.max_a_error << ","
                 << non_interface_summary.mean_phi_error << ","
                 << non_interface_summary.max_phi_error << ","
                 << non_interface_summary.mean_exact_abs_j << ","
                 << non_interface_summary.mean_abs_j << ","
                 << non_interface_summary.max_abs_j << ","
                 << non_interface_summary.mean_abs_j_error << ","
                 << non_interface_summary.mean_rel_abs_j_error << ","
                 << non_interface_summary.mean_exact_joule << ","
                 << non_interface_summary.mean_joule << ","
                 << non_interface_summary.max_joule << ","
                 << non_interface_summary.mean_joule_error << ","
                 << non_interface_summary.mean_rel_joule_error << "\n";

    std::ofstream mode_file(
        io_environment.OutputFolder() + "/em_aphi_laplace_eigen_mode_comparison.csv",
        std::ios::out | std::ios::trunc);
    mode_file << std::setprecision(12);
    mode_file << "mode,laplace_energy,weak_curlcurl_energy,laplace_energy_norm,weak_curlcurl_energy_norm,"
              << "divergence_l2,curl_l2,max_abs_a,weak_to_laplace_ratio\n";
    auto write_mode_summary = [&](const std::string &mode_name,
                                  const electromagnetics::DiscreteModeEnergySummary &summary)
    {
        Real laplace_energy_norm = normalized_energy(summary, summary.laplace_energy);
        Real weak_energy_norm = normalized_energy(summary, summary.weak_curlcurl_energy);
        Real ratio = weak_energy_norm / (laplace_energy_norm + TinyReal);
        mode_file << mode_name << ","
                  << summary.laplace_energy << ","
                  << summary.weak_curlcurl_energy << ","
                  << laplace_energy_norm << ","
                  << weak_energy_norm << ","
                  << summary.divergence_l2 << ","
                  << summary.curl_l2 << ","
                  << summary.max_abs_a << ","
                  << ratio << "\n";
    };
    write_mode_summary("smooth_divfree", smooth_mode_summary);
    write_mode_summary("checkerboard_ax", checker_mode_summary);

    if (write_particle_csv)
    {
        std::ofstream particle_file(
            io_environment.OutputFolder() + "/em_aphi_laplace_eigen_particles.csv",
            std::ios::out | std::ios::trunc);
        particle_file << std::setprecision(12);
        particle_file << "particle_id,x,y,z,is_boundary,is_upper_material,is_interface_shell,interface_distance,sigma,nu,"
                      << "ax_exact_real,ax_exact_imag,ax_solved_real,ax_solved_imag,"
                      << "ay_exact_real,ay_exact_imag,ay_solved_real,ay_solved_imag,"
                      << "az_exact_real,az_exact_imag,az_solved_real,az_solved_imag,"
                      << "phi_exact_real,phi_exact_imag,phi_solved_real,phi_solved_imag,"
                      << "exact_abs_j,solved_abs_j,j_error_abs,j_error_relative,"
                      << "exact_joule,solved_joule,joule_error_abs,joule_error_relative,"
                      << "a_error_abs,phi_error_abs\n";

        for (size_t i = 0; i != total_particles; ++i)
        {
            Real a_error = std::sqrt(std::norm(exact_fields.ax[i] - solved_fields.ax[i]) +
                                     std::norm(exact_fields.ay[i] - solved_fields.ay[i]) +
                                     std::norm(exact_fields.az[i] - solved_fields.az[i]));
            Real phi_error = std::abs(exact_fields.phi[i] - solved_fields.phi[i]);
            Real exact_abs_j = exact_fields.current_density.empty() ? 0.0 :
                static_cast<Real>(exact_fields.current_density[i].norm());
            Real solved_abs_j = solved_fields.current_density.empty() ? 0.0 :
                static_cast<Real>(solved_fields.current_density[i].norm());
            Real j_error_abs = std::abs(solved_abs_j - exact_abs_j);
            Real j_error_relative = j_error_abs / (exact_abs_j + TinyReal);
            Real exact_joule = exact_fields.joule_heat.empty() ? 0.0 : exact_fields.joule_heat[i];
            Real solved_joule = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];
            Real joule_error_abs = std::abs(solved_joule - exact_joule);
            Real joule_error_relative = joule_error_abs / (std::abs(exact_joule) + TinyReal);
            bool upper_material = enable_layered_material &&
                                  is_upper_material(positions[i], material_axis, split_coordinate);
            Real particle_interface_distance = enable_layered_material ?
                interface_distance(positions[i], material_axis, split_coordinate) : -1.0;
            bool in_interface_shell = enable_layered_material &&
                                      is_interface_shell_particle(positions[i], material_axis,
                                                                  split_coordinate, interface_shell_thickness);
            particle_file << i << ","
                          << positions[i][0] << ","
                          << positions[i][1] << ","
                          << positions[i][2] << ","
                          << (boundary_condition.is_dirichlet[i] ? 1 : 0) << ","
                          << (upper_material ? 1 : 0) << ","
                          << (in_interface_shell ? 1 : 0) << ","
                          << particle_interface_distance << ","
                          << sigma[i] << ","
                          << nu[i] << ","
                          << exact_fields.ax[i].real() << ","
                          << exact_fields.ax[i].imag() << ","
                          << solved_fields.ax[i].real() << ","
                          << solved_fields.ax[i].imag() << ","
                          << exact_fields.ay[i].real() << ","
                          << exact_fields.ay[i].imag() << ","
                          << solved_fields.ay[i].real() << ","
                          << solved_fields.ay[i].imag() << ","
                          << exact_fields.az[i].real() << ","
                          << exact_fields.az[i].imag() << ","
                          << solved_fields.az[i].real() << ","
                          << solved_fields.az[i].imag() << ","
                          << exact_fields.phi[i].real() << ","
                          << exact_fields.phi[i].imag() << ","
                          << solved_fields.phi[i].real() << ","
                          << solved_fields.phi[i].imag() << ","
                          << exact_abs_j << ","
                          << solved_abs_j << ","
                          << j_error_abs << ","
                          << j_error_relative << ","
                          << exact_joule << ","
                          << solved_joule << ","
                          << joule_error_abs << ","
                          << joule_error_relative << ","
                          << a_error << ","
                          << phi_error << "\n";
        }
    }

    std::cout << "[em-aphi-laplace] " << solver_message
              << ", used_phi_block_scale=" << solver.UsedPhiBlockScale()
              << ", used_diagonal_equilibration_iterations="
              << solver.UsedDiagonalEquilibrationIterations()
              << ", solver_backend=" << solver.LastSolverBackend()
              << ", solver_iterations=" << solver.LastSolverIterations()
              << ", solver_estimated_error=" << solver.LastSolverEstimatedError()
              << ", enable_layered_material=" << enable_layered_material
              << ", lower_material_particles=" << material_summary.lower_particles
              << ", upper_material_particles=" << material_summary.upper_particles
              << ", interface_shell_particles=" << interface_shell_summary.particles
              << ", interface_shell_mean_a_error=" << interface_shell_summary.mean_a_error
              << ", exact_linear_residual_norm=" << exact_linear_residual_norm
              << ", mean_a_error=" << error_summary.mean_a_error
              << ", l2_a_error=" << error_summary.l2_a_error
              << ", max_a_error=" << error_summary.max_a_error
              << ", solved_phi_equation_residual_l2=" << solved_diagnostics.phi_equation_residual_l2
              << ", solved_physical_current_continuity_l2=" << solved_diagnostics.physical_current_continuity_l2
              << ", smooth_mode_ratio=" << smooth_weak_to_laplace_ratio
              << ", checker_mode_ratio=" << checker_weak_to_laplace_ratio
              << ", solved_linear_residual_norm=" << solved_diagnostics.linear_residual_norm
              << std::endl;

    return 0;
}
