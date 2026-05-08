/**
 * @file  particle_generation_em.cpp
 * @brief TEAM7-inspired single-way A-phi electromagnetic heating scaffold.
 * @details
 * 1) Two-stage workflow:
 *    - Stage A (--relax=1): generate particles, relax coil/plate, write reload files, then exit.
 *    - Stage B (--reload=1): load particles from reload and run A-phi + thermal solve.
 * 2) Particle setup: coil + plate with independent single-resolution spacing by default,
 *    air with adaptive multi-resolution around coil/plate. Optional TEAM7_USE_SMALL_AIR_BOX,
 *    TEAM7_AIR_BOX_* corners, TEAM7_UNIFORM_RESOLUTION, TEAM7_SPH_REF_DP (see parse_team7_mesh_params_pre_system).
 * 3) Electromagnetic-thermal solve (current stage): coil-plate one-way scaffold,
 *    with harmonic equivalent source drive and phi constraints as
 *    "reference point + specified boundary strip".
 *
 * NOTE:
 * Coil/plate receive inner+contact electromagnetic terms. The formulation is still
 * time-domain harmonic-equivalent (not full complex-frequency A-phi).
 */

#include "sphinxsys.h"
#include "em_adaptive_cell_linked_list.h"
#include "electromagnetic_component_hessian_ck.hpp"
#include "electromagnetic_multiturn_coil_drive.h"
#include "electromagnetic_team7_aphi_dynamics.hpp"
#include "electromagnetic_team7_aphi_frequency_dynamics.hpp"
#include "general_gradient.h"
#include "hessian_correction_ck.h"
#include "kernel_correction_ck.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace SPH;

namespace
{
//----------------------------------------------------------------------
//  STL paths
//----------------------------------------------------------------------
const std::string path_coil_stl = "./input/coil.stl";
const std::string path_plate_stl = "./input/plate.stl";

//----------------------------------------------------------------------
//  Domain and resolution defaults (overridden at runtime; see parse_team7_mesh_params_pre_system)
//----------------------------------------------------------------------
constexpr Real kDefaultSphRefDp = 2.8;
constexpr Real kDefaultDpCoil = 4.0;
constexpr Real kDefaultDpPlate = 6.0;
constexpr Real kDefaultDpAirFinest = 3.0;
constexpr int kDefaultAirRefinementLevels = 4;

struct Team7MeshParams
{
    Vec3d air_lower{-200.0, -200.0, -200.0};
    Vec3d air_upper{500.0, 500.0, 300.0};
    Real dp_ref = kDefaultSphRefDp;
    Real dp_coil = kDefaultDpCoil;
    Real dp_plate = kDefaultDpPlate;
    Real dp_air_finest = kDefaultDpAirFinest;
    int air_levels = kDefaultAirRefinementLevels;
    bool used_small_air_box_preset = false;
    bool uniform_resolution = false;
};

inline bool parse_bool_env_pre(const char *name, bool default_value)
{
    const char *env_value = std::getenv(name);
    if (env_value == nullptr)
    {
        return default_value;
    }
    std::string token(env_value);
    for (char &ch : token)
    {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
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

inline Real parse_real_env_pre(const char *name, Real default_value)
{
    const char *env_value = std::getenv(name);
    if (env_value == nullptr)
    {
        return default_value;
    }
    char *parse_end = nullptr;
    Real parsed = static_cast<Real>(std::strtod(env_value, &parse_end));
    if (parse_end == env_value || !std::isfinite(parsed) || parsed <= TinyReal)
    {
        return default_value;
    }
    return parsed;
}

inline Real parse_real_env_any_pre(const char *name, Real default_value)
{
    const char *env_value = std::getenv(name);
    if (env_value == nullptr)
    {
        return default_value;
    }
    char *parse_end = nullptr;
    Real parsed = static_cast<Real>(std::strtod(env_value, &parse_end));
    if (parse_end == env_value || !std::isfinite(parsed))
    {
        return default_value;
    }
    return parsed;
}

inline int parse_nonnegative_int_env_pre(const char *name, int default_value)
{
    const char *env_value = std::getenv(name);
    if (env_value == nullptr)
    {
        return default_value;
    }
    char *parse_end = nullptr;
    long parsed = std::strtol(env_value, &parse_end, 10);
    if (parse_end == env_value || parsed < 0L)
    {
        return default_value;
    }
    return static_cast<int>(parsed);
}

/**
 * Read air box corners and discretization before SPHSystem construction.
 * TEAM7_USE_SMALL_AIR_BOX: preset smaller air domain (~370x440x220 mm) to limit particle count;
 *   override any corner with TEAM7_AIR_BOX_LOWER_* / TEAM7_AIR_BOX_UPPER_* (mm).
 * TEAM7_UNIFORM_RESOLUTION: set dp_coil=dp_plate=dp_air_finest=TEAM7_SPH_REF_DP and air refinement 0
 *   (single air spacing; AdaptiveNearSurface reduces to uniform when finest==coarsest).
 */
Team7MeshParams parse_team7_mesh_params_pre_system()
{
    Team7MeshParams m;
    if (parse_bool_env_pre("TEAM7_USE_SMALL_AIR_BOX", true))
    {
        m.used_small_air_box_preset = true;
        m.air_lower = Vec3d(-30.0, -80.0, -40.0);
        m.air_upper = Vec3d(340.0, 360.0, 180.0);
    }
    m.air_lower[0] = parse_real_env_any_pre("TEAM7_AIR_BOX_LOWER_X", m.air_lower[0]);
    m.air_lower[1] = parse_real_env_any_pre("TEAM7_AIR_BOX_LOWER_Y", m.air_lower[1]);
    m.air_lower[2] = parse_real_env_any_pre("TEAM7_AIR_BOX_LOWER_Z", m.air_lower[2]);
    m.air_upper[0] = parse_real_env_any_pre("TEAM7_AIR_BOX_UPPER_X", m.air_upper[0]);
    m.air_upper[1] = parse_real_env_any_pre("TEAM7_AIR_BOX_UPPER_Y", m.air_upper[1]);
    m.air_upper[2] = parse_real_env_any_pre("TEAM7_AIR_BOX_UPPER_Z", m.air_upper[2]);

    m.dp_ref = parse_real_env_pre("TEAM7_SPH_REF_DP", kDefaultSphRefDp);
    if (parse_bool_env_pre("TEAM7_UNIFORM_RESOLUTION", true))
    {
        m.uniform_resolution = true;
        m.dp_coil = m.dp_plate = m.dp_air_finest = m.dp_ref;
        m.air_levels = 0;
    }
    else
    {
        m.dp_coil = parse_real_env_pre("TEAM7_DP_COIL", kDefaultDpCoil);
        m.dp_plate = parse_real_env_pre("TEAM7_DP_PLATE", kDefaultDpPlate);
        m.dp_air_finest = parse_real_env_pre("TEAM7_DP_AIR_FINEST", kDefaultDpAirFinest);
        m.air_levels =
            parse_nonnegative_int_env_pre("TEAM7_AIR_REFINEMENT_LEVELS", kDefaultAirRefinementLevels);
    }
    return m;
}

//----------------------------------------------------------------------
//  Electromagnetic / thermal scaffold parameters
//----------------------------------------------------------------------
const Real sigma_plate = 3.526e7;
const Real rho_cp_plate = 3.6e6;
const Real sigma_coil = 6.0e7;
const Real rho_cp_coil = 3.4e6;
// NOTE:
// For the current pseudo-time A-phi scaffold, near-zero air conductivity causes
// excessively stiff A updates (rhs / sigma). Use a regularized value for
// numerical stability before introducing a dedicated magnetic-only air solver.
const Real sigma_air = 1.0;
const Real rho_cp_air = 1.0;
const Real magnetic_reluctivity = 1.0 / (4.0 * Pi * 1.0e-7);

const Vec3d harmonic_source_amplitude(0.0, 0.0, 8.0e5);
const Real harmonic_frequency_hz = 50.0;
const Real harmonic_phase = 0.0;
const Real harmonic_angular_frequency = 2.0 * Pi * harmonic_frequency_hz;

const Real thermal_diffusivity = 5.0e-6;
const Real initial_temperature = 293.15;

const Real end_time_default = 2.0e-3;
const size_t output_steps_default = 20;
const Real dt_thermal_max = 1.0e-4;
const Real dt_thermal_min = 1.0e-6;
const Real dt_em_default = 1.0e-6;
const size_t em_substeps_per_thermal_default = 1;
const size_t phi_relax_iterations_default = 8;
const size_t a_relax_iterations_default = 1;
const Real coil_a_relaxation_scaling_default = 1.0;
const Real plate_a_relaxation_scaling_default = 1.0;
const Real air_a_relaxation_scaling_default = 1.0;
const Real coil_a_rate_limit_default = 1.0e4;
const Real plate_a_rate_limit_default = 1.0e4;
const Real air_a_rate_limit_default = 1.0e4;
const Real phi_abs_limit_default = 1.0e2;
const int particle_relaxation_steps = 200;
const int particle_relaxation_output_interval = 100;
const Real particle_randomization_factor = 0.25;
const bool debug_coil_nan_detection = false;

//----------------------------------------------------------------------
//  Shapes
//----------------------------------------------------------------------
class CoilShape : public ComplexShape
{
  public:
    explicit CoilShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0);
    }
};

class PlateShape : public ComplexShape
{
  public:
    explicit PlateShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0);
    }
};

class AirShape : public ComplexShape
{
  public:
    explicit AirShape(const std::string &shape_name, const Vec3d &box_lower, const Vec3d &box_upper)
        : ComplexShape(shape_name)
    {
        Vecd halfsize = 0.5 * (box_upper - box_lower);
        Vecd center = 0.5 * (box_lower + box_upper);
        Transform translation(center);
        add<GeometricShapeBox>(translation, halfsize, "OuterBoundary");
        subtract<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
        subtract<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
    }
};

class InnerBoundaryShape : public ComplexShape
{
  public:
    explicit InnerBoundaryShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
    }
};

class AirOuterBoundaryShell : public ComplexShape
{
  public:
    explicit AirOuterBoundaryShell(const std::string &shape_name,
                                   const Vec3d &box_lower,
                                   const Vec3d &box_upper,
                                   Real shell_thickness)
        : ComplexShape(shape_name)
    {
        Vecd outer_halfsize = 0.5 * (box_upper - box_lower);
        Vecd center = 0.5 * (box_lower + box_upper);
        Vecd inner_halfsize = outer_halfsize - Vecd::Ones() * shell_thickness;
        add<GeometricShapeBox>(Transform(center), outer_halfsize, "AirOuterBoundary");
        if (inner_halfsize.minCoeff() > TinyReal)
        {
            subtract<GeometricShapeBox>(Transform(center), inner_halfsize, "AirInnerCore");
        }
    }
};

class AdaptiveNearInnerSurface
    : public extra_em_mesh::AdaptiveWithTolerantCellLinkedList<AdaptiveNearSurface>
{
    Shape *inner_shape_;

  public:
    AdaptiveNearInnerSurface(Real global_resolution, Real h_spacing_ratio,
                             Real refinement_to_global, int local_refinement_level,
                             Shape *inner_shape)
        : extra_em_mesh::AdaptiveWithTolerantCellLinkedList<AdaptiveNearSurface>(
              global_resolution, h_spacing_ratio,
              refinement_to_global, local_refinement_level),
          inner_shape_(inner_shape) {}

    Real getLocalSpacing(Shape &shape, const Vecd &position) override
    {
        Real phi = fabs(inner_shape_->findSignedDistance(position));
        return smoothedSpacing(phi, spacing_ref_);
    }
};

//----------------------------------------------------------------------
//  Thermal source/initialization helpers
//----------------------------------------------------------------------
class TemperatureInitialization : public LocalDynamics
{
  public:
    explicit TemperatureInitialization(SPHBody &sph_body, Real initial_temperature)
        : LocalDynamics(sph_body),
          initial_temperature_(initial_temperature),
          temperature_(particles_->registerStateVariableData<Real>("Temperature"))
    {
        particles_->addVariableToWrite<Real>("Temperature");
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        temperature_[index_i] = initial_temperature_;
    }

  protected:
    Real initial_temperature_;
    Real *temperature_;
};

class AddJouleHeatToTemperature : public LocalDynamics
{
  public:
    explicit AddJouleHeatToTemperature(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          temperature_(particles_->getVariableDataByName<Real>("Temperature")),
          temperature_change_rate_by_joule_(particles_->getVariableDataByName<Real>("TemperatureChangeRateByJoule")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        temperature_[index_i] += dt * temperature_change_rate_by_joule_[index_i];
    }

  protected:
    Real *temperature_;
    Real *temperature_change_rate_by_joule_;
};

//----------------------------------------------------------------------
//  Harmonic equivalent source: J_s(t) = J0 * sin(2*pi*f*t + phase)
//----------------------------------------------------------------------
class HarmonicSourceCurrentDensity : public LocalDynamics
{
  public:
    explicit HarmonicSourceCurrentDensity(SPHBody &sph_body,
                                          const Vecd &amplitude,
                                          Real frequency_hz,
                                          Real &physical_time,
                                          Real phase = 0.0)
        : LocalDynamics(sph_body),
          amplitude_(amplitude),
          frequency_hz_(frequency_hz),
          phase_(phase),
          physical_time_(physical_time),
          source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Real harmonic_factor = sin(2.0 * Pi * frequency_hz_ * physical_time_ + phase_);
        source_current_density_[index_i] = amplitude_ * harmonic_factor;
    }

  protected:
    Vecd amplitude_;
    Real frequency_hz_;
    Real phase_;
    Real &physical_time_;
    Vecd *source_current_density_;
};

class ScaleSourceCurrentDensityByRuntimeFactor : public LocalDynamics
{
  public:
    explicit ScaleSourceCurrentDensityByRuntimeFactor(SPHBody &sph_body, Real &source_scale)
        : LocalDynamics(sph_body),
          source_scale_(source_scale),
          source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        source_current_density_[index_i] *= source_scale_;
    }

  protected:
    Real &source_scale_;
    Vecd *source_current_density_;
};

class ScaleComplexSourceCurrentDensityByRuntimeFactor : public LocalDynamics
{
  public:
    explicit ScaleComplexSourceCurrentDensityByRuntimeFactor(SPHBody &sph_body, Real &source_scale)
        : LocalDynamics(sph_body),
          source_scale_(source_scale),
          source_current_density_real_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityReal")),
          source_current_density_imag_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityImag")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        (void)dt;
        source_current_density_real_[index_i] *= source_scale_;
        source_current_density_imag_[index_i] *= source_scale_;
    }

  protected:
    Real &source_scale_;
    Vecd *source_current_density_real_;
    Vecd *source_current_density_imag_;
};

using ThermalRelaxationInner =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;

struct PlateDiagnostics
{
    Real total_volume = 0.0;
    Real total_joule_power = 0.0;
    Real avg_joule_heat = 0.0;
    Real max_joule_heat = 0.0;
    Real avg_current_density_magnitude = 0.0;
    Real max_current_density_magnitude = 0.0;
    Real avg_temperature = 0.0;
    Real max_temperature = 0.0;
    Real max_a_rate = 0.0;
    Real max_a_dot = 0.0;
    Real max_grad_phi = 0.0;
    Real max_electric_field = 0.0;
};

struct EmRateBreakdown
{
    Real plate = 0.0;
    Real coil = 0.0;
    Real air = 0.0;
};

struct EmResidualBreakdown
{
    Real plate = 0.0;
    Real coil = 0.0;
    Real air = 0.0;
};

struct EmBodyTermDiagnostics
{
    Real total_volume = 0.0;
    Real avg_source = 0.0;
    Real avg_curl_nu_b = 0.0;
    Real avg_sigma_grad_phi = 0.0;
    Real avg_omega_sigma_a = 0.0;
    Real avg_residual = 0.0;
    Real max_source = 0.0;
    Real max_curl_nu_b = 0.0;
    Real max_sigma_grad_phi = 0.0;
    Real max_omega_sigma_a = 0.0;
    Real max_residual = 0.0;
};

struct EmTermBreakdown
{
    EmBodyTermDiagnostics plate;
    EmBodyTermDiagnostics coil;
    EmBodyTermDiagnostics air;
};

struct EmInterfaceShellBreakdown
{
    EmBodyTermDiagnostics plate_inner;
    EmBodyTermDiagnostics plate_air;
    EmBodyTermDiagnostics coil_inner;
    EmBodyTermDiagnostics coil_air;
};

/** Volume-weighted |A| and |B| on coil inner shell vs air shell near coil (frequency A-phi). */
struct EmCoilAirInterfaceShellFieldDiagnostics
{
    Real coil_inner_volume = 0.0;
    Real coil_inner_avg_a_mag = 0.0;
    Real coil_inner_avg_b_mag = 0.0;
    Real coil_air_volume = 0.0;
    Real coil_air_avg_a_mag = 0.0;
    Real coil_air_avg_b_mag = 0.0;
};

struct EmMagneticDiagonalBodyDiagnostics
{
    Real total_volume = 0.0;
    Real avg_local_diagonal = 0.0;
    Real avg_contact_diagonal = 0.0;
    Real avg_jacobi_diagonal = 0.0;
    Real avg_conservative_diagonal = 0.0;
    Real avg_balanced_diagonal = 0.0;
    Real avg_contact_ratio = 0.0;
    Real avg_conservative_to_jacobi = 0.0;
    Real avg_balanced_to_jacobi = 0.0;
    Real max_contact_ratio = 0.0;
};

struct EmMagneticDiagonalBreakdown
{
    EmMagneticDiagonalBodyDiagnostics plate;
    EmMagneticDiagonalBodyDiagnostics coil;
    EmMagneticDiagonalBodyDiagnostics air;
    EmMagneticDiagonalBodyDiagnostics plate_inner;
    EmMagneticDiagonalBodyDiagnostics plate_air;
    EmMagneticDiagonalBodyDiagnostics coil_inner;
    EmMagneticDiagonalBodyDiagnostics coil_air;
};

struct FrequencyPlateAcceptanceDiagnostics
{
    Real body_avg_residual_ratio = 0.0;
    Real inner_shell_avg_residual_ratio = 0.0;
    Real air_shell_avg_residual_ratio = 0.0;
};

struct FrequencyEmAcceptanceDiagnostics
{
    Real conductor_metric = 0.0;
    Real plate_shell_metric = 0.0;
    Real air_metric = 0.0;
    Real weighted_metric = 0.0;
};

struct FrequencyEmHardGuardDiagnostics
{
    Real weighted_metric = 0.0;
    Real probe_b_metric = 0.0;
    Real probe_ref_ratio_metric = 0.0;
};

struct FrequencyEmStateBackup
{
    std::vector<Vecd> plate_a_real;
    std::vector<Vecd> plate_a_imag;
    std::vector<Real> plate_phi_real;
    std::vector<Real> plate_phi_imag;
    std::vector<Vecd> coil_a_real;
    std::vector<Vecd> coil_a_imag;
    std::vector<Real> coil_phi_real;
    std::vector<Real> coil_phi_imag;
    std::vector<Vecd> air_a_real;
    std::vector<Vecd> air_a_imag;
    std::vector<Real> air_phi_real;
    std::vector<Real> air_phi_imag;
};

struct CoilSourceDiagnostics
{
    Real total_volume_si = 0.0;
    Real path_length_si = 0.0;
    Real configured_ampere_turns_real = 0.0;
    Real configured_ampere_turns_imag = 0.0;
    Real equivalent_ampere_turns_real = 0.0;
    Real equivalent_ampere_turns_imag = 0.0;
    Real equivalent_ampere_turns_magnitude = 0.0;
    Real avg_source_density_magnitude = 0.0;
    Real max_source_density_magnitude = 0.0;
};

struct ProbeDescriptor
{
    std::string name;
    size_t particle_index = 0;
    Vec3d target_position = Vec3d::Zero();
    Vec3d sampled_position = Vec3d::Zero();
};

struct Team7CurveReferencePoint
{
    Real x_mm = 0.0;
    Real ref_50 = std::numeric_limits<Real>::quiet_NaN();
    Real ref_200 = std::numeric_limits<Real>::quiet_NaN();
};

struct Team7CurveSamplePoint
{
    std::string quantity_name;
    size_t sample_index = 0;
    Real x_mm = 0.0;
    Vec3d target_position = Vec3d::Zero();
    Vec3d sampled_position = Vec3d::Zero();
    Real simulated_value = std::numeric_limits<Real>::quiet_NaN();
    Real reference_50 = std::numeric_limits<Real>::quiet_NaN();
    Real reference_200 = std::numeric_limits<Real>::quiet_NaN();
    Real reference_selected = std::numeric_limits<Real>::quiet_NaN();
    Real absolute_error = std::numeric_limits<Real>::quiet_NaN();
    Real relative_error = std::numeric_limits<Real>::quiet_NaN();
    size_t support_point_count = 0;
    Real support_weight_sum = 0.0;
    Real support_nearest_distance = std::numeric_limits<Real>::quiet_NaN();
    Real support_weighted_mean_distance = std::numeric_limits<Real>::quiet_NaN();
    Real boundary_clearance_x = std::numeric_limits<Real>::quiet_NaN();
    Real boundary_clearance_y = std::numeric_limits<Real>::quiet_NaN();
    Real boundary_clearance_z = std::numeric_limits<Real>::quiet_NaN();
    bool high_quality = false;
    bool pass = false;
};

struct Team7CurveValidationSummary
{
    std::string quantity_name;
    std::string complex_mode_name;
    size_t total_points = 0;
    size_t checked_points = 0;
    size_t passed_points = 0;
    Real mean_abs_error = 0.0;
    Real rmse = 0.0;
    Real max_abs_error = 0.0;
    Real mean_rel_error = 0.0;
    /** Mean |reference| and |simulated| over checked points (same set as mean_abs_error). */
    Real mean_abs_reference = 0.0;
    Real mean_abs_simulated = 0.0;
    /** mean_abs_reference / (mean_abs_simulated + TinyReal); scale check vs COMSOL curve. */
    Real mean_abs_ref_over_mean_abs_sim = std::numeric_limits<Real>::quiet_NaN();
};

struct CurveInterpolationSupportInfo
{
    size_t support_point_count = 0;
    Real weight_sum = 0.0;
    Real nearest_distance = std::numeric_limits<Real>::quiet_NaN();
    Real weighted_mean_distance = std::numeric_limits<Real>::quiet_NaN();
};

struct ScalarValidationSummary
{
    std::string quantity_name;
    size_t checked_points = 0;
    Real mean_abs_error = 0.0;
    Real rmse = 0.0;
    Real max_abs_error = 0.0;
    Real mean_rel_error = 0.0;
};

struct OperatorVerificationRow
{
    std::string test_name;
    std::string quantity_name;
    size_t checked_particles = 0;
    size_t invalid_particles = 0;
    Real mean_abs_error = 0.0;
    Real max_abs_error = 0.0;
    Real mean_rel_error = 0.0;
    Real mean_projection_ratio = std::numeric_limits<Real>::quiet_NaN();
    Real mean_alignment = std::numeric_limits<Real>::quiet_NaN();
};

std::string trim_copy(const std::string &value)
{
    size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return "";
    }
    size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

Real parse_numeric_expression(const std::string &token)
{
    std::string expression = trim_copy(token);
    if (expression.empty())
    {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    expression.erase(
        std::remove_if(expression.begin(), expression.end(),
                       [](char ch) { return ch == ',' || std::isspace(static_cast<unsigned char>(ch)); }),
        expression.end());

    if (!expression.empty() &&
        (expression.back() == ';' || expression.back() == ','))
    {
        expression.pop_back();
    }

    size_t mul_pos = expression.find('*');
    if (mul_pos != std::string::npos)
    {
        Real product = 1.0;
        size_t start = 0;
        while (start < expression.size())
        {
            size_t sep = expression.find('*', start);
            std::string factor = expression.substr(start, sep == std::string::npos ? std::string::npos : sep - start);
            Real factor_value = parse_numeric_expression(factor);
            if (!std::isfinite(factor_value))
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
            product *= factor_value;
            if (sep == std::string::npos)
            {
                break;
            }
            start = sep + 1;
        }
        return product;
    }

    size_t pow_pos = expression.find('^');
    if (pow_pos != std::string::npos)
    {
        Real base = parse_numeric_expression(expression.substr(0, pow_pos));
        Real exponent = parse_numeric_expression(expression.substr(pow_pos + 1));
        if (!std::isfinite(base) || !std::isfinite(exponent))
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        return std::pow(base, exponent);
    }

    char *parse_end = nullptr;
    Real parsed_value = static_cast<Real>(std::strtod(expression.c_str(), &parse_end));
    if (parse_end == expression.c_str() || *parse_end != '\0' || !std::isfinite(parsed_value))
    {
        return std::numeric_limits<Real>::quiet_NaN();
    }
    return parsed_value;
}

bool load_team7_curve_reference_table(const std::string &table_file_path,
                                      std::vector<Team7CurveReferencePoint> &reference_points)
{
    reference_points.clear();
    std::ifstream table_file(table_file_path);
    if (!table_file.is_open())
    {
        return false;
    }

    std::string line;
    while (std::getline(table_file, line))
    {
        std::string trimmed = trim_copy(line);
        if (trimmed.empty() || trimmed[0] == '%')
        {
            continue;
        }

        std::istringstream line_stream(trimmed);
        std::string token_x, token_50, token_200;
        if (!(line_stream >> token_x >> token_50 >> token_200))
        {
            continue;
        }

        Team7CurveReferencePoint point;
        point.x_mm = parse_numeric_expression(token_x);
        point.ref_50 = parse_numeric_expression(token_50);
        point.ref_200 = parse_numeric_expression(token_200);
        if (std::isfinite(point.x_mm) &&
            (std::isfinite(point.ref_50) || std::isfinite(point.ref_200)))
        {
            reference_points.push_back(point);
        }
    }
    return !reference_points.empty();
}

Real select_reference_value_by_frequency(const Team7CurveReferencePoint &reference_point, Real reference_frequency_hz)
{
    bool has_ref_50 = std::isfinite(reference_point.ref_50);
    bool has_ref_200 = std::isfinite(reference_point.ref_200);
    if (!has_ref_50 && !has_ref_200)
    {
        return std::numeric_limits<Real>::quiet_NaN();
    }
    if (!has_ref_50)
    {
        return reference_point.ref_200;
    }
    if (!has_ref_200)
    {
        return reference_point.ref_50;
    }
    return fabs(reference_frequency_hz - 50.0) <= fabs(reference_frequency_hz - 200.0)
               ? reference_point.ref_50
               : reference_point.ref_200;
}

enum class ComplexCurveComponentMode
{
    RealPart,
    ImagPart,
    Magnitude,
    SignedMagnitude
};

ComplexCurveComponentMode parse_complex_curve_component_mode(
    const std::string &mode_token,
    ComplexCurveComponentMode default_mode)
{
    std::string mode = trim_copy(mode_token);
    if (mode.empty())
    {
        return default_mode;
    }
    std::transform(mode.begin(), mode.end(), mode.begin(), [](char ch)
                   { return static_cast<char>(std::tolower(static_cast<unsigned char>(ch))); });
    if (mode == "real" || mode == "real_part" || mode == "re")
    {
        return ComplexCurveComponentMode::RealPart;
    }
    if (mode == "imag" || mode == "imag_part" || mode == "im")
    {
        return ComplexCurveComponentMode::ImagPart;
    }
    if (mode == "magnitude" || mode == "mag" || mode == "abs")
    {
        return ComplexCurveComponentMode::Magnitude;
    }
    if (mode == "signed_magnitude" || mode == "signed-mag" || mode == "signedmag")
    {
        return ComplexCurveComponentMode::SignedMagnitude;
    }
    return default_mode;
}

const char *complex_curve_component_mode_name(ComplexCurveComponentMode mode)
{
    switch (mode)
    {
    case ComplexCurveComponentMode::RealPart:
        return "real";
    case ComplexCurveComponentMode::ImagPart:
        return "imag";
    case ComplexCurveComponentMode::Magnitude:
        return "magnitude";
    case ComplexCurveComponentMode::SignedMagnitude:
        return "signed_magnitude";
    default:
        return "signed_magnitude";
    }
}

Real signed_complex_component_magnitude(Real real_component, Real imag_component)
{
    Real magnitude = sqrt(real_component * real_component + imag_component * imag_component);
    return real_component >= 0.0 ? magnitude : -magnitude;
}

Real project_complex_component(Real real_component,
                               Real imag_component,
                               ComplexCurveComponentMode mode)
{
    switch (mode)
    {
    case ComplexCurveComponentMode::RealPart:
        return real_component;
    case ComplexCurveComponentMode::ImagPart:
        return imag_component;
    case ComplexCurveComponentMode::Magnitude:
        return sqrt(real_component * real_component + imag_component * imag_component);
    case ComplexCurveComponentMode::SignedMagnitude:
        return signed_complex_component_magnitude(real_component, imag_component);
    default:
        return signed_complex_component_magnitude(real_component, imag_component);
    }
}

} // namespace

int main(int ac, char *av[])
{
    const Team7MeshParams team7_mesh = parse_team7_mesh_params_pre_system();
    const Vec3d air_box_lower = team7_mesh.air_lower;
    const Vec3d air_box_upper = team7_mesh.air_upper;
    const BoundingBoxd system_domain_bounds(air_box_lower, air_box_upper);
    const Real dp_0 = team7_mesh.dp_ref;
    const Real dp_coil = team7_mesh.dp_coil;
    const Real dp_plate = team7_mesh.dp_plate;
    const Real dp_air_finest = team7_mesh.dp_air_finest;
    const int air_refinement_levels = team7_mesh.air_levels;
    const Real dp_air_coarsest =
        dp_air_finest * std::pow(2.0, static_cast<Real>(air_refinement_levels));

    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    auto get_env_real = [](const char *name, Real default_value) -> Real
    {
        const char *env_value = std::getenv(name);
        if (env_value == nullptr)
        {
            return default_value;
        }
        char *parse_end = nullptr;
        Real parsed_value = static_cast<Real>(std::strtod(env_value, &parse_end));
        if (parse_end == env_value || !std::isfinite(parsed_value) || parsed_value <= TinyReal)
        {
            return default_value;
        }
        return parsed_value;
    };
    auto get_env_real_any = [](const char *name, Real default_value) -> Real
    {
        const char *env_value = std::getenv(name);
        if (env_value == nullptr)
        {
            return default_value;
        }
        char *parse_end = nullptr;
        Real parsed_value = static_cast<Real>(std::strtod(env_value, &parse_end));
        if (parse_end == env_value || !std::isfinite(parsed_value))
        {
            return default_value;
        }
        return parsed_value;
    };
    auto get_env_size_t = [](const char *name, size_t default_value) -> size_t
    {
        const char *env_value = std::getenv(name);
        if (env_value == nullptr)
        {
            return default_value;
        }
        char *parse_end = nullptr;
        unsigned long long parsed_value = std::strtoull(env_value, &parse_end, 10);
        if (parse_end == env_value || parsed_value == 0ULL)
        {
            return default_value;
        }
        return static_cast<size_t>(parsed_value);
    };
    auto get_env_bool = [](const char *name, bool default_value) -> bool
    {
        const char *env_value = std::getenv(name);
        if (env_value == nullptr)
        {
            return default_value;
        }
        std::string parsed_value(env_value);
        for (char &ch : parsed_value)
        {
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
        if (parsed_value == "1" || parsed_value == "true" || parsed_value == "yes" || parsed_value == "on")
        {
            return true;
        }
        if (parsed_value == "0" || parsed_value == "false" || parsed_value == "no" || parsed_value == "off")
        {
            return false;
        }
        return default_value;
    };
    auto get_env_vec3 = [&](const std::string &prefix, const Vec3d &default_value) -> Vec3d
    {
        std::string name_x = prefix + "_X";
        std::string name_y = prefix + "_Y";
        std::string name_z = prefix + "_Z";
        return Vec3d(get_env_real_any(name_x.c_str(), default_value[0]),
                     get_env_real_any(name_y.c_str(), default_value[1]),
                     get_env_real_any(name_z.c_str(), default_value[2]));
    };
    auto get_env_string = [](const char *name, const std::string &default_value) -> std::string
    {
        const char *env_value = std::getenv(name);
        if (env_value == nullptr)
        {
            return default_value;
        }
        std::string parsed_value(env_value);
        if (parsed_value.empty())
        {
            return default_value;
        }
        return parsed_value;
    };
    auto parse_vector_component_index =
        [](const std::string &token, int default_index) -> int
    {
        std::string normalized = trim_copy(token);
        std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](char ch)
                       { return static_cast<char>(std::tolower(static_cast<unsigned char>(ch))); });
        if (normalized == "x" || normalized == "0")
        {
            return 0;
        }
        if (normalized == "y" || normalized == "1")
        {
            return 1;
        }
        if (normalized == "z" || normalized == "2")
        {
            return 2;
        }
        return default_index;
    };

    Real dt_em = get_env_real("TEAM7_DT_EM", dt_em_default);
    size_t em_substeps_per_thermal = get_env_size_t("TEAM7_EM_SUBSTEPS", em_substeps_per_thermal_default);
    size_t phi_relax_iterations = get_env_size_t("TEAM7_PHI_RELAX_ITER", phi_relax_iterations_default);
    size_t a_relax_iterations = get_env_size_t("TEAM7_A_RELAX_ITER", a_relax_iterations_default);
    Real coil_a_relaxation_scaling =
        get_env_real("TEAM7_COIL_A_RELAX_SCALE", coil_a_relaxation_scaling_default);
    Real plate_a_relaxation_scaling =
        get_env_real("TEAM7_PLATE_A_RELAX_SCALE", plate_a_relaxation_scaling_default);
    Real air_a_relaxation_scaling =
        get_env_real("TEAM7_AIR_A_RELAX_SCALE", air_a_relaxation_scaling_default);
    Real coil_a_rate_limit = get_env_real("TEAM7_COIL_A_RATE_LIMIT", coil_a_rate_limit_default);
    Real plate_a_rate_limit = get_env_real("TEAM7_PLATE_A_RATE_LIMIT", plate_a_rate_limit_default);
    Real air_a_rate_limit = get_env_real("TEAM7_AIR_A_RATE_LIMIT", air_a_rate_limit_default);
    Real phi_abs_limit = get_env_real("TEAM7_PHI_ABS_LIMIT", phi_abs_limit_default);
    Real freq_sigma_relaxation_scaling =
        get_env_real("TEAM7_FREQ_SIGMA_RELAX_SCALE", 1.0);
    Real freq_sigma_relaxation_floor =
        get_env_real("TEAM7_FREQ_SIGMA_RELAX_FLOOR", 1.0);
    freq_sigma_relaxation_scaling = SMAX(freq_sigma_relaxation_scaling, TinyReal);
    freq_sigma_relaxation_floor = SMAX(freq_sigma_relaxation_floor, TinyReal);
    Real harmonic_frequency_hz_runtime =
        get_env_real("TEAM7_HARMONIC_FREQUENCY_HZ", harmonic_frequency_hz);
    Real harmonic_phase_runtime =
        get_env_real_any("TEAM7_HARMONIC_PHASE", harmonic_phase);
    Real harmonic_angular_frequency_runtime =
        2.0 * Pi * harmonic_frequency_hz_runtime;
    bool use_frequency_aphi = get_env_bool("TEAM7_USE_FREQ_APHI", true);
    Real frequency_source_imag_scale = get_env_real_any("TEAM7_FREQ_SOURCE_IMAG_SCALE", 0.0);
    bool use_frequency_coupled_implicit_solver =
        get_env_bool("TEAM7_FREQ_COUPLED_IMPLICIT_SOLVER", true);
    bool use_frequency_pseudo_steady = get_env_bool("TEAM7_FREQ_PSEUDO_STEADY", true);
    bool frequency_em_solve_once = get_env_bool("TEAM7_FREQ_EM_SOLVE_ONCE", true);
    size_t frequency_pseudo_min_iterations =
        get_env_size_t("TEAM7_FREQ_PSEUDO_MIN_ITERS", 40);
    size_t frequency_pseudo_max_iterations =
        get_env_size_t("TEAM7_FREQ_PSEUDO_MAX_ITERS", 80);
    size_t frequency_block_gs_sweeps =
        get_env_size_t("TEAM7_FREQ_BLOCK_GS_SWEEPS", 3);
    size_t frequency_air_block_sweeps =
        get_env_size_t("TEAM7_FREQ_AIR_BLOCK_SWEEPS", use_frequency_aphi ? 2 : 1);
    size_t frequency_plate_block_sweeps =
        get_env_size_t("TEAM7_FREQ_PLATE_BLOCK_SWEEPS", 2);
    bool use_frequency_adaptive_plate_block_sweeps =
        get_env_bool("TEAM7_FREQ_ADAPTIVE_PLATE_BLOCK_SWEEPS", false);
    Real em_rate_convergence_tolerance =
        get_env_real("TEAM7_EM_RATE_CONVERGENCE_TOL", 10.0);
    Real em_rate_divergence_threshold =
        get_env_real("TEAM7_EM_RATE_DIVERGENCE_THRESHOLD", 1.0e6);
    bool use_enhanced_em_convergence =
        get_env_bool("TEAM7_EM_USE_ENHANCED_CONVERGENCE", true);
    Real em_delta_a_convergence_tolerance =
        get_env_real_any("TEAM7_EM_DELTA_A_CONVERGENCE_TOL", 1.0e-4);
    Real em_joule_rel_convergence_tolerance =
        get_env_real_any("TEAM7_EM_JOULE_REL_CONVERGENCE_TOL", 1.0e-4);
    bool em_stop_on_saturation =
        get_env_bool("TEAM7_EM_STOP_ON_SATURATION", true);
    Real em_saturation_ratio =
        get_env_real_any("TEAM7_EM_SATURATION_RATIO", 0.98);
    size_t em_saturation_patience =
        get_env_size_t("TEAM7_EM_SATURATION_PATIENCE", frequency_pseudo_max_iterations);
    bool use_em_adaptive_dt = get_env_bool("TEAM7_EM_ADAPTIVE_DT", true);
    Real dt_em_min = get_env_real("TEAM7_DT_EM_MIN", 1.0e-9);
    Real dt_em_max = get_env_real("TEAM7_DT_EM_MAX", dt_em);
    Real em_dt_reduce_factor =
        get_env_real_any("TEAM7_EM_DT_REDUCE_FACTOR", 0.5);
    Real em_dt_increase_factor =
        get_env_real_any("TEAM7_EM_DT_INCREASE_FACTOR", 1.02);
    size_t em_dt_instability_patience =
        get_env_size_t("TEAM7_EM_DT_INSTABILITY_PATIENCE", use_frequency_aphi ? 3 : 1);
    size_t em_dt_max_reductions_per_solve =
        get_env_size_t("TEAM7_EM_DT_MAX_REDUCTIONS_PER_SOLVE", use_frequency_aphi ? 64 : 1000000);
    Real em_dt_recovery_factor =
        get_env_real_any("TEAM7_EM_DT_RECOVERY_FACTOR", use_frequency_aphi ? 1.10 : em_dt_increase_factor);
    Real em_dt_recovery_min_factor =
        get_env_real_any("TEAM7_EM_DT_RECOVERY_MIN_FACTOR", use_frequency_aphi ? 8.0 : 1.0);
    Real em_dt_recovery_relative_threshold =
        get_env_real_any("TEAM7_EM_DT_RECOVERY_RELATIVE_THRESHOLD",
                         use_frequency_aphi ? 0.10 : 0.0);
    bool use_em_relative_residual =
        get_env_bool("TEAM7_EM_USE_RELATIVE_RESIDUAL", use_frequency_aphi);
    Real em_relative_residual_convergence_tolerance =
        get_env_real_any("TEAM7_EM_RELATIVE_RESIDUAL_TOL", 1.0e-2);
    Real em_rate_growth_trigger =
        get_env_real_any("TEAM7_EM_RATE_GROWTH_TRIGGER", use_frequency_aphi ? 1.5 : std::numeric_limits<Real>::infinity());
    Real em_residual_growth_trigger =
        get_env_real_any("TEAM7_EM_RESIDUAL_GROWTH_TRIGGER", use_frequency_aphi ? 1.01 : std::numeric_limits<Real>::infinity());
    Real em_rate_stability_limit_factor =
        get_env_real_any("TEAM7_EM_RATE_STABILITY_LIMIT_FACTOR", use_frequency_aphi ? 10.0 : std::numeric_limits<Real>::infinity());
    size_t em_adaptive_dt_warmup_iterations =
        get_env_size_t("TEAM7_EM_ADAPTIVE_DT_WARMUP_ITERS", use_frequency_aphi ? 10 : 0);
    dt_em_min = SMAX(dt_em_min, TinyReal);
    dt_em_max = SMAX(dt_em_max, dt_em_min);
    if (!(em_dt_reduce_factor > 0.0 && em_dt_reduce_factor < 1.0))
    {
        em_dt_reduce_factor = 0.5;
    }
    if (em_dt_increase_factor <= 1.0)
    {
        em_dt_increase_factor = 1.02;
    }
    em_dt_instability_patience =
        SMAX(static_cast<size_t>(1), em_dt_instability_patience);
    em_dt_max_reductions_per_solve =
        SMAX(static_cast<size_t>(1), em_dt_max_reductions_per_solve);
    if (em_dt_recovery_factor <= 1.0)
    {
        em_dt_recovery_factor = em_dt_increase_factor;
    }
    em_dt_recovery_factor = SMAX(em_dt_recovery_factor, em_dt_increase_factor);
    em_dt_recovery_min_factor =
        SMAX(static_cast<Real>(1.0), em_dt_recovery_min_factor);
    em_dt_recovery_relative_threshold =
        SMAX(static_cast<Real>(0.0), em_dt_recovery_relative_threshold);
    if (frequency_pseudo_min_iterations > frequency_pseudo_max_iterations)
    {
        frequency_pseudo_min_iterations = frequency_pseudo_max_iterations;
    }
    frequency_block_gs_sweeps =
        SMAX(static_cast<size_t>(1), frequency_block_gs_sweeps);
    frequency_plate_block_sweeps =
        SMAX(static_cast<size_t>(1), frequency_plate_block_sweeps);
    em_delta_a_convergence_tolerance = fabs(em_delta_a_convergence_tolerance);
    em_joule_rel_convergence_tolerance = fabs(em_joule_rel_convergence_tolerance);
    em_relative_residual_convergence_tolerance =
        fabs(em_relative_residual_convergence_tolerance);
    if (!(em_rate_growth_trigger > 1.0))
    {
        em_rate_growth_trigger = std::numeric_limits<Real>::infinity();
    }
    if (!(em_residual_growth_trigger > 1.0))
    {
        em_residual_growth_trigger = std::numeric_limits<Real>::infinity();
    }
    if (!(em_rate_stability_limit_factor > 1.0))
    {
        em_rate_stability_limit_factor = std::numeric_limits<Real>::infinity();
    }
    em_saturation_ratio =
        SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), em_saturation_ratio));
    bool use_voltage_driven_source =
        get_env_bool("TEAM7_USE_VOLTAGE_DRIVEN_SOURCE", false);
    bool use_current_driven_source =
        get_env_bool("TEAM7_USE_CURRENT_DRIVEN_SOURCE",
                     use_frequency_aphi && !use_voltage_driven_source);
    bool use_lumped_coil_source =
        use_current_driven_source || use_voltage_driven_source;
    bool use_frequency_coil_scalar_potential =
        get_env_bool("TEAM7_FREQ_USE_COIL_SCALAR_POTENTIAL",
                     use_frequency_aphi && !use_lumped_coil_source);
    bool use_frequency_coil_magnetic_only =
        get_env_bool("TEAM7_FREQ_USE_MULTITURN_COIL_MAGNETIC_ONLY",
                     use_frequency_aphi && use_lumped_coil_source &&
                         !use_frequency_coil_scalar_potential);
    Real coil_current_rms =
        get_env_real_any("TEAM7_COIL_CURRENT_RMS", use_frequency_aphi ? 1.0 : 0.0);
    Real coil_voltage_rms =
        get_env_real_any("TEAM7_COIL_VOLTAGE_RMS", 0.0);
    Real coil_turns =
        get_env_real("TEAM7_COIL_TURNS", use_frequency_aphi ? 2742.0 : 1.0);
    Real coil_effective_area_input = get_env_real_any("TEAM7_COIL_EFFECTIVE_AREA", -1.0);
    Real coil_current_peak_factor = get_env_real("TEAM7_COIL_CURRENT_PEAK_FACTOR", sqrt(2.0));
    Real coil_circuit_resistance_ohm =
        get_env_real_any("TEAM7_COIL_CIRCUIT_RESISTANCE_OHM", 0.0);
    Real coil_circuit_inductance_h =
        get_env_real_any("TEAM7_COIL_CIRCUIT_INDUCTANCE_H", 0.0);
    Real coil_circuit_series_resistance_ohm =
        get_env_real_any("TEAM7_COIL_CIRCUIT_SERIES_RESISTANCE_OHM", 0.0);
    Real coil_circuit_series_inductance_h =
        get_env_real_any("TEAM7_COIL_CIRCUIT_SERIES_INDUCTANCE_H", 0.0);
    Real geom_length_to_m = get_env_real_any("TEAM7_GEOM_LENGTH_TO_M", 1.0e-3);
    if (geom_length_to_m <= TinyReal)
    {
        geom_length_to_m = 1.0e-3;
    }
    coil_current_rms = SMAX(static_cast<Real>(0.0), coil_current_rms);
    coil_voltage_rms = SMAX(static_cast<Real>(0.0), coil_voltage_rms);
    coil_turns = SMAX(coil_turns, TinyReal);
    coil_current_peak_factor = SMAX(coil_current_peak_factor, static_cast<Real>(1.0));
    coil_circuit_resistance_ohm = SMAX(static_cast<Real>(0.0), coil_circuit_resistance_ohm);
    coil_circuit_inductance_h = SMAX(static_cast<Real>(0.0), coil_circuit_inductance_h);
    coil_circuit_series_resistance_ohm =
        SMAX(static_cast<Real>(0.0), coil_circuit_series_resistance_ohm);
    coil_circuit_series_inductance_h =
        SMAX(static_cast<Real>(0.0), coil_circuit_series_inductance_h);
    Real geom_area_to_m2 = geom_length_to_m * geom_length_to_m;
    Real geom_volume_to_m3 = geom_area_to_m2 * geom_length_to_m;
    Real differential_operator_scaling = 1.0 / geom_length_to_m;
    Real second_order_operator_scaling =
        differential_operator_scaling * differential_operator_scaling;
    Real curl_operator_scaling =
        get_env_real_any("TEAM7_CURL_OPERATOR_SCALING", differential_operator_scaling);
    Real frequency_magnetic_diagonal_scaling =
        get_env_real_any("TEAM7_FREQ_MAG_DIAG_SCALE", curl_operator_scaling * curl_operator_scaling);
    bool use_frequency_balanced_contact_magnetic_diagonal =
        get_env_bool("TEAM7_FREQ_USE_BALANCED_CONTACT_MAG_DIAGONAL", false);
    bool use_frequency_balanced_contact_magnetic_diagonal_plate =
        get_env_bool("TEAM7_FREQ_USE_BALANCED_CONTACT_MAG_DIAGONAL_PLATE",
                     use_frequency_balanced_contact_magnetic_diagonal);
    bool use_frequency_balanced_contact_magnetic_diagonal_coil_air =
        get_env_bool("TEAM7_FREQ_USE_BALANCED_CONTACT_MAG_DIAGONAL_COIL_AIR",
                     use_frequency_balanced_contact_magnetic_diagonal);
    Real frequency_balanced_contact_mag_diagonal_weight =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_BALANCE_WEIGHT",
                         use_frequency_balanced_contact_magnetic_diagonal ? 1.0 : 0.0);
    Real frequency_balanced_contact_mag_diagonal_plate_weight =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_BALANCE_WEIGHT_PLATE",
                         use_frequency_balanced_contact_magnetic_diagonal_plate ? 1.0 : 0.0);
    Real frequency_balanced_contact_mag_diagonal_coil_air_weight =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_BALANCE_WEIGHT_COIL_AIR",
                         use_frequency_balanced_contact_magnetic_diagonal_coil_air ? 1.0 : 0.0);
    Real frequency_balanced_contact_mag_diagonal_weight_requested =
        frequency_balanced_contact_mag_diagonal_weight;
    Real frequency_balanced_contact_mag_diagonal_plate_weight_requested =
        frequency_balanced_contact_mag_diagonal_plate_weight;
    Real frequency_balanced_contact_mag_diagonal_coil_air_weight_requested =
        frequency_balanced_contact_mag_diagonal_coil_air_weight;
    frequency_balanced_contact_mag_diagonal_weight =
        SMIN(static_cast<Real>(1.0),
             SMAX(static_cast<Real>(0.0), frequency_balanced_contact_mag_diagonal_weight));
    frequency_balanced_contact_mag_diagonal_plate_weight =
        SMIN(static_cast<Real>(1.0),
             SMAX(static_cast<Real>(0.0), frequency_balanced_contact_mag_diagonal_plate_weight));
    frequency_balanced_contact_mag_diagonal_coil_air_weight =
        SMIN(static_cast<Real>(1.0),
             SMAX(static_cast<Real>(0.0), frequency_balanced_contact_mag_diagonal_coil_air_weight));
    bool use_frequency_contact_mag_diag_auto_cap =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP", false);
    bool use_frequency_contact_mag_diag_auto_cap_plate =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP_PLATE",
                     use_frequency_contact_mag_diag_auto_cap);
    bool use_frequency_contact_mag_diag_auto_cap_coil_air =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP_COIL_AIR",
                     use_frequency_aphi || use_frequency_contact_mag_diag_auto_cap);
    Real frequency_contact_mag_diag_ratio_cap =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_RATIO_CAP", 0.0);
    Real frequency_contact_mag_diag_ratio_cap_plate =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_RATIO_CAP_PLATE",
                         frequency_contact_mag_diag_ratio_cap);
    Real frequency_contact_mag_diag_ratio_cap_coil_air =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_RATIO_CAP_COIL_AIR",
                         frequency_contact_mag_diag_ratio_cap);
    Real frequency_contact_mag_diag_adaptive_cap_strength =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_ADAPTIVE_CAP_STRENGTH", 0.0);
    Real frequency_contact_mag_diag_adaptive_cap_strength_plate =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_ADAPTIVE_CAP_STRENGTH_PLATE",
                         frequency_contact_mag_diag_adaptive_cap_strength);
    Real frequency_contact_mag_diag_adaptive_cap_strength_coil_air =
        get_env_real_any("TEAM7_FREQ_CONTACT_MAG_DIAG_ADAPTIVE_CAP_STRENGTH_COIL_AIR",
                         frequency_contact_mag_diag_adaptive_cap_strength);
    bool frequency_contact_mag_diag_auto_cap_prefer_conservative =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP_PREFER_CONSERVATIVE", false);
    bool frequency_contact_mag_diag_auto_cap_prefer_conservative_plate =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP_PREFER_CONSERVATIVE_PLATE",
                     frequency_contact_mag_diag_auto_cap_prefer_conservative);
    bool frequency_contact_mag_diag_auto_cap_prefer_conservative_coil_air =
        get_env_bool("TEAM7_FREQ_CONTACT_MAG_DIAG_AUTO_CAP_PREFER_CONSERVATIVE_COIL_AIR",
                     use_frequency_aphi || frequency_contact_mag_diag_auto_cap_prefer_conservative);
    if (use_frequency_contact_mag_diag_auto_cap &&
        !(frequency_contact_mag_diag_ratio_cap > static_cast<Real>(1.0)))
    {
        frequency_contact_mag_diag_ratio_cap = 6.0;
    }
    if (use_frequency_contact_mag_diag_auto_cap_plate &&
        !(frequency_contact_mag_diag_ratio_cap_plate > static_cast<Real>(1.0)))
    {
        frequency_contact_mag_diag_ratio_cap_plate = 6.0;
    }
    if (use_frequency_contact_mag_diag_auto_cap_coil_air &&
        !(frequency_contact_mag_diag_ratio_cap_coil_air > static_cast<Real>(1.0)))
    {
        frequency_contact_mag_diag_ratio_cap_coil_air = 6.0;
    }
    frequency_contact_mag_diag_ratio_cap =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_ratio_cap);
    frequency_contact_mag_diag_ratio_cap_plate =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_ratio_cap_plate);
    frequency_contact_mag_diag_ratio_cap_coil_air =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_ratio_cap_coil_air);
    frequency_contact_mag_diag_adaptive_cap_strength =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_adaptive_cap_strength);
    frequency_contact_mag_diag_adaptive_cap_strength_plate =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_adaptive_cap_strength_plate);
    frequency_contact_mag_diag_adaptive_cap_strength_coil_air =
        SMAX(static_cast<Real>(0.0), frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    if (use_frequency_contact_mag_diag_auto_cap &&
        frequency_contact_mag_diag_auto_cap_prefer_conservative)
    {
        frequency_balanced_contact_mag_diagonal_weight = 0.0;
    }
    if (use_frequency_contact_mag_diag_auto_cap_plate &&
        frequency_contact_mag_diag_auto_cap_prefer_conservative_plate)
    {
        frequency_balanced_contact_mag_diagonal_plate_weight = 0.0;
    }
    if (use_frequency_contact_mag_diag_auto_cap_coil_air &&
        frequency_contact_mag_diag_auto_cap_prefer_conservative_coil_air)
    {
        frequency_balanced_contact_mag_diagonal_coil_air_weight = 0.0;
    }
    if (frequency_balanced_contact_mag_diagonal_weight_requested !=
        frequency_balanced_contact_mag_diagonal_weight)
    {
        std::cout << "[team7-em] overriding freq_contact_mag_diag_balance_weight from "
                  << frequency_balanced_contact_mag_diagonal_weight_requested
                  << " to " << frequency_balanced_contact_mag_diagonal_weight
                  << " due to auto-cap conservative preference." << std::endl;
    }
    if (frequency_balanced_contact_mag_diagonal_plate_weight_requested !=
        frequency_balanced_contact_mag_diagonal_plate_weight)
    {
        std::cout << "[team7-em] overriding freq_contact_mag_diag_balance_weight_plate from "
                  << frequency_balanced_contact_mag_diagonal_plate_weight_requested
                  << " to " << frequency_balanced_contact_mag_diagonal_plate_weight
                  << " due to auto-cap conservative preference." << std::endl;
    }
    if (frequency_balanced_contact_mag_diagonal_coil_air_weight_requested !=
        frequency_balanced_contact_mag_diagonal_coil_air_weight)
    {
        std::cout << "[team7-em] overriding freq_contact_mag_diag_balance_weight_coil_air from "
                  << frequency_balanced_contact_mag_diagonal_coil_air_weight_requested
                  << " to " << frequency_balanced_contact_mag_diagonal_coil_air_weight
                  << " due to auto-cap conservative preference." << std::endl;
    }
    Real phi_divergence_scaling =
        get_env_real_any("TEAM7_PHI_DIV_SCALING", differential_operator_scaling);
    Real phi_laplacian_scaling =
        get_env_real_any("TEAM7_PHI_LAPLACIAN_SCALING", second_order_operator_scaling);
    Real phi_gradient_scaling =
        get_env_real_any("TEAM7_PHI_GRAD_SCALING", differential_operator_scaling);
    Real plate_interface_inner_shell_thickness_input =
        get_env_real_any("TEAM7_PLATE_INTERFACE_INNER_SHELL_THICKNESS", 1.5 * dp_plate);
    Real plate_interface_air_shell_thickness_input =
        get_env_real_any("TEAM7_PLATE_INTERFACE_AIR_SHELL_THICKNESS", 2.0 * dp_air_finest);
    Real coil_interface_inner_shell_thickness_input =
        get_env_real_any("TEAM7_COIL_INTERFACE_INNER_SHELL_THICKNESS", 1.5 * dp_coil);
    Real coil_interface_air_shell_thickness_input =
        get_env_real_any("TEAM7_COIL_INTERFACE_AIR_SHELL_THICKNESS", 2.0 * dp_air_finest);
    Real b_curve_interpolation_radius =
        get_env_real_any("TEAM7_B_CURVE_INTERP_RADIUS", 2.5 * dp_air_finest);
    Real j_curve_interpolation_radius =
        get_env_real_any("TEAM7_J_CURVE_INTERP_RADIUS", 2.0 * dp_plate);
    Real coil_source_dir_x = get_env_real_any("TEAM7_COIL_SOURCE_DIR_X", harmonic_source_amplitude[0]);
    Real coil_source_dir_y = get_env_real_any("TEAM7_COIL_SOURCE_DIR_Y", harmonic_source_amplitude[1]);
    Real coil_source_dir_z = get_env_real_any("TEAM7_COIL_SOURCE_DIR_Z", harmonic_source_amplitude[2]);
    bool use_circular_coil_source =
        get_env_bool("TEAM7_USE_CIRCULAR_COIL_SOURCE", use_frequency_aphi);
    Vecd coil_source_axis(get_env_real_any("TEAM7_COIL_SOURCE_AXIS_X", 0.0),
                          get_env_real_any("TEAM7_COIL_SOURCE_AXIS_Y", 0.0),
                          get_env_real_any("TEAM7_COIL_SOURCE_AXIS_Z", 1.0));
    Real runtime_source_scale = get_env_real_any("TEAM7_RUNTIME_SOURCE_SCALE", 1.0);
    bool auto_normalize_source = get_env_bool("TEAM7_AUTO_NORMALIZE_SOURCE", false);
    size_t auto_normalize_interval = get_env_size_t("TEAM7_AUTO_NORMALIZE_INTERVAL", 1);
    Real auto_normalize_gain =
        get_env_real_any("TEAM7_AUTO_NORMALIZE_GAIN", 0.25);
    Real auto_normalize_min_factor =
        get_env_real_any("TEAM7_AUTO_NORMALIZE_MIN_FACTOR", 0.8);
    Real auto_normalize_max_factor =
        get_env_real_any("TEAM7_AUTO_NORMALIZE_MAX_FACTOR", 1.25);
    Real runtime_source_scale_min =
        get_env_real_any("TEAM7_RUNTIME_SOURCE_SCALE_MIN", 1.0e-4);
    Real runtime_source_scale_max =
        get_env_real_any("TEAM7_RUNTIME_SOURCE_SCALE_MAX", 1.0e4);
    bool auto_normalize_use_j_curve =
        get_env_bool("TEAM7_AUTO_NORMALIZE_USE_J_CURVE", true);
    bool enable_circuit_ni_feedback =
        get_env_bool("TEAM7_ENABLE_CIRCUIT_NI_FEEDBACK", use_voltage_driven_source);
    size_t circuit_ni_feedback_interval =
        get_env_size_t("TEAM7_CIRCUIT_NI_FEEDBACK_INTERVAL", 1);
    Real circuit_ni_feedback_gain =
        get_env_real_any("TEAM7_CIRCUIT_NI_FEEDBACK_GAIN", 0.25);
    Real circuit_ni_feedback_min_factor =
        get_env_real_any("TEAM7_CIRCUIT_NI_FEEDBACK_MIN_FACTOR", 0.8);
    Real circuit_ni_feedback_max_factor =
        get_env_real_any("TEAM7_CIRCUIT_NI_FEEDBACK_MAX_FACTOR", 1.25);
    bool em_convergence_require_rate =
        get_env_bool("TEAM7_EM_CONVERGENCE_REQUIRE_RATE", true);
    bool em_convergence_require_enhanced =
        get_env_bool("TEAM7_EM_CONVERGENCE_REQUIRE_ENHANCED", false);

    auto_normalize_gain = SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), auto_normalize_gain));
    auto_normalize_min_factor = SMAX(auto_normalize_min_factor, TinyReal);
    auto_normalize_max_factor = SMAX(auto_normalize_max_factor, auto_normalize_min_factor);
    runtime_source_scale_min = SMAX(runtime_source_scale_min, TinyReal);
    runtime_source_scale_max = SMAX(runtime_source_scale_max, runtime_source_scale_min);
    runtime_source_scale = SMIN(runtime_source_scale_max, SMAX(runtime_source_scale_min, runtime_source_scale));
    auto_normalize_interval = SMAX(static_cast<size_t>(1), auto_normalize_interval);
    circuit_ni_feedback_gain =
        SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), circuit_ni_feedback_gain));
    circuit_ni_feedback_min_factor = SMAX(circuit_ni_feedback_min_factor, TinyReal);
    circuit_ni_feedback_max_factor =
        SMAX(circuit_ni_feedback_max_factor, circuit_ni_feedback_min_factor);
    circuit_ni_feedback_interval = SMAX(static_cast<size_t>(1), circuit_ni_feedback_interval);
    bool constrain_phi_boundary_plate_enabled = get_env_bool("TEAM7_CONSTRAIN_PHI_BOUNDARY_PLATE", false);
    bool constrain_phi_boundary_coil_enabled = get_env_bool("TEAM7_CONSTRAIN_PHI_BOUNDARY_COIL", false);
    bool constrain_phi_boundary_air_enabled = get_env_bool("TEAM7_CONSTRAIN_PHI_BOUNDARY_AIR", true);
    bool constrain_phi_all_air_enabled = get_env_bool("TEAM7_CONSTRAIN_PHI_ALL_AIR", false);
    bool constrain_a_reference_plate_enabled = get_env_bool("TEAM7_CONSTRAIN_A_REFERENCE_PLATE", true);
    bool constrain_a_reference_coil_enabled = get_env_bool("TEAM7_CONSTRAIN_A_REFERENCE_COIL", true);
    bool constrain_a_reference_air_enabled = get_env_bool("TEAM7_CONSTRAIN_A_REFERENCE_AIR", true);
    bool constrain_a_boundary_air_enabled = get_env_bool("TEAM7_CONSTRAIN_A_BOUNDARY_AIR", true);
    bool use_frequency_air_magnetic_insulation_boundary =
        get_env_bool("TEAM7_FREQ_USE_AIR_MAGNETIC_INSULATION_BOUNDARY", use_frequency_aphi);
    bool use_frequency_scalar_contact_coupling =
        get_env_bool("TEAM7_FREQ_USE_SCALAR_CONTACT_COUPLING", false);
    bool use_frequency_air_scalar_potential =
        get_env_bool("TEAM7_FREQ_USE_AIR_SCALAR_POTENTIAL", false);
    bool use_frequency_air_post_sweep =
        get_env_bool("TEAM7_FREQ_AIR_POST_SWEEP", false);
    bool use_frequency_split_air_sweep =
        get_env_bool("TEAM7_FREQ_SPLIT_AIR_SWEEP", use_frequency_aphi);
    bool use_frequency_operator_coil_air_only =
        get_env_bool("TEAM7_FREQ_OPERATOR_COIL_AIR_ONLY", false);
    bool use_frequency_air_inner_term =
        get_env_bool("TEAM7_FREQ_USE_AIR_INNER_TERM", true);
    bool use_frequency_air_contact_term =
        get_env_bool("TEAM7_FREQ_USE_AIR_CONTACT_TERM", true);
    bool use_frequency_coil_magnetic_block_solver =
        get_env_bool("TEAM7_FREQ_USE_COIL_MAGNETIC_BLOCK_SOLVER",
                     use_frequency_aphi && use_frequency_coil_magnetic_only);
    bool use_frequency_plate_component_hessian_inner =
        get_env_bool("TEAM7_USE_PLATE_COMPONENT_HESSIAN_INNER", false);
    bool use_frequency_coil_component_hessian_inner =
        get_env_bool("TEAM7_USE_COIL_COMPONENT_HESSIAN_INNER", false);
    bool enable_coil_component_hessian_diagnostics =
        get_env_bool("TEAM7_COIL_COMPONENT_HESSIAN_DIAGNOSTICS", false);
    bool enable_coil_component_hessian_core_quality_gate =
        get_env_bool("TEAM7_COIL_COMPONENT_HESSIAN_CORE_QUALITY_GATE", true);
    Real coil_component_hessian_core_gate_rel_tol =
        get_env_real_any("TEAM7_COIL_COMPONENT_HESSIAN_CORE_GATE_REL_TOL", 2.0);
    Real coil_component_hessian_core_gate_abs_ref =
        get_env_real_any("TEAM7_COIL_COMPONENT_HESSIAN_CORE_GATE_ABS_REF", 1.0e4);
    Real coil_component_hessian_outlier_ratio_threshold =
        get_env_real_any("TEAM7_COIL_COMPONENT_HESSIAN_OUTLIER_RATIO", 1.0e3);
    Real coil_component_hessian_outlier_abs_threshold =
        get_env_real_any("TEAM7_COIL_COMPONENT_HESSIAN_OUTLIER_ABS", 1.0e3);
    bool use_frequency_plate_backtracking =
        get_env_bool("TEAM7_FREQ_PLATE_BACKTRACKING", false);
    size_t frequency_plate_backtracking_max_trials =
        get_env_size_t("TEAM7_FREQ_PLATE_BACKTRACK_MAX_TRIALS", 4);
    bool use_frequency_global_backtracking =
        get_env_bool("TEAM7_FREQ_GLOBAL_BACKTRACKING", false);
    size_t frequency_global_backtracking_max_trials =
        get_env_size_t("TEAM7_FREQ_GLOBAL_BACKTRACK_MAX_TRIALS", 4);
    Real frequency_global_backtracking_air_penalty =
        get_env_real_any("TEAM7_FREQ_GLOBAL_BACKTRACK_AIR_PENALTY",
                         use_frequency_aphi ? 0.05 : 0.0);
    Real frequency_global_backtracking_metric_growth_limit =
        get_env_real_any("TEAM7_FREQ_GLOBAL_BACKTRACK_METRIC_GROWTH_LIMIT", 1.0);
    bool use_frequency_hard_guard =
        get_env_bool("TEAM7_FREQ_HARD_GUARD", use_frequency_aphi);
    size_t frequency_hard_guard_warmup_iterations =
        get_env_size_t("TEAM7_FREQ_HARD_GUARD_WARMUP_ITERS",
                       em_adaptive_dt_warmup_iterations);
    size_t frequency_hard_guard_max_trials =
        get_env_size_t("TEAM7_FREQ_HARD_GUARD_MAX_TRIALS", 3);
    size_t frequency_hard_guard_failfast_streak =
        get_env_size_t("TEAM7_FREQ_HARD_GUARD_FAILFAST_STREAK", 6);
    Real frequency_hard_guard_metric_growth_limit =
        get_env_real_any("TEAM7_FREQ_HARD_GUARD_METRIC_GROWTH_LIMIT", 1.05);
    Real frequency_hard_guard_probe_growth_limit =
        get_env_real_any("TEAM7_FREQ_HARD_GUARD_PROBE_GROWTH_LIMIT", 2.0);
    Real frequency_hard_guard_probe_abs_tol =
        get_env_real_any("TEAM7_FREQ_HARD_GUARD_PROBE_ABS_TOL", 1.0e-6);
    Real frequency_hard_guard_dt_reduce_factor =
        get_env_real_any("TEAM7_FREQ_HARD_GUARD_DT_REDUCE_FACTOR", 0.5);
    bool frequency_hard_guard_use_probe_ref_cap =
        get_env_bool("TEAM7_FREQ_HARD_GUARD_USE_PROBE_REF_CAP", true);
    Real frequency_hard_guard_probe_ref_ratio_limit =
        get_env_real_any("TEAM7_FREQ_HARD_GUARD_PROBE_REF_RATIO_LIMIT", 50.0);
    bool frequency_failfast_divergence_guard =
        get_env_bool("TEAM7_FREQ_FAILFAST_DIVERGENCE_GUARD", use_frequency_aphi);
    Real frequency_failfast_residual_ratio_limit =
        get_env_real_any("TEAM7_FREQ_FAILFAST_RESIDUAL_RATIO_LIMIT", 5.0);
    Real frequency_failfast_em_rate_limit =
        get_env_real_any("TEAM7_FREQ_FAILFAST_EM_RATE_LIMIT",
                         10.0 * em_rate_divergence_threshold);
    Real frequency_failfast_air_residual_limit =
        get_env_real_any("TEAM7_FREQ_FAILFAST_AIR_RESIDUAL_LIMIT", 5.0);
    Real frequency_failfast_dt_reduce_factor =
        get_env_real_any("TEAM7_FREQ_FAILFAST_DT_REDUCE_FACTOR",
                         frequency_hard_guard_dt_reduce_factor);
    bool em_rate_use_effective_update =
        get_env_bool("TEAM7_EM_RATE_USE_EFFECTIVE_UPDATE", true);
    bool em_rate_use_pseudo_dt_scaling =
        get_env_bool("TEAM7_EM_RATE_USE_PSEUDO_DT_SCALING", use_frequency_aphi);
    bool em_rate_include_air =
        get_env_bool("TEAM7_EM_RATE_INCLUDE_AIR", use_frequency_aphi);
    Real end_time = get_env_real("TEAM7_END_TIME", end_time_default);
    size_t output_steps = get_env_size_t("TEAM7_OUTPUT_STEPS", output_steps_default);
    Real output_interval = end_time / static_cast<Real>(output_steps);
    bool record_last_step_only =
        get_env_bool("TEAM7_RECORD_LAST_STEP_ONLY", false);
    frequency_plate_backtracking_max_trials =
        SMAX(static_cast<size_t>(1), frequency_plate_backtracking_max_trials);
    frequency_global_backtracking_max_trials =
        SMAX(static_cast<size_t>(1), frequency_global_backtracking_max_trials);
    frequency_hard_guard_warmup_iterations =
        SMAX(static_cast<size_t>(0), frequency_hard_guard_warmup_iterations);
    frequency_hard_guard_max_trials =
        SMAX(static_cast<size_t>(1), frequency_hard_guard_max_trials);
    frequency_hard_guard_failfast_streak =
        SMAX(static_cast<size_t>(1), frequency_hard_guard_failfast_streak);
    frequency_global_backtracking_air_penalty =
        SMAX(static_cast<Real>(0.0), frequency_global_backtracking_air_penalty);
    frequency_hard_guard_metric_growth_limit =
        SMAX(static_cast<Real>(1.0), frequency_hard_guard_metric_growth_limit);
    frequency_hard_guard_probe_growth_limit =
        SMAX(static_cast<Real>(1.0), frequency_hard_guard_probe_growth_limit);
    frequency_hard_guard_probe_abs_tol =
        SMAX(static_cast<Real>(0.0), frequency_hard_guard_probe_abs_tol);
    if (!(frequency_hard_guard_dt_reduce_factor > 0.0 &&
          frequency_hard_guard_dt_reduce_factor < 1.0))
    {
        frequency_hard_guard_dt_reduce_factor = 0.5;
    }
    if (!(frequency_failfast_dt_reduce_factor > 0.0 &&
          frequency_failfast_dt_reduce_factor < 1.0))
    {
        frequency_failfast_dt_reduce_factor = frequency_hard_guard_dt_reduce_factor;
    }
    frequency_hard_guard_probe_ref_ratio_limit =
        SMAX(static_cast<Real>(1.0), frequency_hard_guard_probe_ref_ratio_limit);
    frequency_failfast_residual_ratio_limit =
        SMAX(static_cast<Real>(0.0), frequency_failfast_residual_ratio_limit);
    frequency_failfast_em_rate_limit =
        SMAX(static_cast<Real>(0.0), frequency_failfast_em_rate_limit);
    frequency_failfast_air_residual_limit =
        SMAX(static_cast<Real>(0.0), frequency_failfast_air_residual_limit);
    coil_component_hessian_outlier_ratio_threshold =
        SMAX(static_cast<Real>(1.0), coil_component_hessian_outlier_ratio_threshold);
    coil_component_hessian_outlier_abs_threshold =
        SMAX(static_cast<Real>(0.0), coil_component_hessian_outlier_abs_threshold);
    coil_component_hessian_core_gate_rel_tol =
        SMAX(static_cast<Real>(0.0), coil_component_hessian_core_gate_rel_tol);
    coil_component_hessian_core_gate_abs_ref =
        SMAX(static_cast<Real>(0.0), coil_component_hessian_core_gate_abs_ref);
    frequency_global_backtracking_metric_growth_limit =
        SMAX(static_cast<Real>(1.0), frequency_global_backtracking_metric_growth_limit);
    bool enable_validation_output = get_env_bool("TEAM7_ENABLE_VALIDATION_OUTPUT", true);
    Real validation_rel_tolerance = get_env_real("TEAM7_VALIDATION_REL_TOL", 0.10);
    Real validation_abs_tolerance = get_env_real("TEAM7_VALIDATION_ABS_TOL", 1.0e-8);
    Real nan_reference = std::numeric_limits<Real>::quiet_NaN();
    Real ref_total_joule_power = get_env_real_any("TEAM7_REF_TOTAL_JOULE_POWER", nan_reference);
    Real ref_avg_temperature = get_env_real_any("TEAM7_REF_AVG_TEMPERATURE", nan_reference);
    Real ref_max_temperature = get_env_real_any("TEAM7_REF_MAX_TEMPERATURE", nan_reference);
    Real ref_b_air_above_mag = get_env_real_any("TEAM7_REF_B_AIR_ABOVE_MAG", nan_reference);
    Real ref_b_air_below_mag = get_env_real_any("TEAM7_REF_B_AIR_BELOW_MAG", nan_reference);
    Real ref_b_plate_center_mag = get_env_real_any("TEAM7_REF_B_PLATE_CENTER_MAG", nan_reference);
    Real ref_t_plate_center = get_env_real_any("TEAM7_REF_T_PLATE_CENTER", nan_reference);
    Real ref_t_plate_edge_xplus = get_env_real_any("TEAM7_REF_T_PLATE_EDGE_XPLUS", nan_reference);
    Real ref_delta_t_plate_center = get_env_real_any("TEAM7_REF_DELTA_T_PLATE_CENTER", nan_reference);
    Real ref_delta_t_plate_edge_xplus = get_env_real_any("TEAM7_REF_DELTA_T_PLATE_EDGE_XPLUS", nan_reference);
    bool enable_curve_validation = get_env_bool("TEAM7_ENABLE_CURVE_VALIDATION", true);
    Real curve_validation_rel_tolerance = get_env_real("TEAM7_CURVE_VALIDATION_REL_TOL", 0.10);
    Real curve_validation_abs_tolerance = get_env_real("TEAM7_CURVE_VALIDATION_ABS_TOL", 1.0e-6);
    Real curve_reference_frequency_hz =
        get_env_real("TEAM7_CURVE_REF_FREQUENCY_HZ", harmonic_frequency_hz_runtime);
    Real b_curve_line_y = get_env_real_any("TEAM7_B_CURVE_LINE_Y", 72.0);
    Real b_curve_line_z = get_env_real_any("TEAM7_B_CURVE_LINE_Z", 34.0);
    Real j_curve_line_y = get_env_real_any("TEAM7_J_CURVE_LINE_Y", 72.0);
    Real j_curve_line_z = get_env_real_any("TEAM7_J_CURVE_LINE_Z", 19.0);
    /** Squared-distance metric for J-curve plate fallback: dx^2+dy^2+w*dz^2 (larger w favors target z). */
    Real j_curve_plate_z_distance_weight =
        get_env_real("TEAM7_J_CURVE_PLATE_Z_DISTANCE_WEIGHT", 25.0);
    /** Initial |dz| half-width (mm) for primary J-curve search; <=0 uses dp_plate, then doubles until cap. */
    Real j_curve_plate_dz_band_start_mm =
        get_env_real_any("TEAM7_J_CURVE_PLATE_DZ_BAND_START_MM", -1.0);
    Real j_curve_plate_dz_band_cap_mm =
        get_env_real_any("TEAM7_J_CURVE_PLATE_DZ_BAND_CAP_MM", 40.0);
    Real b_curve_unit_scale = get_env_real("TEAM7_B_CURVE_UNIT_SCALE", 1.0e3);
    Real j_curve_unit_scale = get_env_real("TEAM7_J_CURVE_UNIT_SCALE", 1.0);
    bool enable_j_curve_decomposition_output =
        get_env_bool("TEAM7_ENABLE_J_CURVE_DECOMPOSITION_OUTPUT", false);
    int b_curve_component_index = parse_vector_component_index(
        get_env_string("TEAM7_B_CURVE_COMPONENT", "z"), 2);
    int j_curve_component_index = parse_vector_component_index(
        get_env_string("TEAM7_J_CURVE_COMPONENT", "y"), 1);
    std::string default_curve_complex_mode_token =
        get_env_string("TEAM7_CURVE_COMPLEX_MODE", "signed_magnitude");
    ComplexCurveComponentMode b_curve_complex_mode =
        parse_complex_curve_component_mode(
            get_env_string("TEAM7_B_CURVE_COMPLEX_MODE", default_curve_complex_mode_token),
            ComplexCurveComponentMode::SignedMagnitude);
    ComplexCurveComponentMode j_curve_complex_mode =
        parse_complex_curve_component_mode(
            get_env_string("TEAM7_J_CURVE_COMPLEX_MODE", default_curve_complex_mode_token),
            ComplexCurveComponentMode::SignedMagnitude);
    bool em_only_diagnostics_mode = get_env_bool("TEAM7_EM_ONLY_DIAGNOSTICS", false);
    size_t em_only_max_steps = get_env_size_t("TEAM7_EM_ONLY_MAX_STEPS", 2000);
    size_t em_only_guarded_min_steps =
        get_env_size_t("TEAM7_EM_ONLY_GUARDED_MIN_STEPS", 100);
    bool em_only_force_stability_guards =
        get_env_bool("TEAM7_EM_ONLY_FORCE_STABILITY_GUARDS",
                     use_frequency_aphi &&
                         frequency_contact_mag_diag_auto_cap_prefer_conservative_coil_air);
    bool run_operator_verification =
        get_env_bool("TEAM7_RUN_OPERATOR_VERIFICATION", false);
    bool allow_experimental_plate_component_hessian =
        get_env_bool("TEAM7_ALLOW_EXPERIMENTAL_PLATE_COMPONENT_HESSIAN",
                     run_operator_verification);
    Real operator_verification_margin_dp =
        get_env_real("TEAM7_OPERATOR_VERIFY_MARGIN_DP", 3.0);
    em_only_guarded_min_steps = SMAX(static_cast<size_t>(1), em_only_guarded_min_steps);
    bool em_only_long_run_needs_guard =
        use_frequency_aphi && em_only_diagnostics_mode &&
        em_only_max_steps >= em_only_guarded_min_steps &&
        !use_frequency_hard_guard && !frequency_failfast_divergence_guard;
    if (em_only_long_run_needs_guard)
    {
        if (em_only_force_stability_guards)
        {
            use_frequency_hard_guard = true;
            frequency_failfast_divergence_guard = true;
            std::cout << "[team7-em] enabled hard-guard + fail-fast for long EM-only diagnostics "
                      << "(TEAM7_EM_ONLY_FORCE_STABILITY_GUARDS=1, max_steps="
                      << em_only_max_steps << ")." << std::endl;
        }
        else
        {
            std::cout << "[team7-em-warning] long EM-only diagnostics without stability guards "
                      << "(max_steps=" << em_only_max_steps
                      << ") may diverge. Consider enabling TEAM7_FREQ_HARD_GUARD=1 and "
                      << "TEAM7_FREQ_FAILFAST_DIVERGENCE_GUARD=1." << std::endl;
        }
    }
    if (use_frequency_aphi &&
        use_frequency_plate_component_hessian_inner &&
        !run_operator_verification &&
        !allow_experimental_plate_component_hessian)
    {
        std::cout << "[team7-em] disabling TEAM7_USE_PLATE_COMPONENT_HESSIAN_INNER "
                  << "for runtime solve (experimental branch can suppress plate induction). "
                  << "Set TEAM7_ALLOW_EXPERIMENTAL_PLATE_COMPONENT_HESSIAN=1 to force-enable."
                  << std::endl;
        use_frequency_plate_component_hessian_inner = false;
    }
    bool enable_coil_air_vacuum_mode =
        get_env_bool("TEAM7_COIL_AIR_VACUUM_MODE", false);
    bool enable_coil_air_probe_feedback =
        get_env_bool("TEAM7_ENABLE_COIL_AIR_PROBE_FEEDBACK",
                     enable_coil_air_vacuum_mode && use_frequency_aphi);
    size_t coil_air_probe_feedback_interval =
        get_env_size_t("TEAM7_COIL_AIR_PROBE_FEEDBACK_INTERVAL", 1);
    Real coil_air_probe_feedback_target_ratio =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_FEEDBACK_TARGET_RATIO", 1.0);
    Real coil_air_probe_feedback_gain =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_FEEDBACK_GAIN", 0.5);
    Real coil_air_probe_feedback_min_factor =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_FEEDBACK_MIN_FACTOR", 0.1);
    Real coil_air_probe_feedback_max_factor =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_FEEDBACK_MAX_FACTOR", 1.05);
    bool enable_vacuum_reluctivity_feedback =
        get_env_bool("TEAM7_ENABLE_VACUUM_RELUCTIVITY_FEEDBACK",
                     enable_coil_air_vacuum_mode && use_frequency_aphi);
    size_t vacuum_reluctivity_feedback_interval =
        get_env_size_t("TEAM7_VACUUM_RELUCTIVITY_FEEDBACK_INTERVAL", 1);
    Real vacuum_reluctivity_feedback_target_ratio =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_TARGET_SIM_OVER_REF", 1.0);
    Real vacuum_reluctivity_feedback_gain =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_FEEDBACK_GAIN", 0.5);
    Real vacuum_reluctivity_feedback_min_factor =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_FEEDBACK_MIN_FACTOR", 0.8);
    Real vacuum_reluctivity_feedback_max_factor =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_FEEDBACK_MAX_FACTOR", 1.25);
    Real vacuum_reluctivity_feedback_max_residual_ratio =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_FEEDBACK_MAX_RESIDUAL_RATIO", 5.0);
    Real vacuum_reluctivity_scale_min =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_SCALE_MIN", 1.0e-3);
    Real vacuum_reluctivity_scale_max =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_SCALE_MAX", 1.0e3);
    Real vacuum_reluctivity_scale_initial =
        get_env_real_any("TEAM7_VACUUM_RELUCTIVITY_SCALE_INITIAL", 1.0);
    bool enable_coil_air_transfer_diagnostics =
        get_env_bool("TEAM7_ENABLE_COIL_AIR_TRANSFER_DIAGNOSTICS",
                     use_frequency_aphi);
    bool enable_contact_magnetic_diagonal_diagnostics =
        get_env_bool("TEAM7_ENABLE_CONTACT_MAG_DIAG_DIAGNOSTICS",
                     use_frequency_aphi);
    bool enable_coil_air_interface_b_diag =
        get_env_bool("TEAM7_ENABLE_COIL_AIR_INTERFACE_B_DIAG", false);
    size_t coil_air_transfer_profile_bins =
        get_env_size_t("TEAM7_COIL_AIR_TRANSFER_PROFILE_BINS", 12);
    bool enable_coil_air_probe_path_diagnostics =
        get_env_bool("TEAM7_ENABLE_COIL_AIR_PROBE_PATH_DIAGNOSTICS",
                     use_frequency_aphi);
    size_t coil_air_probe_path_points =
        get_env_size_t("TEAM7_COIL_AIR_PROBE_PATH_POINTS", 21);
    Real coil_air_probe_path_interp_radius =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_PATH_INTERP_RADIUS",
                         b_curve_interpolation_radius);
    Real coil_air_probe_path_start_offset =
        get_env_real_any("TEAM7_COIL_AIR_PROBE_PATH_START_OFFSET", dp_air_finest);
    bool enable_coil_air_biot_savart_validation =
        get_env_bool("TEAM7_ENABLE_COIL_AIR_BIOT_SAVART_VALIDATION",
                     use_frequency_aphi);
    if (use_frequency_operator_coil_air_only && use_frequency_aphi)
    {
        use_frequency_coil_magnetic_only = true;
        use_frequency_coil_scalar_potential = false;
        use_frequency_air_scalar_potential = false;
        use_frequency_scalar_contact_coupling = false;
        use_frequency_plate_component_hessian_inner = false;
        use_frequency_coil_component_hessian_inner = false;
        use_frequency_plate_backtracking = false;
        use_frequency_adaptive_plate_block_sweeps = false;
        frequency_plate_block_sweeps = 1;
        std::cout << "[team7-em] TEAM7_FREQ_OPERATOR_COIL_AIR_ONLY=1 enabled: "
                  << "coil-air contact operator debug mode (plate EM updates disabled)."
                  << std::endl;
    }
    if (!use_frequency_air_inner_term && !use_frequency_air_contact_term)
    {
        std::cout << "[team7-em] TEAM7_FREQ_USE_AIR_INNER_TERM=0 and "
                  << "TEAM7_FREQ_USE_AIR_CONTACT_TERM=0: air magnetic operator is disabled."
                  << std::endl;
    }
    Real sigma_air_default = use_frequency_air_scalar_potential ? sigma_air : 0.0;
    Real sigma_air_runtime = get_env_real_any("TEAM7_SIGMA_AIR", sigma_air_default);
    Real sigma_plate_default = enable_coil_air_vacuum_mode ? 0.0 : sigma_plate;
    Real sigma_plate_runtime =
        get_env_real_any("TEAM7_SIGMA_PLATE", sigma_plate_default);
    Real coil_air_biot_savart_rel_tolerance =
        get_env_real("TEAM7_COIL_AIR_BIOT_SAVART_REL_TOL", curve_validation_rel_tolerance);
    Real coil_air_biot_savart_abs_tolerance =
        get_env_real("TEAM7_COIL_AIR_BIOT_SAVART_ABS_TOL", curve_validation_abs_tolerance);
    if (use_frequency_operator_coil_air_only)
    {
        sigma_plate_runtime = 0.0;
    }
    sigma_air_runtime = SMAX(static_cast<Real>(0.0), sigma_air_runtime);
    sigma_plate_runtime = SMAX(static_cast<Real>(0.0), sigma_plate_runtime);
    coil_air_probe_feedback_interval =
        SMAX(static_cast<size_t>(1), coil_air_probe_feedback_interval);
    coil_air_probe_feedback_target_ratio =
        SMAX(static_cast<Real>(TinyReal), coil_air_probe_feedback_target_ratio);
    coil_air_probe_feedback_gain =
        SMIN(static_cast<Real>(1.0), SMAX(static_cast<Real>(0.0), coil_air_probe_feedback_gain));
    coil_air_probe_feedback_min_factor =
        SMAX(static_cast<Real>(TinyReal), coil_air_probe_feedback_min_factor);
    coil_air_probe_feedback_max_factor =
        SMAX(coil_air_probe_feedback_min_factor, coil_air_probe_feedback_max_factor);
    vacuum_reluctivity_feedback_interval =
        SMAX(static_cast<size_t>(1), vacuum_reluctivity_feedback_interval);
    vacuum_reluctivity_feedback_target_ratio =
        SMAX(static_cast<Real>(TinyReal), vacuum_reluctivity_feedback_target_ratio);
    vacuum_reluctivity_feedback_gain =
        SMIN(static_cast<Real>(1.0),
             SMAX(static_cast<Real>(0.0), vacuum_reluctivity_feedback_gain));
    vacuum_reluctivity_feedback_min_factor =
        SMAX(static_cast<Real>(TinyReal), vacuum_reluctivity_feedback_min_factor);
    vacuum_reluctivity_feedback_max_factor =
        SMAX(vacuum_reluctivity_feedback_min_factor,
             vacuum_reluctivity_feedback_max_factor);
    vacuum_reluctivity_feedback_max_residual_ratio =
        SMAX(static_cast<Real>(TinyReal),
             vacuum_reluctivity_feedback_max_residual_ratio);
    vacuum_reluctivity_scale_min =
        SMAX(static_cast<Real>(TinyReal), vacuum_reluctivity_scale_min);
    vacuum_reluctivity_scale_max =
        SMAX(vacuum_reluctivity_scale_min, vacuum_reluctivity_scale_max);
    vacuum_reluctivity_scale_initial =
        SMIN(vacuum_reluctivity_scale_max,
             SMAX(vacuum_reluctivity_scale_min, vacuum_reluctivity_scale_initial));
    coil_air_transfer_profile_bins =
        SMAX(static_cast<size_t>(4), coil_air_transfer_profile_bins);
    coil_air_probe_path_points =
        SMAX(static_cast<size_t>(3), coil_air_probe_path_points);
    coil_air_probe_path_interp_radius =
        SMAX(coil_air_probe_path_interp_radius, TinyReal);
    coil_air_probe_path_start_offset =
        SMAX(static_cast<Real>(0.0), coil_air_probe_path_start_offset);
    coil_air_biot_savart_rel_tolerance = SMAX(static_cast<Real>(0.0), coil_air_biot_savart_rel_tolerance);
    coil_air_biot_savart_abs_tolerance = SMAX(static_cast<Real>(0.0), coil_air_biot_savart_abs_tolerance);
    std::string team7_table1_path = get_env_string(
        "TEAM7_TABLE1_PATH",
        "./input/comsol_team7_raw/multiturn_coil_asymmetric_conductor_table1.txt");
    std::string team7_table2_path = get_env_string(
        "TEAM7_TABLE2_PATH",
        "./input/comsol_team7_raw/multiturn_coil_asymmetric_conductor_table2.txt");

    std::cout << std::scientific << std::setprecision(3)
              << "[team7-config] air_box_lower=(" << team7_mesh.air_lower[0] << ","
              << team7_mesh.air_lower[1] << "," << team7_mesh.air_lower[2] << ")"
              << ", air_box_upper=(" << team7_mesh.air_upper[0] << ","
              << team7_mesh.air_upper[1] << "," << team7_mesh.air_upper[2] << ")"
              << ", sph_ref_dp=" << team7_mesh.dp_ref << ", dp_coil=" << team7_mesh.dp_coil
              << ", dp_plate=" << team7_mesh.dp_plate << ", dp_air_finest=" << team7_mesh.dp_air_finest
              << ", air_refinement_levels=" << team7_mesh.air_levels
              << ", uniform_resolution=" << static_cast<int>(team7_mesh.uniform_resolution)
              << ", small_air_box_preset=" << static_cast<int>(team7_mesh.used_small_air_box_preset)
              << "\n"
              << "[team7-config] end_time=" << end_time
              << ", output_steps=" << output_steps
              << ", output_interval=" << output_interval
              << ", record_last_step_only=" << static_cast<int>(record_last_step_only)
              << "\n"
              << "[team7-config] dt_em=" << dt_em
              << ", em_substeps=" << em_substeps_per_thermal
              << ", phi_relax_iterations=" << phi_relax_iterations
              << ", a_relax_iterations=" << a_relax_iterations
              << ", coil_a_relax_scale=" << coil_a_relaxation_scaling
              << ", plate_a_relax_scale=" << plate_a_relaxation_scaling
              << ", air_a_relax_scale=" << air_a_relaxation_scaling
              << ", coil_a_rate_limit=" << coil_a_rate_limit
              << ", plate_a_rate_limit=" << plate_a_rate_limit
              << ", air_a_rate_limit=" << air_a_rate_limit
              << ", phi_abs_limit=" << phi_abs_limit
              << ", freq_sigma_relax_scale=" << freq_sigma_relaxation_scaling
              << ", freq_sigma_relax_floor=" << freq_sigma_relaxation_floor
              << ", harmonic_frequency_hz=" << harmonic_frequency_hz_runtime
              << ", harmonic_phase=" << harmonic_phase_runtime
              << ", use_freq_aphi=" << use_frequency_aphi
              << ", freq_source_imag_scale=" << frequency_source_imag_scale
              << ", freq_coupled_implicit_solver=" << use_frequency_coupled_implicit_solver
              << ", use_freq_pseudo_steady=" << use_frequency_pseudo_steady
              << ", freq_em_solve_once=" << frequency_em_solve_once
              << ", freq_pseudo_min_iters=" << frequency_pseudo_min_iterations
              << ", freq_pseudo_max_iters=" << frequency_pseudo_max_iterations
              << ", freq_block_gs_sweeps=" << frequency_block_gs_sweeps
              << ", freq_plate_block_sweeps=" << frequency_plate_block_sweeps
              << ", freq_adaptive_plate_block_sweeps=" << use_frequency_adaptive_plate_block_sweeps
              << ", freq_coil_mag_block_solver=" << use_frequency_coil_magnetic_block_solver
              << ", em_rate_convergence_tol=" << em_rate_convergence_tolerance
              << ", em_rate_divergence_threshold=" << em_rate_divergence_threshold
              << ", use_enhanced_em_convergence=" << use_enhanced_em_convergence
              << ", em_delta_a_convergence_tol=" << em_delta_a_convergence_tolerance
              << ", em_joule_rel_convergence_tol=" << em_joule_rel_convergence_tolerance
              << ", em_stop_on_saturation=" << em_stop_on_saturation
              << ", em_saturation_ratio=" << em_saturation_ratio
              << ", em_saturation_patience=" << em_saturation_patience
              << ", use_em_adaptive_dt=" << use_em_adaptive_dt
              << ", dt_em_min=" << dt_em_min
              << ", dt_em_max=" << dt_em_max
              << ", em_dt_reduce_factor=" << em_dt_reduce_factor
              << ", em_dt_increase_factor=" << em_dt_increase_factor
              << ", em_dt_instability_patience=" << em_dt_instability_patience
              << ", em_dt_max_reductions_per_solve=" << em_dt_max_reductions_per_solve
              << ", em_dt_recovery_factor=" << em_dt_recovery_factor
              << ", em_dt_recovery_min_factor=" << em_dt_recovery_min_factor
              << ", em_dt_recovery_relative_threshold=" << em_dt_recovery_relative_threshold
              << ", use_em_relative_residual=" << use_em_relative_residual
              << ", em_relative_residual_tol=" << em_relative_residual_convergence_tolerance
              << ", em_rate_growth_trigger=" << em_rate_growth_trigger
              << ", em_residual_growth_trigger=" << em_residual_growth_trigger
              << ", em_rate_stability_limit_factor=" << em_rate_stability_limit_factor
              << ", em_adaptive_dt_warmup_iters=" << em_adaptive_dt_warmup_iterations
              << ", use_current_driven_source=" << use_current_driven_source
              << ", use_voltage_driven_source=" << use_voltage_driven_source
              << ", freq_use_coil_scalar_potential=" << use_frequency_coil_scalar_potential
              << ", freq_use_multiturn_coil_magnetic_only=" << use_frequency_coil_magnetic_only
              << ", freq_global_backtracking=" << use_frequency_global_backtracking
              << ", freq_global_backtrack_max_trials="
              << frequency_global_backtracking_max_trials
              << ", freq_global_backtrack_air_penalty="
              << frequency_global_backtracking_air_penalty
              << ", freq_global_backtrack_metric_growth_limit="
              << frequency_global_backtracking_metric_growth_limit
              << ", freq_hard_guard=" << use_frequency_hard_guard
              << ", freq_hard_guard_warmup_iters="
              << frequency_hard_guard_warmup_iterations
              << ", freq_hard_guard_max_trials="
              << frequency_hard_guard_max_trials
              << ", freq_hard_guard_failfast_streak="
              << frequency_hard_guard_failfast_streak
              << ", freq_hard_guard_metric_growth_limit="
              << frequency_hard_guard_metric_growth_limit
              << ", freq_hard_guard_probe_growth_limit="
              << frequency_hard_guard_probe_growth_limit
              << ", freq_hard_guard_probe_abs_tol="
              << frequency_hard_guard_probe_abs_tol
              << ", freq_hard_guard_dt_reduce_factor="
              << frequency_hard_guard_dt_reduce_factor
              << ", freq_hard_guard_use_probe_ref_cap="
              << frequency_hard_guard_use_probe_ref_cap
              << ", freq_hard_guard_probe_ref_ratio_limit="
              << frequency_hard_guard_probe_ref_ratio_limit
              << ", freq_failfast_guard="
              << frequency_failfast_divergence_guard
              << ", freq_failfast_residual_ratio_limit="
              << frequency_failfast_residual_ratio_limit
              << ", freq_failfast_em_rate_limit="
              << frequency_failfast_em_rate_limit
              << ", freq_failfast_air_residual_limit="
              << frequency_failfast_air_residual_limit
              << ", freq_failfast_dt_reduce_factor="
              << frequency_failfast_dt_reduce_factor
              << ", freq_plate_backtracking=" << use_frequency_plate_backtracking
              << ", freq_plate_backtrack_max_trials="
              << frequency_plate_backtracking_max_trials
              << ", plate_if_inner_shell=" << plate_interface_inner_shell_thickness_input
              << ", plate_if_air_shell=" << plate_interface_air_shell_thickness_input
              << ", coil_if_inner_shell=" << coil_interface_inner_shell_thickness_input
              << ", coil_if_air_shell=" << coil_interface_air_shell_thickness_input
              << ", coil_current_rms=" << coil_current_rms
              << ", coil_voltage_rms=" << coil_voltage_rms
              << ", coil_turns=" << coil_turns
              << ", coil_effective_area_input=" << coil_effective_area_input
              << ", coil_current_peak_factor=" << coil_current_peak_factor
              << ", coil_circuit_R_ohm=" << coil_circuit_resistance_ohm
              << ", coil_circuit_L_h=" << coil_circuit_inductance_h
              << ", coil_circuit_series_R_ohm=" << coil_circuit_series_resistance_ohm
              << ", coil_circuit_series_L_h=" << coil_circuit_series_inductance_h
              << ", geom_length_to_m=" << geom_length_to_m
              << ", differential_operator_scaling=" << differential_operator_scaling
              << ", second_order_operator_scaling=" << second_order_operator_scaling
              << ", curl_operator_scaling=" << curl_operator_scaling
              << ", freq_mag_diag_scale=" << frequency_magnetic_diagonal_scaling
              << ", freq_use_balanced_contact_mag_diagonal="
              << use_frequency_balanced_contact_magnetic_diagonal
              << ", freq_use_balanced_contact_mag_diagonal_plate="
              << use_frequency_balanced_contact_magnetic_diagonal_plate
              << ", freq_use_balanced_contact_mag_diagonal_coil_air="
              << use_frequency_balanced_contact_magnetic_diagonal_coil_air
              << ", freq_contact_mag_diag_balance_weight="
              << frequency_balanced_contact_mag_diagonal_weight
              << ", freq_contact_mag_diag_balance_weight_requested="
              << frequency_balanced_contact_mag_diagonal_weight_requested
              << ", freq_contact_mag_diag_balance_weight_plate="
              << frequency_balanced_contact_mag_diagonal_plate_weight
              << ", freq_contact_mag_diag_balance_weight_plate_requested="
              << frequency_balanced_contact_mag_diagonal_plate_weight_requested
              << ", freq_contact_mag_diag_balance_weight_coil_air="
              << frequency_balanced_contact_mag_diagonal_coil_air_weight
              << ", freq_contact_mag_diag_balance_weight_coil_air_requested="
              << frequency_balanced_contact_mag_diagonal_coil_air_weight_requested
              << ", freq_contact_mag_diag_auto_cap="
              << use_frequency_contact_mag_diag_auto_cap
              << ", freq_contact_mag_diag_auto_cap_plate="
              << use_frequency_contact_mag_diag_auto_cap_plate
              << ", freq_contact_mag_diag_auto_cap_coil_air="
              << use_frequency_contact_mag_diag_auto_cap_coil_air
              << ", freq_contact_mag_diag_auto_cap_prefer_conservative="
              << frequency_contact_mag_diag_auto_cap_prefer_conservative
              << ", freq_contact_mag_diag_auto_cap_prefer_conservative_plate="
              << frequency_contact_mag_diag_auto_cap_prefer_conservative_plate
              << ", freq_contact_mag_diag_auto_cap_prefer_conservative_coil_air="
              << frequency_contact_mag_diag_auto_cap_prefer_conservative_coil_air
              << ", freq_contact_mag_diag_ratio_cap="
              << frequency_contact_mag_diag_ratio_cap
              << ", freq_contact_mag_diag_ratio_cap_plate="
              << frequency_contact_mag_diag_ratio_cap_plate
              << ", freq_contact_mag_diag_ratio_cap_coil_air="
              << frequency_contact_mag_diag_ratio_cap_coil_air
              << ", freq_contact_mag_diag_adaptive_cap_strength="
              << frequency_contact_mag_diag_adaptive_cap_strength
              << ", freq_contact_mag_diag_adaptive_cap_strength_plate="
              << frequency_contact_mag_diag_adaptive_cap_strength_plate
              << ", freq_contact_mag_diag_adaptive_cap_strength_coil_air="
              << frequency_contact_mag_diag_adaptive_cap_strength_coil_air
              << ", phi_div_scaling=" << phi_divergence_scaling
              << ", phi_laplacian_scaling=" << phi_laplacian_scaling
              << ", phi_grad_scaling=" << phi_gradient_scaling
              << ", use_circular_coil_source=" << use_circular_coil_source
              << ", coil_source_axis=(" << coil_source_axis[0]
              << "," << coil_source_axis[1]
              << "," << coil_source_axis[2] << ")"
              << ", runtime_source_scale=" << runtime_source_scale
              << ", auto_normalize_source=" << auto_normalize_source
              << ", auto_normalize_interval=" << auto_normalize_interval
              << ", auto_normalize_gain=" << auto_normalize_gain
              << ", auto_normalize_min_factor=" << auto_normalize_min_factor
              << ", auto_normalize_max_factor=" << auto_normalize_max_factor
              << ", runtime_source_scale_min=" << runtime_source_scale_min
              << ", runtime_source_scale_max=" << runtime_source_scale_max
              << ", auto_normalize_use_j_curve=" << auto_normalize_use_j_curve
              << ", enable_circuit_ni_feedback=" << enable_circuit_ni_feedback
              << ", circuit_ni_feedback_interval=" << circuit_ni_feedback_interval
              << ", circuit_ni_feedback_gain=" << circuit_ni_feedback_gain
              << ", circuit_ni_feedback_min_factor=" << circuit_ni_feedback_min_factor
              << ", circuit_ni_feedback_max_factor=" << circuit_ni_feedback_max_factor
              << ", em_convergence_require_rate=" << em_convergence_require_rate
              << ", em_convergence_require_enhanced=" << em_convergence_require_enhanced
              << ", constrain_phi_boundary_plate=" << constrain_phi_boundary_plate_enabled
              << ", constrain_phi_boundary_coil=" << constrain_phi_boundary_coil_enabled
              << ", constrain_phi_boundary_air=" << constrain_phi_boundary_air_enabled
              << ", constrain_phi_all_air=" << constrain_phi_all_air_enabled
              << ", constrain_a_reference_plate=" << constrain_a_reference_plate_enabled
              << ", constrain_a_reference_coil=" << constrain_a_reference_coil_enabled
              << ", constrain_a_reference_air=" << constrain_a_reference_air_enabled
              << ", constrain_a_boundary_air=" << constrain_a_boundary_air_enabled
              << ", freq_air_magnetic_insulation_boundary="
              << use_frequency_air_magnetic_insulation_boundary
              << ", freq_use_scalar_contact_coupling=" << use_frequency_scalar_contact_coupling
              << ", freq_use_air_scalar_potential=" << use_frequency_air_scalar_potential
              << ", freq_air_post_sweep=" << use_frequency_air_post_sweep
              << ", freq_split_air_sweep=" << use_frequency_split_air_sweep
              << ", freq_operator_coil_air_only="
              << use_frequency_operator_coil_air_only
              << ", freq_use_air_inner_term="
              << use_frequency_air_inner_term
              << ", freq_use_air_contact_term="
              << use_frequency_air_contact_term
              << ", freq_air_block_sweeps=" << frequency_air_block_sweeps
              << ", use_plate_component_hessian_inner=" << use_frequency_plate_component_hessian_inner
              << ", allow_experimental_plate_component_hessian="
              << allow_experimental_plate_component_hessian
              << ", use_coil_component_hessian_inner=" << use_frequency_coil_component_hessian_inner
              << ", coil_component_hessian_diagnostics=" << enable_coil_component_hessian_diagnostics
              << ", coil_component_hessian_core_quality_gate=" << enable_coil_component_hessian_core_quality_gate
              << ", coil_component_hessian_core_gate_rel_tol=" << coil_component_hessian_core_gate_rel_tol
              << ", coil_component_hessian_core_gate_abs_ref=" << coil_component_hessian_core_gate_abs_ref
              << ", coil_component_hessian_outlier_ratio=" << coil_component_hessian_outlier_ratio_threshold
              << ", coil_component_hessian_outlier_abs=" << coil_component_hessian_outlier_abs_threshold
              << ", em_rate_use_effective_update=" << em_rate_use_effective_update
              << ", em_rate_use_pseudo_dt_scaling=" << em_rate_use_pseudo_dt_scaling
              << ", em_rate_include_air=" << em_rate_include_air
              << ", enable_validation_output=" << enable_validation_output
              << ", validation_rel_tol=" << validation_rel_tolerance
              << ", validation_abs_tol=" << validation_abs_tolerance
              << ", enable_curve_validation=" << enable_curve_validation
              << ", curve_validation_rel_tol=" << curve_validation_rel_tolerance
              << ", curve_validation_abs_tol=" << curve_validation_abs_tolerance
              << ", curve_ref_freq_hz=" << curve_reference_frequency_hz
              << ", b_curve_line=(y=" << b_curve_line_y << ",z=" << b_curve_line_z << ")"
              << ", j_curve_line=(y=" << j_curve_line_y << ",z=" << j_curve_line_z << ")"
              << ", b_curve_unit_scale=" << b_curve_unit_scale
              << ", j_curve_unit_scale=" << j_curve_unit_scale
              << ", b_curve_interp_radius=" << b_curve_interpolation_radius
              << ", j_curve_interp_radius=" << j_curve_interpolation_radius
              << ", j_curve_plate_z_distance_weight=" << j_curve_plate_z_distance_weight
              << ", j_curve_plate_dz_band_start_mm=" << j_curve_plate_dz_band_start_mm
              << ", j_curve_plate_dz_band_cap_mm=" << j_curve_plate_dz_band_cap_mm
              << ", b_curve_component=" << b_curve_component_index
              << ", j_curve_component=" << j_curve_component_index
              << ", b_curve_complex_mode=" << complex_curve_component_mode_name(b_curve_complex_mode)
              << ", j_curve_complex_mode=" << complex_curve_component_mode_name(j_curve_complex_mode)
              << ", enable_j_curve_decomposition_output=" << enable_j_curve_decomposition_output
              << ", em_only_diagnostics_mode=" << em_only_diagnostics_mode
              << ", em_only_max_steps=" << em_only_max_steps
              << ", em_only_guarded_min_steps=" << em_only_guarded_min_steps
              << ", em_only_force_stability_guards=" << em_only_force_stability_guards
              << ", run_operator_verification=" << run_operator_verification
              << ", operator_verify_margin_dp=" << operator_verification_margin_dp
              << ", coil_air_vacuum_mode=" << enable_coil_air_vacuum_mode
              << ", coil_air_probe_feedback=" << enable_coil_air_probe_feedback
              << ", coil_air_probe_feedback_interval=" << coil_air_probe_feedback_interval
              << ", coil_air_probe_feedback_target_ratio=" << coil_air_probe_feedback_target_ratio
              << ", coil_air_probe_feedback_gain=" << coil_air_probe_feedback_gain
              << ", coil_air_probe_feedback_min_factor=" << coil_air_probe_feedback_min_factor
              << ", coil_air_probe_feedback_max_factor=" << coil_air_probe_feedback_max_factor
              << ", vacuum_reluctivity_feedback=" << enable_vacuum_reluctivity_feedback
              << ", vacuum_reluctivity_feedback_interval=" << vacuum_reluctivity_feedback_interval
              << ", vacuum_reluctivity_target_sim_over_ref="
              << vacuum_reluctivity_feedback_target_ratio
              << ", vacuum_reluctivity_feedback_gain=" << vacuum_reluctivity_feedback_gain
              << ", vacuum_reluctivity_feedback_min_factor="
              << vacuum_reluctivity_feedback_min_factor
              << ", vacuum_reluctivity_feedback_max_factor="
              << vacuum_reluctivity_feedback_max_factor
              << ", vacuum_reluctivity_feedback_max_residual_ratio="
              << vacuum_reluctivity_feedback_max_residual_ratio
              << ", vacuum_reluctivity_scale_min=" << vacuum_reluctivity_scale_min
              << ", vacuum_reluctivity_scale_max=" << vacuum_reluctivity_scale_max
              << ", vacuum_reluctivity_scale_initial=" << vacuum_reluctivity_scale_initial
              << ", enable_coil_air_transfer_diagnostics="
              << enable_coil_air_transfer_diagnostics
              << ", enable_contact_magnetic_diagonal_diagnostics="
              << enable_contact_magnetic_diagonal_diagnostics
              << ", enable_coil_air_interface_b_diag="
              << enable_coil_air_interface_b_diag
              << ", coil_air_transfer_profile_bins="
              << coil_air_transfer_profile_bins
              << ", enable_coil_air_probe_path_diagnostics="
              << enable_coil_air_probe_path_diagnostics
              << ", coil_air_probe_path_points="
              << coil_air_probe_path_points
              << ", coil_air_probe_path_interp_radius="
              << coil_air_probe_path_interp_radius
              << ", coil_air_probe_path_start_offset="
              << coil_air_probe_path_start_offset
              << ", enable_coil_air_biot_savart_validation="
              << enable_coil_air_biot_savart_validation
              << ", coil_air_biot_savart_rel_tol="
              << coil_air_biot_savart_rel_tolerance
              << ", coil_air_biot_savart_abs_tol="
              << coil_air_biot_savart_abs_tolerance
              << ", sigma_plate=" << sigma_plate_runtime
              << ", sigma_air=" << sigma_air_runtime
              << ", table1_path=" << team7_table1_path
              << ", table2_path=" << team7_table2_path
              << std::fixed << std::setprecision(6)
              << std::endl;
    if (auto_normalize_source && enable_circuit_ni_feedback)
    {
        std::cout << "[team7-config] both TEAM7_AUTO_NORMALIZE_SOURCE and "
                  << "TEAM7_ENABLE_CIRCUIT_NI_FEEDBACK are enabled."
                  << " Both controllers will update TEAM7_RUNTIME_SOURCE_SCALE."
                  << std::endl;
    }

    //------------------------------------------------------------------
    //  Particle generation: coil + plate + adaptive air
    //------------------------------------------------------------------
    SolidBody coil_body(sph_system, makeShared<CoilShape>("Coil"));
    coil_body.defineAdaptation<SPHAdaptation>(1.15, dp_0 / dp_coil);
    coil_body.defineMaterial<Solid>();
    coil_body.defineBodyLevelSetShape();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? coil_body.generateParticles<BaseParticles, Reload>(coil_body.getName())
        : coil_body.generateParticles<BaseParticles, Lattice>();

    SolidBody plate_body(sph_system, makeShared<PlateShape>("Plate"));
    plate_body.defineAdaptation<SPHAdaptation>(1.15, dp_0 / dp_plate);
    plate_body.defineClosure<Solid, IsotropicDiffusion>(
        Solid(), ConstructArgs(std::string("Temperature"), thermal_diffusivity));
    plate_body.defineBodyLevelSetShape();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? plate_body.generateParticles<BaseParticles, Reload>(plate_body.getName())
        : plate_body.generateParticles<BaseParticles, Lattice>();

    auto &inner_boundary_shape = sph_system.addShape<InnerBoundaryShape>("InnerBoundary");
    AdaptiveNearInnerSurface air_adaptation(
        dp_air_coarsest, 1.15, 1.0, air_refinement_levels, &inner_boundary_shape);
    auto &air_body = sph_system.addAdaptiveBody<FluidBody, AdaptiveNearInnerSurface>(
        air_adaptation, makeShared<AirShape>("Air", air_box_lower, air_box_upper));
    air_body.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? air_body.generateParticles<BaseParticles, Reload>(air_body.getName())
        : air_body.generateParticles<BaseParticles, Lattice>();

    if (team7_mesh.uniform_resolution && air_refinement_levels == 0)
    {
        BaseParticles &air_particles = air_body.getBaseParticles();
        Real *air_h_ratio = air_particles.getVariableDataByName<Real>("SmoothingLengthRatio");
        const size_t air_total_particles = air_particles.TotalRealParticles();
        const size_t air_particles_bound = air_particles.ParticlesBound();
        size_t corrected_air_h_ratio_count = 0;
        size_t non_finite_air_h_ratio_count = 0;
        Real max_air_h_ratio_before_clamp = 0.0;
        for (size_t i = 0; i < air_particles_bound; ++i)
        {
            const Real ratio_value = air_h_ratio[i];
            if (!std::isfinite(ratio_value))
            {
                air_h_ratio[i] = 1.0;
                non_finite_air_h_ratio_count++;
                corrected_air_h_ratio_count++;
                continue;
            }
            max_air_h_ratio_before_clamp = std::max(max_air_h_ratio_before_clamp, ratio_value);
            if (ratio_value > 1.0)
            {
                air_h_ratio[i] = 1.0;
                corrected_air_h_ratio_count++;
            }
        }
        if (corrected_air_h_ratio_count > 0)
        {
            std::cout << "[team7-air-h-ratio] corrected=" << corrected_air_h_ratio_count
                      << ", non_finite=" << non_finite_air_h_ratio_count
                      << ", real_particles=" << air_total_particles
                      << ", particles_bound=" << air_particles_bound
                      << ", max_before_clamp=" << max_air_h_ratio_before_clamp
                      << ", mode=uniform_air_levels0" << std::endl;
        }
    }

    const size_t total_coil_particles_generated =
        coil_body.getBaseParticles().TotalRealParticles();
    const size_t total_plate_particles_generated =
        plate_body.getBaseParticles().TotalRealParticles();
    const size_t total_air_particles_generated =
        air_body.getBaseParticles().TotalRealParticles();
    const size_t total_particles_generated =
        total_coil_particles_generated + total_plate_particles_generated + total_air_particles_generated;
    std::cout << "[team7-particles] mode="
              << (sph_system.RunParticleRelaxation() ? "relax" : "solve")
              << ", reload=" << (sph_system.ReloadParticles() ? 1 : 0)
              << ", coil=" << total_coil_particles_generated
              << ", plate=" << total_plate_particles_generated
              << ", air=" << total_air_particles_generated
              << ", total=" << total_particles_generated
              << std::endl;

    if (sph_system.RunParticleRelaxation())
    {
        int particle_relaxation_steps_runtime =
            static_cast<int>(get_env_size_t("TEAM7_PARTICLE_RELAX_STEPS", static_cast<size_t>(particle_relaxation_steps)));
        int particle_relaxation_output_interval_runtime =
            static_cast<int>(get_env_size_t("TEAM7_PARTICLE_RELAX_OUTPUT_INTERVAL",
                                            static_cast<size_t>(particle_relaxation_output_interval)));
        particle_relaxation_steps_runtime = std::max(1, particle_relaxation_steps_runtime);
        particle_relaxation_output_interval_runtime = std::max(1, particle_relaxation_output_interval_runtime);
        std::cout << "[team7-relax] steps=" << particle_relaxation_steps_runtime
                  << ", output_interval=" << particle_relaxation_output_interval_runtime
                  << std::endl;

        using namespace relax_dynamics;
        InnerRelation coil_inner_relax(coil_body);
        InnerRelation plate_inner_relax(plate_body);

        SimpleDynamics<RandomizeParticlePosition> randomize_coil_particles(coil_body);
        SimpleDynamics<RandomizeParticlePosition> randomize_plate_particles(plate_body);
        RelaxationStepLevelSetCorrectionInner relax_coil_step(coil_inner_relax);
        RelaxationStepLevelSetCorrectionInner relax_plate_step(plate_inner_relax);

        BodyStatesRecordingToVtp write_relax_states(sph_system);
        write_relax_states.addToWrite<Real>(air_body, "SmoothingLengthRatio");
        ReloadParticleIO write_particle_reload_files({&coil_body, &plate_body, &air_body});
        write_particle_reload_files.addToReload<Real>(air_body, "SmoothingLengthRatio");

        randomize_coil_particles.exec(particle_randomization_factor);
        randomize_plate_particles.exec(particle_randomization_factor);
        relax_coil_step.SurfaceBounding().exec();
        relax_plate_step.SurfaceBounding().exec();
        write_relax_states.writeToFile(0);

        int ite_p = 0;
        while (ite_p < particle_relaxation_steps_runtime)
        {
            relax_coil_step.exec();
            relax_plate_step.exec();
            ite_p++;
            if (ite_p % particle_relaxation_output_interval_runtime == 0)
            {
                std::cout << std::fixed << std::setprecision(6)
                          << "Particle relaxation N = " << ite_p << std::endl;
                write_relax_states.writeToFile(ite_p);
            }
        }
        write_particle_reload_files.writeToFile(0);
        std::cout << "Particle relaxation finished. Reload files are ready." << std::endl;
        return 0;
    }

    const electromagnetics::ComponentVariableNames<> operator_verify_a_real_component_names = {
        "OperatorVerifyAxReal", "OperatorVerifyAyReal", "OperatorVerifyAzReal"};
    const electromagnetics::ComponentVariableNames<> plate_a_real_component_names = {
        "PlateAxReal", "PlateAyReal", "PlateAzReal"};
    const electromagnetics::ComponentVariableNames<> plate_a_imag_component_names = {
        "PlateAxImag", "PlateAyImag", "PlateAzImag"};
    plate_body.getBaseParticles().registerStateVariableData<Real>(operator_verify_a_real_component_names[0], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(operator_verify_a_real_component_names[1], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(operator_verify_a_real_component_names[2], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[0], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[1], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[2], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[0], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[1], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[2], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[0], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[1], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[2], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[0], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[1], Real(0));
    coil_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[2], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[0], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[1], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_real_component_names[2], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[0], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[1], Real(0));
    air_body.getBaseParticles().registerStateVariableData<Real>(plate_a_imag_component_names[2], Real(0));
    plate_body.getBaseParticles().registerStateVariableData<Vecd>(
        "OperatorVerifyCurlNuBFromComponentHessiansReal", ZeroData<Vecd>::value);

    //------------------------------------------------------------------
    //  Electromagnetic solve setup (coil-plate with contact coupling)
    //------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    InnerRelation coil_inner(coil_body);
    InnerRelation plate_inner(plate_body);
    AdaptiveInnerRelation air_inner(air_body);
    AdaptiveContactRelation coil_contact(coil_body, {&plate_body, &air_body});
    AdaptiveContactRelation plate_contact(plate_body, {&coil_body, &air_body});
    AdaptiveContactRelation air_contact(air_body, {&coil_body, &plate_body});
    AdaptiveContactRelation coil_contact_air_only(coil_body, {&air_body});
    AdaptiveContactRelation air_contact_coil_only(air_body, {&coil_body});
    Inner<> plate_inner_ck(plate_body);
    using PlateCoilContactCK = Contact<Relation<SolidBody, SolidBody>>;
    using AirBodyType = std::remove_reference_t<decltype(air_body)>;
    using PlateAirContactCK = Contact<Relation<SolidBody, AirBodyType>>;
    PlateCoilContactCK plate_coil_contact_ck(plate_body, {&coil_body});
    PlateAirContactCK plate_air_contact_ck(plate_body, {&air_body});
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_plate_cell_linked_list_ck(plate_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_plate_inner_ck(plate_inner_ck);
    UpdateRelation<MainExecutionPolicy,
                   Contact<Relation<SolidBody, SolidBody>>,
                   Contact<Relation<SolidBody, AirBodyType>>>
        update_plate_contact_ck(plate_coil_contact_ck, plate_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        plate_verify_linear_correction_ck(DynamicsArgs(plate_inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>>>
        plate_verify_displacement_matrix_gradient_ck(plate_inner_ck);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>>>
        plate_verify_hessian_correction_ck(DynamicsArgs(plate_inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearCorrectionMatrix<Inner<WithUpdate>,
                                                 Contact<Relation<SolidBody, SolidBody>>,
                                                 Contact<Relation<SolidBody, AirBodyType>>>>
        plate_component_linear_correction_ck(
            DynamicsArgs(plate_inner_ck, 0.0), plate_coil_contact_ck, plate_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy,
                          DisplacementMatrixGradient<Inner<>,
                                                     Contact<Relation<SolidBody, SolidBody>>,
                                                     Contact<Relation<SolidBody, AirBodyType>>>>
        plate_component_displacement_matrix_gradient_ck(
            plate_inner_ck, plate_coil_contact_ck, plate_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy,
                          HessianCorrectionMatrix<Inner<WithUpdate>,
                                                  Contact<Relation<SolidBody, SolidBody>>,
                                                  Contact<Relation<SolidBody, AirBodyType>>>>
        plate_component_hessian_correction_ck(
            DynamicsArgs(plate_inner_ck, 0.0), plate_coil_contact_ck, plate_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_verify_scalar_gradient_ax_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAxReal")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_verify_scalar_hessian_ax_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAxReal")));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_verify_scalar_gradient_ay_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAyReal")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_verify_scalar_hessian_ay_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAyReal")));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_verify_scalar_gradient_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAzReal")));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_verify_scalar_hessian_ck(DynamicsArgs(plate_inner_ck, std::string("OperatorVerifyAzReal")));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_real_gradient_x_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_real_hessian_x_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_gradient_x_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[0]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[0]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_hessian_x_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[0]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[0]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_real_gradient_y_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_real_hessian_y_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_gradient_y_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[1]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[1]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_hessian_y_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[1]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[1]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_real_gradient_z_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_real_hessian_z_ck(DynamicsArgs(plate_inner_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_gradient_z_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[2]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[2]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_real_hessian_z_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_real_component_names[2]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_real_component_names[2]),
            DynamicsArgs(plate_air_contact_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_imag_gradient_x_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_imag_hessian_x_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_gradient_x_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_hessian_x_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_imag_gradient_y_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_imag_hessian_y_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_gradient_y_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_hessian_y_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        plate_a_imag_gradient_z_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        plate_a_imag_hessian_z_ck(DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_gradient_z_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        plate_a_imag_hessian_z_contact_ck(
            DynamicsArgs(plate_inner_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(plate_coil_contact_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(plate_air_contact_ck, plate_a_imag_component_names[2]));

    Inner<> coil_inner_ck(coil_body);
    using CoilPlateContactCK = Contact<Relation<SolidBody, SolidBody>>;
    using CoilAirContactCK = Contact<Relation<SolidBody, AirBodyType>>;
    CoilPlateContactCK coil_plate_contact_ck(coil_body, {&plate_body});
    CoilAirContactCK coil_air_contact_ck(coil_body, {&air_body});
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_coil_cell_linked_list_ck(coil_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_coil_inner_ck(coil_inner_ck);
    UpdateRelation<MainExecutionPolicy,
                   Contact<Relation<SolidBody, SolidBody>>,
                   Contact<Relation<SolidBody, AirBodyType>>>
        update_coil_contact_ck(coil_plate_contact_ck, coil_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        coil_verify_linear_correction_ck(DynamicsArgs(coil_inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy, DisplacementMatrixGradient<Inner<>>>
        coil_verify_displacement_matrix_gradient_ck(coil_inner_ck);
    InteractionDynamicsCK<MainExecutionPolicy, HessianCorrectionMatrix<Inner<WithUpdate>>>
        coil_verify_hessian_correction_ck(DynamicsArgs(coil_inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearCorrectionMatrix<Inner<WithUpdate>,
                                                 Contact<Relation<SolidBody, SolidBody>>,
                                                 Contact<Relation<SolidBody, AirBodyType>>>>
        coil_component_linear_correction_ck(
            DynamicsArgs(coil_inner_ck, 0.0), coil_plate_contact_ck, coil_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy,
                          DisplacementMatrixGradient<Inner<>,
                                                     Contact<Relation<SolidBody, SolidBody>>,
                                                     Contact<Relation<SolidBody, AirBodyType>>>>
        coil_component_displacement_matrix_gradient_ck(
            coil_inner_ck, coil_plate_contact_ck, coil_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy,
                          HessianCorrectionMatrix<Inner<WithUpdate>,
                                                  Contact<Relation<SolidBody, SolidBody>>,
                                                  Contact<Relation<SolidBody, AirBodyType>>>>
        coil_component_hessian_correction_ck(
            DynamicsArgs(coil_inner_ck, 0.0), coil_plate_contact_ck, coil_air_contact_ck);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_real_gradient_x_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_real_hessian_x_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_gradient_x_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[0]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[0]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_hessian_x_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[0]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[0]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_real_gradient_y_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_real_hessian_y_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_gradient_y_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[1]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[1]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_hessian_y_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[1]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[1]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_real_gradient_z_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_real_hessian_z_ck(DynamicsArgs(coil_inner_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_gradient_z_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[2]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[2]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_real_hessian_z_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_real_component_names[2]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_real_component_names[2]),
            DynamicsArgs(coil_air_contact_ck, plate_a_real_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_imag_gradient_x_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_imag_hessian_x_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_gradient_x_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_hessian_x_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[0]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[0]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_imag_gradient_y_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_imag_hessian_y_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_gradient_y_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_hessian_y_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[1]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[1]));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>>
        coil_a_imag_gradient_z_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy, Hessian<Inner<Real>>>
        coil_a_imag_hessian_z_ck(DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          LinearGradient<Inner<Real>,
                                         Contact<Real, Relation<SolidBody, SolidBody>>,
                                         Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_gradient_z_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[2]));
    InteractionDynamicsCK<MainExecutionPolicy,
                          Hessian<Inner<Real>,
                                  Contact<Real, Relation<SolidBody, SolidBody>>,
                                  Contact<Real, Relation<SolidBody, AirBodyType>>>>
        coil_a_imag_hessian_z_contact_ck(
            DynamicsArgs(coil_inner_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(coil_plate_contact_ck, plate_a_imag_component_names[2]),
            DynamicsArgs(coil_air_contact_ck, plate_a_imag_component_names[2]));

    BaseParticles &plate_particles_for_bc = plate_body.getBaseParticles();
    BaseParticles &coil_particles_for_bc = coil_body.getBaseParticles();
    BaseParticles &air_particles_for_bc = air_body.getBaseParticles();
    Vecd *plate_positions = plate_particles_for_bc.getVariableDataByName<Vecd>("Position");
    Vecd *coil_positions = coil_particles_for_bc.getVariableDataByName<Vecd>("Position");
    Vecd *air_positions = air_particles_for_bc.getVariableDataByName<Vecd>("Position");

    BoundingBoxd plate_bounds = plate_body.getSPHBodyBounds();
    Vec3d plate_center = 0.5 * (plate_bounds.lower_ + plate_bounds.upper_);
    Vec3d plate_size = plate_bounds.upper_ - plate_bounds.lower_;
    Vec3d plate_halfsize = 0.5 * plate_size;
    Vec3d plate_phi_ref_point = plate_positions[0];
    BodyRegionByParticle plate_phi_reference_region(
        plate_body,
        makeShared<GeometricShapeBall>(plate_phi_ref_point, 1.5 * dp_plate, "PlatePhiReferencePoint"));
    BodyRegionByParticle plate_a_reference_region(
        plate_body,
        makeShared<GeometricShapeBall>(plate_phi_ref_point, 1.5 * dp_plate, "PlateAReferencePoint"));
    Vec3d plate_phi_boundary_halfsize(0.75 * dp_plate, 0.5 * plate_size[1], 0.5 * plate_size[2]);
    Vec3d plate_phi_boundary_center = plate_center;
    plate_phi_boundary_center[0] = plate_bounds.lower_[0] + plate_phi_boundary_halfsize[0];
    BodyRegionByParticle plate_phi_boundary_region(
        plate_body,
        makeShared<GeometricShapeBox>(Transform(plate_phi_boundary_center), plate_phi_boundary_halfsize, "PlatePhiBoundaryStrip"));

    BoundingBoxd coil_bounds = coil_body.getSPHBodyBounds();
    Vec3d coil_center = 0.5 * (coil_bounds.lower_ + coil_bounds.upper_);
    Vec3d coil_size = coil_bounds.upper_ - coil_bounds.lower_;
    Vec3d coil_halfsize = 0.5 * coil_size;
    Vec3d coil_phi_ref_point = coil_positions[0];
    BodyRegionByParticle coil_phi_reference_region(
        coil_body,
        makeShared<GeometricShapeBall>(coil_phi_ref_point, 1.5 * dp_coil, "CoilPhiReferencePoint"));
    BodyRegionByParticle coil_a_reference_region(
        coil_body,
        makeShared<GeometricShapeBall>(coil_phi_ref_point, 1.5 * dp_coil, "CoilAReferencePoint"));
    Vec3d coil_phi_boundary_halfsize(0.75 * dp_coil, 0.5 * coil_size[1], 0.5 * coil_size[2]);
    Vec3d coil_phi_boundary_center = coil_center;
    coil_phi_boundary_center[0] = coil_bounds.lower_[0] + coil_phi_boundary_halfsize[0];
    BodyRegionByParticle coil_phi_boundary_region(
        coil_body,
        makeShared<GeometricShapeBox>(Transform(coil_phi_boundary_center), coil_phi_boundary_halfsize, "CoilPhiBoundaryStrip"));

    BoundingBoxd air_bounds = air_body.getSPHBodyBounds();
    Vec3d air_center = 0.5 * (air_bounds.lower_ + air_bounds.upper_);
    Vec3d air_size = air_bounds.upper_ - air_bounds.lower_;
    auto clamp_inner_shell_thickness = [](Real thickness, const Vec3d &halfsize) -> Real
    {
        Real max_thickness = 0.95 * halfsize.minCoeff();
        Real abs_thickness = static_cast<Real>(std::fabs(thickness));
        return SMAX(TinyReal, SMIN(abs_thickness, max_thickness));
    };
    Real plate_interface_inner_shell_thickness =
        clamp_inner_shell_thickness(plate_interface_inner_shell_thickness_input, plate_halfsize);
    Real coil_interface_inner_shell_thickness =
        clamp_inner_shell_thickness(coil_interface_inner_shell_thickness_input, coil_halfsize);
    Real plate_interface_air_shell_thickness =
        SMAX(static_cast<Real>(std::fabs(plate_interface_air_shell_thickness_input)), TinyReal);
    Real coil_interface_air_shell_thickness =
        SMAX(static_cast<Real>(std::fabs(coil_interface_air_shell_thickness_input)), TinyReal);
    Real air_a_boundary_shell_thickness = 1.5 * dp_air_coarsest;
    Vec3d air_phi_ref_point = air_positions[0];
    BodyRegionByParticle air_phi_reference_region(
        air_body,
        makeShared<GeometricShapeBall>(air_phi_ref_point, 1.5 * dp_air_coarsest, "AirPhiReferencePoint"));
    BodyRegionByParticle air_a_reference_region(
        air_body,
        makeShared<GeometricShapeBall>(air_phi_ref_point, 1.5 * dp_air_coarsest, "AirAReferencePoint"));
    BodyRegionByParticle air_a_boundary_region(
        air_body,
        makeShared<AirOuterBoundaryShell>("AirABoundaryShell", air_box_lower, air_box_upper,
                                          air_a_boundary_shell_thickness));
    Vec3d air_phi_boundary_halfsize(0.75 * dp_air_coarsest, 0.5 * air_size[1], 0.5 * air_size[2]);
    Vec3d air_phi_boundary_center = air_center;
    air_phi_boundary_center[0] = air_bounds.lower_[0] + air_phi_boundary_halfsize[0];
    BodyRegionByParticle air_phi_boundary_region(
        air_body,
        makeShared<GeometricShapeBox>(Transform(air_phi_boundary_center), air_phi_boundary_halfsize, "AirPhiBoundaryStrip"));
    BodyRegionByParticle air_phi_all_region(
        air_body,
        makeShared<GeometricShapeBox>(Transform(air_center), 0.5 * air_size, "AirPhiAll"));

    extra_electromagnetics::MultiturnCoilDriveConfig coil_drive_config;
    coil_drive_config.use_current_driven_source = use_current_driven_source;
    coil_drive_config.use_voltage_driven_source = use_voltage_driven_source;
    coil_drive_config.use_circular_source = use_circular_coil_source;
    coil_drive_config.current_rms = coil_current_rms;
    coil_drive_config.voltage_rms = coil_voltage_rms;
    coil_drive_config.turns = coil_turns;
    coil_drive_config.current_peak_factor = coil_current_peak_factor;
    coil_drive_config.imag_scale = frequency_source_imag_scale;
    coil_drive_config.circuit_resistance_ohm = coil_circuit_resistance_ohm;
    coil_drive_config.circuit_inductance_h = coil_circuit_inductance_h;
    coil_drive_config.external_series_resistance_ohm = coil_circuit_series_resistance_ohm;
    coil_drive_config.external_series_inductance_h = coil_circuit_series_inductance_h;
    coil_drive_config.effective_area_input_geom = coil_effective_area_input;
    coil_drive_config.geom_length_to_m = geom_length_to_m;
    coil_drive_config.source_direction = Vecd(coil_source_dir_x, coil_source_dir_y, coil_source_dir_z);
    coil_drive_config.source_axis = coil_source_axis;
    coil_drive_config.fallback_harmonic_amplitude = harmonic_source_amplitude;
    extra_electromagnetics::MultiturnCoilDriveState coil_drive_state =
        extra_electromagnetics::BuildMultiturnCoilDriveState(coil_particles_for_bc,
                                                             coil_bounds,
                                                             coil_center,
                                                             harmonic_angular_frequency_runtime,
                                                             coil_drive_config);

    Vecd harmonic_source_amplitude_used = coil_drive_state.harmonic_source_amplitude_used;
    Vecd frequency_source_real = coil_drive_state.frequency_source_real;
    Vecd frequency_source_imag = coil_drive_state.frequency_source_imag;
    Real coil_effective_area_used = coil_drive_state.effective_area_used_geom;
    Real coil_effective_area_used_si = coil_drive_state.effective_area_used_si;
    Real coil_volume_total = coil_drive_state.volume_total_geom;
    Real coil_projection_length = coil_drive_state.projection_length_geom;
    Real coil_mean_radius = coil_drive_state.mean_radius_geom;
    Real coil_volume_total_si = coil_drive_state.volume_total_si;
    Real coil_projection_length_si = coil_drive_state.projection_length_si;
    Vecd coil_source_direction_used = coil_drive_state.source_direction_used;
    Vecd coil_source_axis_used = coil_drive_state.source_axis_used;

    std::cout << std::scientific << std::setprecision(6)
              << "[team7-source] Js_real=(" << frequency_source_real[0] << ", "
              << frequency_source_real[1] << ", "
              << frequency_source_real[2] << ")"
              << ", Js_imag=(" << frequency_source_imag[0] << ", "
              << frequency_source_imag[1] << ", "
              << frequency_source_imag[2] << ")"
              << ", coil_volume_total_geom=" << coil_volume_total
              << ", coil_projection_length_geom=" << coil_projection_length
              << ", coil_mean_radius_geom=" << coil_mean_radius
              << ", coil_effective_area_used_geom=" << coil_effective_area_used
              << ", coil_volume_total_m3=" << coil_volume_total_si
              << ", coil_projection_length_m=" << coil_projection_length_si
              << ", coil_effective_area_used_m2=" << coil_effective_area_used_si
              << std::endl;
    if (coil_drive_state.used_current_driven_source ||
        coil_drive_state.used_voltage_driven_source)
    {
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-source-circuit] current_driven="
                  << coil_drive_state.used_current_driven_source
                  << ", voltage_driven=" << coil_drive_state.used_voltage_driven_source
                  << ", I_rms=" << coil_drive_state.circuit_current_rms
                  << ", I_peak=" << coil_drive_state.circuit_current_peak
                  << ", I_phase_rad=" << coil_drive_state.circuit_current_phase_rad
                  << ", R_total_ohm=" << coil_drive_state.circuit_total_resistance_ohm
                  << ", X_total_ohm="
                  << coil_drive_state.circuit_total_inductive_reactance_ohm
                  << ", |Z|_ohm=" << coil_drive_state.circuit_impedance_magnitude_ohm
                  << std::endl;
    }

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real harmonic_source_magnitude_used = coil_drive_state.harmonic_source_magnitude_used;
    Real frequency_source_real_magnitude = coil_drive_state.frequency_source_real_magnitude;
    Real frequency_source_imag_magnitude = coil_drive_state.frequency_source_imag_magnitude;
    Real target_coil_ampere_turns_real =
        frequency_source_real_magnitude * coil_effective_area_used_si;
    Real target_coil_ampere_turns_imag =
        frequency_source_imag_magnitude * coil_effective_area_used_si;
    Real target_coil_ampere_turns_magnitude =
        sqrt(target_coil_ampere_turns_real * target_coil_ampere_turns_real +
             target_coil_ampere_turns_imag * target_coil_ampere_turns_imag);
    bool apply_runtime_source_scale =
        enable_coil_air_probe_feedback ||
        auto_normalize_source || enable_circuit_ni_feedback ||
        fabs(runtime_source_scale - static_cast<Real>(1.0)) > TinyReal;
    if (frequency_source_real.norm() > TinyReal &&
        frequency_source_real.dot(frequency_source_imag) < 0.0)
    {
        frequency_source_imag_magnitude = -frequency_source_imag_magnitude;
    }

    SimpleDynamics<TemperatureInitialization> temperature_initialization(plate_body, initial_temperature);

    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_aphi_plate(plate_body, sigma_plate_runtime, rho_cp_plate, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_aphi_coil(coil_body, sigma_coil, rho_cp_coil, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiElectromagneticVariables>
        initialize_aphi_air(air_body, sigma_air_runtime, rho_cp_air, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_aphi_frequency_plate(plate_body, sigma_plate_runtime, rho_cp_plate, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_aphi_frequency_coil(coil_body, sigma_coil, rho_cp_coil, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::InitializeAphiFrequencyElectromagneticVariables>
        initialize_aphi_frequency_air(air_body, sigma_air_runtime, rho_cp_air, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_operator_verify_a_real_components(
            plate_body, "VectorPotentialReal", operator_verify_a_real_component_names);
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_operator_verify_curl_nu_b_real_from_components(
            plate_body, operator_verify_a_real_component_names,
            second_order_operator_scaling,
            "OperatorVerifyCurlNuBFromComponentHessiansReal");
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_plate_a_real_components(
            plate_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_coil_a_real_components_for_plate(
            coil_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_air_a_real_components_for_plate(
            air_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_plate_curl_nu_b_real_from_components(
            plate_body, plate_a_real_component_names,
            second_order_operator_scaling,
            "CurlNuBReal");
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_coil_a_real_components(
            coil_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_plate_a_real_components_for_coil(
            plate_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_air_a_real_components_for_coil(
            air_body, "VectorPotentialReal", plate_a_real_component_names);
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_coil_curl_nu_b_real_from_components(
            coil_body, plate_a_real_component_names,
            second_order_operator_scaling,
            "CurlNuBReal");
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_plate_a_imag_components(
            plate_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_coil_a_imag_components_for_plate(
            coil_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_air_a_imag_components_for_plate(
            air_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_plate_curl_nu_b_imag_from_components(
            plate_body, plate_a_imag_component_names,
            second_order_operator_scaling,
            "CurlNuBImag");
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_coil_a_imag_components(
            coil_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_plate_a_imag_components_for_coil(
            plate_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::CopyVectorFieldComponentsToScalarVariables>
        copy_air_a_imag_components_for_coil(
            air_body, "VectorPotentialImag", plate_a_imag_component_names);
    SimpleDynamics<electromagnetics::ReconstructCurlNuBFromScalarComponentHessians>
        reconstruct_coil_curl_nu_b_imag_from_components(
            coil_body, plate_a_imag_component_names,
            second_order_operator_scaling,
            "CurlNuBImag");

    SimpleDynamics<electromagnetics::SetConstantElectromagneticMaterialProperties>
        set_em_material_plate(plate_body, sigma_plate_runtime, rho_cp_plate, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::SetConstantElectromagneticMaterialProperties>
        set_em_material_coil(coil_body, sigma_coil, rho_cp_coil, magnetic_reluctivity);
    SimpleDynamics<electromagnetics::SetConstantElectromagneticMaterialProperties>
        set_em_material_air(air_body, sigma_air_runtime, rho_cp_air, magnetic_reluctivity);

    SimpleDynamics<HarmonicSourceCurrentDensity>
        set_harmonic_source_on_coil(coil_body, harmonic_source_amplitude_used, harmonic_frequency_hz_runtime, physical_time, harmonic_phase_runtime);
    SimpleDynamics<extra_electromagnetics::CircularHarmonicSourceCurrentDensity>
        set_harmonic_source_on_coil_circular(coil_body, harmonic_source_magnitude_used,
                                             harmonic_frequency_hz_runtime, physical_time,
                                             coil_center, coil_source_axis_used, harmonic_phase_runtime);
    SimpleDynamics<electromagnetics::PrescribedComplexSourceCurrentDensity>
        set_frequency_source_on_coil(coil_body, frequency_source_real, frequency_source_imag);
    SimpleDynamics<extra_electromagnetics::CircularFrequencySourceCurrentDensity>
        set_frequency_source_on_coil_circular(coil_body, frequency_source_real_magnitude,
                                              frequency_source_imag_magnitude,
                                              coil_center, coil_source_axis_used);
    SimpleDynamics<ScaleSourceCurrentDensityByRuntimeFactor>
        scale_harmonic_source_on_coil(coil_body, runtime_source_scale);
    SimpleDynamics<ScaleComplexSourceCurrentDensityByRuntimeFactor>
        scale_frequency_source_on_coil(coil_body, runtime_source_scale);

    SimpleDynamics<electromagnetics::UpdateVectorPotentialTimeDerivative> update_a_dot_plate(plate_body);
    SimpleDynamics<electromagnetics::UpdateVectorPotentialTimeDerivative> update_a_dot_coil(coil_body);
    SimpleDynamics<electromagnetics::UpdateVectorPotentialTimeDerivative> update_a_dot_air(air_body);

    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermInner> electric_potential_source_plate(plate_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermInner> electric_potential_source_coil(coil_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermInner> electric_potential_source_air(air_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermContact> electric_potential_source_plate_contact(plate_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermContact> electric_potential_source_coil_contact(coil_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceTermContact> electric_potential_source_air_contact(air_contact);

    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        electric_potential_relaxation_plate_inner(plate_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        electric_potential_relaxation_coil_inner(coil_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationInner>
        electric_potential_relaxation_air_inner(air_inner, 1.0);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationContact>
        electric_potential_relaxation_plate_contact(plate_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationContact>
        electric_potential_relaxation_coil_contact(coil_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialRelaxationContact>
        electric_potential_relaxation_air_contact(air_contact);
    SimpleDynamics<electromagnetics::UpdateElectricPotentialByRelaxationRate>
        electric_potential_relaxation_plate_update(plate_body, 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateElectricPotentialByRelaxationRate>
        electric_potential_relaxation_coil_update(coil_body, 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateElectricPotentialByRelaxationRate>
        electric_potential_relaxation_air_update(air_body, 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_reference_plate(plate_phi_reference_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_boundary_plate(plate_phi_boundary_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_reference_coil(coil_phi_reference_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_boundary_coil(coil_phi_boundary_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_reference_air(air_phi_reference_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_boundary_air(air_phi_boundary_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainElectricPotential>
        constrain_phi_all_air(air_phi_all_region, 0.0);
    SimpleDynamics<electromagnetics::ConstrainVectorPotential>
        constrain_a_reference_plate(plate_a_reference_region, ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorPotential>
        constrain_a_reference_coil(coil_a_reference_region, ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorPotential>
        constrain_a_reference_air(air_a_reference_region, ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorPotential>
        constrain_a_boundary_air_dirichlet(air_a_boundary_region, ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldToAxisAlignedBoxNormalByName>
        constrain_a_boundary_air_magnetic_insulation(air_a_boundary_region, "VectorPotential",
                                                     air_center, 0.5 * air_size);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_plate_inner(plate_inner, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_plate_contact(plate_contact, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_coil_inner(coil_inner, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_coil_contact(coil_contact, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_air_inner(air_inner, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_air_contact(air_contact, "VectorPotential", "VectorPotentialCurl", curl_operator_scaling);

    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_plate_inner(plate_inner, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_plate_contact(plate_contact, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_coil_inner(coil_inner, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_coil_contact(coil_contact, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_air_inner(air_inner, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_air_contact(air_contact, "VectorPotentialCurl", "CurlNuB", curl_operator_scaling);

    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        electric_potential_gradient_plate_inner(plate_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientContact>
        electric_potential_gradient_plate_contact(plate_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        electric_potential_gradient_coil_inner(coil_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientContact>
        electric_potential_gradient_coil_contact(coil_contact);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientInner>
        electric_potential_gradient_air_inner(air_inner);
    InteractionDynamics<electromagnetics::ElectricPotentialGradientContact>
        electric_potential_gradient_air_contact(air_contact);

    InteractionWithUpdate<electromagnetics::VectorPotentialEquationInner>
        vector_potential_equation_plate(plate_inner, plate_a_relaxation_scaling, plate_a_rate_limit);
    InteractionWithUpdate<electromagnetics::VectorPotentialEquationInner>
        vector_potential_equation_coil(coil_inner, coil_a_relaxation_scaling, coil_a_rate_limit);
    InteractionWithUpdate<electromagnetics::VectorPotentialEquationMagneticOnlyInner>
        vector_potential_equation_air(air_inner, air_a_relaxation_scaling, air_a_rate_limit);

    InteractionDynamics<electromagnetics::ElectricFieldCurrentAndJouleHeatInner>
        electric_field_current_heat(plate_inner, false);

    // Frequency-domain (complex A-phi) dynamics.
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_real_plate(plate_inner, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                             harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_real_coil(coil_inner, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                            harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_real_air(air_inner, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                           harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_real_plate_contact(plate_contact, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                                     harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_real_coil_contact(coil_contact, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                                    harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_real_air_contact(air_contact, "VectorPotentialImag", "ElectricPotentialSourceReal",
                                                   harmonic_angular_frequency_runtime, phi_divergence_scaling);

    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_imag_plate(plate_inner, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                             -harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_imag_coil(coil_inner, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                            -harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldInner>
        electric_potential_source_imag_air(air_inner, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                           -harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_imag_plate_contact(plate_contact, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                                     -harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_imag_coil_contact(coil_contact, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                                    -harmonic_angular_frequency_runtime, phi_divergence_scaling);
    InteractionDynamics<electromagnetics::ElectricPotentialSourceFromVectorFieldContact>
        electric_potential_source_imag_air_contact(air_contact, "VectorPotentialReal", "ElectricPotentialSourceImag",
                                                   -harmonic_angular_frequency_runtime, phi_divergence_scaling);

    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_real_relaxation_plate_inner(plate_inner, plate_contact,
                                                       "ElectricPotentialReal", "ElectricPotentialSourceReal",
                                                       "ElectricPotentialChangeRateReal", phi_laplacian_scaling,
                                                       use_frequency_scalar_contact_coupling);
    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_real_relaxation_coil_inner(coil_inner, coil_contact,
                                                      "ElectricPotentialReal", "ElectricPotentialSourceReal",
                                                      "ElectricPotentialChangeRateReal", phi_laplacian_scaling,
                                                      use_frequency_scalar_contact_coupling);
    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_real_relaxation_air_inner(air_inner, air_contact,
                                                     "ElectricPotentialReal", "ElectricPotentialSourceReal",
                                                     "ElectricPotentialChangeRateReal", phi_laplacian_scaling,
                                                     use_frequency_scalar_contact_coupling);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_real_relaxation_plate_update(plate_body, "ElectricPotentialReal", "ElectricPotentialChangeRateReal", 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_real_relaxation_coil_update(coil_body, "ElectricPotentialReal", "ElectricPotentialChangeRateReal", 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_real_relaxation_air_update(air_body, "ElectricPotentialReal", "ElectricPotentialChangeRateReal", 1.0, phi_abs_limit);

    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_imag_relaxation_plate_inner(plate_inner, plate_contact,
                                                       "ElectricPotentialImag", "ElectricPotentialSourceImag",
                                                       "ElectricPotentialChangeRateImag", phi_laplacian_scaling,
                                                       use_frequency_scalar_contact_coupling);
    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_imag_relaxation_coil_inner(coil_inner, coil_contact,
                                                      "ElectricPotentialImag", "ElectricPotentialSourceImag",
                                                      "ElectricPotentialChangeRateImag", phi_laplacian_scaling,
                                                      use_frequency_scalar_contact_coupling);
    InteractionDynamics<electromagnetics::ScalarRelaxationComplexByName>
        electric_potential_imag_relaxation_air_inner(air_inner, air_contact,
                                                     "ElectricPotentialImag", "ElectricPotentialSourceImag",
                                                     "ElectricPotentialChangeRateImag", phi_laplacian_scaling,
                                                     use_frequency_scalar_contact_coupling);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_imag_relaxation_plate_update(plate_body, "ElectricPotentialImag", "ElectricPotentialChangeRateImag", 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_imag_relaxation_coil_update(coil_body, "ElectricPotentialImag", "ElectricPotentialChangeRateImag", 1.0, phi_abs_limit);
    SimpleDynamics<electromagnetics::UpdateScalarByRelaxationRateByName>
        electric_potential_imag_relaxation_air_update(air_body, "ElectricPotentialImag", "ElectricPotentialChangeRateImag", 1.0, phi_abs_limit);

    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_reference_plate(plate_phi_reference_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_boundary_plate(plate_phi_boundary_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_reference_coil(coil_phi_reference_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_boundary_coil(coil_phi_boundary_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_reference_air(air_phi_reference_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_boundary_air(air_phi_boundary_region, "ElectricPotentialReal", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_real_all_air(air_phi_all_region, "ElectricPotentialReal", 0.0);

    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_reference_plate(plate_phi_reference_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_boundary_plate(plate_phi_boundary_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_reference_coil(coil_phi_reference_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_boundary_coil(coil_phi_boundary_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_reference_air(air_phi_reference_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_boundary_air(air_phi_boundary_region, "ElectricPotentialImag", 0.0);
    SimpleDynamics<electromagnetics::ConstrainScalarFieldByName>
        constrain_phi_imag_all_air(air_phi_all_region, "ElectricPotentialImag", 0.0);

    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_real_plate_inner(plate_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                     phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_real_plate_contact(plate_contact, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                       phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_real_coil_inner(coil_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                    phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_real_coil_contact(coil_contact, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                      phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_real_air_inner(air_inner, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                   phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_real_air_contact(air_contact, "ElectricPotentialReal", "ElectricPotentialGradientReal",
                                                     phi_gradient_scaling);

    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_imag_plate_inner(plate_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                     phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_imag_plate_contact(plate_contact, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                       phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_imag_coil_inner(coil_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                    phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_imag_coil_contact(coil_contact, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                      phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientInnerByName>
        electric_potential_gradient_imag_air_inner(air_inner, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                   phi_gradient_scaling);
    InteractionDynamics<electromagnetics::ScalarGradientContactByName>
        electric_potential_gradient_imag_air_contact(air_contact, "ElectricPotentialImag", "ElectricPotentialGradientImag",
                                                     phi_gradient_scaling);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_real_plate_inner(plate_inner, "VectorPotentialReal", "VectorPotentialCurlReal",
                                               curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_real_plate_contact(plate_contact, "VectorPotentialReal", "VectorPotentialCurlReal",
                                                 curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_real_coil_inner(coil_inner, "VectorPotentialReal", "VectorPotentialCurlReal",
                                              curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_real_coil_contact(coil_contact, "VectorPotentialReal", "VectorPotentialCurlReal",
                                                curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_real_air_inner(air_inner, "VectorPotentialReal", "VectorPotentialCurlReal",
                                             curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_real_air_contact(air_contact, "VectorPotentialReal", "VectorPotentialCurlReal",
                                               curl_operator_scaling);

    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_imag_plate_inner(plate_inner, "VectorPotentialImag", "VectorPotentialCurlImag",
                                               curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_imag_plate_contact(plate_contact, "VectorPotentialImag", "VectorPotentialCurlImag",
                                                 curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_imag_coil_inner(coil_inner, "VectorPotentialImag", "VectorPotentialCurlImag",
                                              curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_imag_coil_contact(coil_contact, "VectorPotentialImag", "VectorPotentialCurlImag",
                                                curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlInner>
        vector_potential_curl_imag_air_inner(air_inner, "VectorPotentialImag", "VectorPotentialCurlImag",
                                             curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_imag_air_contact(air_contact, "VectorPotentialImag", "VectorPotentialCurlImag",
                                               curl_operator_scaling);

    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_plate_inner(plate_inner, "VectorPotentialCurlReal", "CurlNuBReal",
                                   curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_plate_contact(plate_contact, "VectorPotentialCurlReal", "CurlNuBReal",
                                     curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_coil_inner(coil_inner, "VectorPotentialCurlReal", "CurlNuBReal",
                                  curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_coil_contact(coil_contact, "VectorPotentialCurlReal", "CurlNuBReal",
                                    curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_real_air_inner(air_inner, "VectorPotentialCurlReal", "CurlNuBReal",
                                 curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_air_contact(air_contact, "VectorPotentialCurlReal", "CurlNuBReal",
                                   curl_operator_scaling);

    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_plate_inner(plate_inner, "VectorPotentialCurlImag", "CurlNuBImag",
                                   curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_plate_contact(plate_contact, "VectorPotentialCurlImag", "CurlNuBImag",
                                     curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_coil_inner(coil_inner, "VectorPotentialCurlImag", "CurlNuBImag",
                                  curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_coil_contact(coil_contact, "VectorPotentialCurlImag", "CurlNuBImag",
                                    curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBInner>
        curl_nu_b_imag_air_inner(air_inner, "VectorPotentialCurlImag", "CurlNuBImag",
                                 curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_air_contact(air_contact, "VectorPotentialCurlImag", "CurlNuBImag",
                                   curl_operator_scaling);

    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_real_plate(plate_inner, plate_contact, harmonic_angular_frequency_runtime, 1.0,
                                             "VectorPotentialReal", "VectorPotentialImag",
                                             "SourceCurrentDensityReal", "ElectricPotentialGradientReal",
                                             "CurlNuBReal", "VectorPotentialChangeRateReal",
                                             freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                             frequency_magnetic_diagonal_scaling, dt_em,
                                             plate_a_relaxation_scaling, plate_a_rate_limit,
                                             frequency_balanced_contact_mag_diagonal_plate_weight,
                                             frequency_contact_mag_diag_ratio_cap_plate,
                                             frequency_contact_mag_diag_adaptive_cap_strength_plate);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_real_coil(coil_inner, coil_contact, harmonic_angular_frequency_runtime, 1.0,
                                            "VectorPotentialReal", "VectorPotentialImag",
                                            "SourceCurrentDensityReal", "ElectricPotentialGradientReal",
                                            "CurlNuBReal", "VectorPotentialChangeRateReal",
                                            freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                            frequency_magnetic_diagonal_scaling, dt_em,
                                            coil_a_relaxation_scaling, coil_a_rate_limit,
                                            frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                            frequency_contact_mag_diag_ratio_cap_coil_air,
                                            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_equation_real_air(air_inner, air_contact,
                                           "VectorPotentialReal",
                                           "SourceCurrentDensityReal",
                                           "CurlNuBReal",
                                           "VectorPotentialChangeRateReal",
                                           frequency_magnetic_diagonal_scaling, dt_em,
                                           air_a_relaxation_scaling, air_a_rate_limit,
                                           frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                           frequency_contact_mag_diag_ratio_cap_coil_air,
                                           frequency_contact_mag_diag_adaptive_cap_strength_coil_air);

    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_imag_plate(plate_inner, plate_contact, harmonic_angular_frequency_runtime, -1.0,
                                             "VectorPotentialImag", "VectorPotentialReal",
                                             "SourceCurrentDensityImag", "ElectricPotentialGradientImag",
                                             "CurlNuBImag", "VectorPotentialChangeRateImag",
                                             freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                             frequency_magnetic_diagonal_scaling, dt_em,
                                             plate_a_relaxation_scaling, plate_a_rate_limit,
                                             frequency_balanced_contact_mag_diagonal_plate_weight,
                                             frequency_contact_mag_diag_ratio_cap_plate,
                                             frequency_contact_mag_diag_adaptive_cap_strength_plate);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_imag_coil(coil_inner, coil_contact, harmonic_angular_frequency_runtime, -1.0,
                                            "VectorPotentialImag", "VectorPotentialReal",
                                            "SourceCurrentDensityImag", "ElectricPotentialGradientImag",
                                            "CurlNuBImag", "VectorPotentialChangeRateImag",
                                            freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                            frequency_magnetic_diagonal_scaling, dt_em,
                                            coil_a_relaxation_scaling, coil_a_rate_limit,
                                            frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                            frequency_contact_mag_diag_ratio_cap_coil_air,
                                            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyEquationComplex>
        vector_potential_magnetic_only_equation_real_coil(coil_inner, coil_contact,
                                                          "VectorPotentialReal",
                                                          "SourceCurrentDensityReal",
                                                          "CurlNuBReal",
                                                          "VectorPotentialChangeRateReal",
                                                          frequency_magnetic_diagonal_scaling, dt_em,
                                                          coil_a_relaxation_scaling, coil_a_rate_limit,
                                                          frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                                          frequency_contact_mag_diag_ratio_cap_coil_air,
                                                          frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_magnetic_only_block_equation_real_coil(coil_inner, coil_contact,
                                                                "VectorPotentialReal",
                                                                "SourceCurrentDensityReal",
                                                                "CurlNuBReal",
                                                                "VectorPotentialChangeRateReal",
                                                                frequency_magnetic_diagonal_scaling, dt_em,
                                                                coil_a_relaxation_scaling, coil_a_rate_limit,
                                                                frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                                                frequency_contact_mag_diag_ratio_cap_coil_air,
                                                                frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyEquationComplex>
        vector_potential_magnetic_only_equation_imag_coil(coil_inner, coil_contact,
                                                          "VectorPotentialImag",
                                                          "SourceCurrentDensityImag",
                                                          "CurlNuBImag",
                                                          "VectorPotentialChangeRateImag",
                                                          frequency_magnetic_diagonal_scaling, dt_em,
                                                          coil_a_relaxation_scaling, coil_a_rate_limit,
                                                          frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                                          frequency_contact_mag_diag_ratio_cap_coil_air,
                                                          frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_magnetic_only_block_equation_imag_coil(coil_inner, coil_contact,
                                                                "VectorPotentialImag",
                                                                "SourceCurrentDensityImag",
                                                                "CurlNuBImag",
                                                                "VectorPotentialChangeRateImag",
                                                                frequency_magnetic_diagonal_scaling, dt_em,
                                                                coil_a_relaxation_scaling, coil_a_rate_limit,
                                                                frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                                                frequency_contact_mag_diag_ratio_cap_coil_air,
                                                                frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyCoupledBlockEquationComplex>
        vector_potential_equation_coupled_plate(plate_inner, plate_contact, harmonic_angular_frequency_runtime,
                                                "VectorPotentialReal", "VectorPotentialImag",
                                                "SourceCurrentDensityReal", "SourceCurrentDensityImag",
                                                "ElectricPotentialGradientReal", "ElectricPotentialGradientImag",
                                                "CurlNuBReal", "CurlNuBImag",
                                                "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
                                                freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                                frequency_magnetic_diagonal_scaling, dt_em,
                                                plate_a_relaxation_scaling, plate_a_rate_limit,
                                                frequency_balanced_contact_mag_diagonal_plate_weight,
                                                frequency_contact_mag_diag_ratio_cap_plate,
                                                frequency_contact_mag_diag_adaptive_cap_strength_plate);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyCoupledEquationComplex>
        vector_potential_equation_coupled_coil(coil_inner, coil_contact, harmonic_angular_frequency_runtime,
                                               "VectorPotentialReal", "VectorPotentialImag",
                                               "SourceCurrentDensityReal", "SourceCurrentDensityImag",
                                               "ElectricPotentialGradientReal", "ElectricPotentialGradientImag",
                                               "CurlNuBReal", "CurlNuBImag",
                                               "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
                                               freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
                                               frequency_magnetic_diagonal_scaling, dt_em,
                                               coil_a_relaxation_scaling, coil_a_rate_limit,
                                               frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                               frequency_contact_mag_diag_ratio_cap_coil_air,
                                               frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_equation_imag_air(air_inner, air_contact,
                                           "VectorPotentialImag",
                                           "SourceCurrentDensityImag",
                                           "CurlNuBImag",
                                           "VectorPotentialChangeRateImag",
                                           frequency_magnetic_diagonal_scaling, dt_em,
                                           air_a_relaxation_scaling, air_a_rate_limit,
                                           frequency_balanced_contact_mag_diagonal_coil_air_weight,
                                           frequency_contact_mag_diag_ratio_cap_coil_air,
                                           frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_real_coil_contact_air_only(
            coil_contact_air_only, "VectorPotentialReal", "VectorPotentialCurlReal",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_imag_coil_contact_air_only(
            coil_contact_air_only, "VectorPotentialImag", "VectorPotentialCurlImag",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_coil_contact_air_only(
            coil_contact_air_only, "VectorPotentialCurlReal", "CurlNuBReal",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_coil_contact_air_only(
            coil_contact_air_only, "VectorPotentialCurlImag", "CurlNuBImag",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_real_air_contact_coil_only(
            air_contact_coil_only, "VectorPotentialReal", "VectorPotentialCurlReal",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::VectorPotentialCurlContact>
        vector_potential_curl_imag_air_contact_coil_only(
            air_contact_coil_only, "VectorPotentialImag", "VectorPotentialCurlImag",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_real_air_contact_coil_only(
            air_contact_coil_only, "VectorPotentialCurlReal", "CurlNuBReal",
            curl_operator_scaling);
    InteractionDynamics<electromagnetics::CurlNuBContact>
        curl_nu_b_imag_air_contact_coil_only(
            air_contact_coil_only, "VectorPotentialCurlImag", "CurlNuBImag",
            curl_operator_scaling);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_real_coil_air_only(
            coil_inner, coil_contact_air_only, harmonic_angular_frequency_runtime, 1.0,
            "VectorPotentialReal", "VectorPotentialImag",
            "SourceCurrentDensityReal", "ElectricPotentialGradientReal",
            "CurlNuBReal", "VectorPotentialChangeRateReal",
            freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyEquationComplex>
        vector_potential_equation_imag_coil_air_only(
            coil_inner, coil_contact_air_only, harmonic_angular_frequency_runtime, -1.0,
            "VectorPotentialImag", "VectorPotentialReal",
            "SourceCurrentDensityImag", "ElectricPotentialGradientImag",
            "CurlNuBImag", "VectorPotentialChangeRateImag",
            freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyEquationComplex>
        vector_potential_magnetic_only_equation_real_coil_air_only(
            coil_inner, coil_contact_air_only,
            "VectorPotentialReal",
            "SourceCurrentDensityReal",
            "CurlNuBReal",
            "VectorPotentialChangeRateReal",
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_magnetic_only_block_equation_real_coil_air_only(
            coil_inner, coil_contact_air_only,
            "VectorPotentialReal",
            "SourceCurrentDensityReal",
            "CurlNuBReal",
            "VectorPotentialChangeRateReal",
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyEquationComplex>
        vector_potential_magnetic_only_equation_imag_coil_air_only(
            coil_inner, coil_contact_air_only,
            "VectorPotentialImag",
            "SourceCurrentDensityImag",
            "CurlNuBImag",
            "VectorPotentialChangeRateImag",
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_magnetic_only_block_equation_imag_coil_air_only(
            coil_inner, coil_contact_air_only,
            "VectorPotentialImag",
            "SourceCurrentDensityImag",
            "CurlNuBImag",
            "VectorPotentialChangeRateImag",
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyCoupledEquationComplex>
        vector_potential_equation_coupled_coil_air_only(
            coil_inner, coil_contact_air_only, harmonic_angular_frequency_runtime,
            "VectorPotentialReal", "VectorPotentialImag",
            "SourceCurrentDensityReal", "SourceCurrentDensityImag",
            "ElectricPotentialGradientReal", "ElectricPotentialGradientImag",
            "CurlNuBReal", "CurlNuBImag",
            "VectorPotentialChangeRateReal", "VectorPotentialChangeRateImag",
            freq_sigma_relaxation_scaling, freq_sigma_relaxation_floor,
            frequency_magnetic_diagonal_scaling, dt_em,
            coil_a_relaxation_scaling, coil_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_equation_real_air_coil_only(
            air_inner, air_contact_coil_only,
            "VectorPotentialReal",
            "SourceCurrentDensityReal",
            "CurlNuBReal",
            "VectorPotentialChangeRateReal",
            frequency_magnetic_diagonal_scaling, dt_em,
            air_a_relaxation_scaling, air_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);
    InteractionWithUpdate<electromagnetics::VectorPotentialFrequencyMagneticOnlyBlockEquationComplex>
        vector_potential_equation_imag_air_coil_only(
            air_inner, air_contact_coil_only,
            "VectorPotentialImag",
            "SourceCurrentDensityImag",
            "CurlNuBImag",
            "VectorPotentialChangeRateImag",
            frequency_magnetic_diagonal_scaling, dt_em,
            air_a_relaxation_scaling, air_a_rate_limit,
            frequency_balanced_contact_mag_diagonal_coil_air_weight,
            frequency_contact_mag_diag_ratio_cap_coil_air,
            frequency_contact_mag_diag_adaptive_cap_strength_coil_air);

    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_real_reference_plate(plate_a_reference_region, "VectorPotentialReal", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_real_reference_coil(coil_a_reference_region, "VectorPotentialReal", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_real_reference_air(air_a_reference_region, "VectorPotentialReal", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_real_boundary_air_dirichlet(air_a_boundary_region, "VectorPotentialReal", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldToAxisAlignedBoxNormalByName>
        constrain_a_real_boundary_air_magnetic_insulation(air_a_boundary_region, "VectorPotentialReal",
                                                          air_center, 0.5 * air_size);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_imag_reference_plate(plate_a_reference_region, "VectorPotentialImag", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_imag_reference_coil(coil_a_reference_region, "VectorPotentialImag", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_imag_reference_air(air_a_reference_region, "VectorPotentialImag", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldByName>
        constrain_a_imag_boundary_air_dirichlet(air_a_boundary_region, "VectorPotentialImag", ZeroData<Vecd>::value);
    SimpleDynamics<electromagnetics::ConstrainVectorFieldToAxisAlignedBoxNormalByName>
        constrain_a_imag_boundary_air_magnetic_insulation(air_a_boundary_region, "VectorPotentialImag",
                                                          air_center, 0.5 * air_size);

    InteractionDynamics<electromagnetics::FrequencyElectricFieldCurrentAndJouleHeatInner>
        electric_field_current_heat_frequency(plate_inner, harmonic_angular_frequency_runtime);

    //------------------------------------------------------------------
    //  Thermal one-way coupling setup
    //------------------------------------------------------------------
    GetDiffusionTimeStepSize get_thermal_time_step(plate_body);
    ThermalRelaxationInner thermal_relaxation(plate_inner);
    SimpleDynamics<AddJouleHeatToTemperature> add_joule_heat_to_temperature(plate_body);

    auto apply_time_domain_air_a_boundary_constraint = [&]()
    {
        if (!constrain_a_boundary_air_enabled)
        {
            return;
        }
        if (use_frequency_air_magnetic_insulation_boundary)
        {
            constrain_a_boundary_air_magnetic_insulation.exec();
        }
        else
        {
            constrain_a_boundary_air_dirichlet.exec();
        }
    };

    auto apply_frequency_air_a_boundary_constraint = [&]()
    {
        if (!constrain_a_boundary_air_enabled)
        {
            return;
        }
        if (use_frequency_air_magnetic_insulation_boundary)
        {
            constrain_a_real_boundary_air_magnetic_insulation.exec();
            constrain_a_imag_boundary_air_magnetic_insulation.exec();
        }
        else
        {
            constrain_a_real_boundary_air_dirichlet.exec();
            constrain_a_imag_boundary_air_dirichlet.exec();
        }
    };

    //------------------------------------------------------------------
    //  Output
    //------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(plate_body, "Temperature");
    write_states.addToWrite<Real>(plate_body, "ElectricPotential");
    write_states.addToWrite<Real>(plate_body, "JouleHeatSource");
    write_states.addToWrite<Vecd>(plate_body, "VectorPotential");
    write_states.addToWrite<AngularVecd>(plate_body, "VectorPotentialCurl");
    write_states.addToWrite<Vecd>(plate_body, "VectorPotentialReal");
    write_states.addToWrite<Vecd>(plate_body, "VectorPotentialImag");
    write_states.addToWrite<AngularVecd>(plate_body, "VectorPotentialCurlReal");
    write_states.addToWrite<AngularVecd>(plate_body, "VectorPotentialCurlImag");
    write_states.addToWrite<Vecd>(plate_body, "CurrentDensity");
    write_states.addToWrite<Vecd>(plate_body, "ElectricField");
    write_states.addToWrite<Vecd>(plate_body, "CurrentDensityReal");
    write_states.addToWrite<Vecd>(plate_body, "CurrentDensityImag");
    write_states.addToWrite<Vecd>(plate_body, "ElectricFieldReal");
    write_states.addToWrite<Vecd>(plate_body, "ElectricFieldImag");
    write_states.addToWrite<Vecd>(coil_body, "VectorPotential");
    write_states.addToWrite<AngularVecd>(coil_body, "VectorPotentialCurl");
    write_states.addToWrite<Vecd>(coil_body, "VectorPotentialReal");
    write_states.addToWrite<Vecd>(coil_body, "VectorPotentialImag");
    write_states.addToWrite<AngularVecd>(coil_body, "VectorPotentialCurlReal");
    write_states.addToWrite<AngularVecd>(coil_body, "VectorPotentialCurlImag");
    write_states.addToWrite<Real>(coil_body, "ElectricPotential");
    write_states.addToWrite<Real>(coil_body, "ElectricPotentialReal");
    write_states.addToWrite<Real>(coil_body, "ElectricPotentialImag");

    //------------------------------------------------------------------
    //  Initialization
    //------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    temperature_initialization.exec();
    initialize_aphi_plate.exec();
    initialize_aphi_coil.exec();
    initialize_aphi_air.exec();
    if (use_frequency_aphi)
    {
        initialize_aphi_frequency_plate.exec();
        initialize_aphi_frequency_coil.exec();
        initialize_aphi_frequency_air.exec();
    }
    set_em_material_plate.exec();
    set_em_material_coil.exec();
    set_em_material_air.exec();
    if (use_frequency_aphi)
    {
        if (use_circular_coil_source)
        {
            set_frequency_source_on_coil_circular.exec();
        }
        else
        {
            set_frequency_source_on_coil.exec();
        }
        if (apply_runtime_source_scale)
        {
            scale_frequency_source_on_coil.exec();
        }

        constrain_phi_real_reference_plate.exec();
        constrain_phi_imag_reference_plate.exec();
        if (constrain_phi_boundary_plate_enabled)
        {
            constrain_phi_real_boundary_plate.exec();
            constrain_phi_imag_boundary_plate.exec();
        }
        constrain_phi_real_reference_coil.exec();
        constrain_phi_imag_reference_coil.exec();
        if (constrain_phi_boundary_coil_enabled)
        {
            constrain_phi_real_boundary_coil.exec();
            constrain_phi_imag_boundary_coil.exec();
        }
        if (use_frequency_air_scalar_potential)
        {
            constrain_phi_real_reference_air.exec();
            constrain_phi_imag_reference_air.exec();
            if (constrain_phi_boundary_air_enabled)
            {
                constrain_phi_real_boundary_air.exec();
                constrain_phi_imag_boundary_air.exec();
            }
            if (constrain_phi_all_air_enabled)
            {
                constrain_phi_real_all_air.exec();
                constrain_phi_imag_all_air.exec();
            }
        }

        constrain_a_real_reference_plate.exec();
        constrain_a_real_reference_coil.exec();
        constrain_a_real_reference_air.exec();
        apply_frequency_air_a_boundary_constraint();
        constrain_a_imag_reference_plate.exec();
        constrain_a_imag_reference_coil.exec();
        constrain_a_imag_reference_air.exec();
    }
    else
    {
        if (use_circular_coil_source)
        {
            set_harmonic_source_on_coil_circular.exec();
        }
        else
        {
            set_harmonic_source_on_coil.exec();
        }
        if (apply_runtime_source_scale)
        {
            scale_harmonic_source_on_coil.exec();
        }
        constrain_phi_reference_plate.exec();
        if (constrain_phi_boundary_plate_enabled)
        {
            constrain_phi_boundary_plate.exec();
        }
        constrain_phi_reference_coil.exec();
        if (constrain_phi_boundary_coil_enabled)
        {
            constrain_phi_boundary_coil.exec();
        }
        constrain_phi_reference_air.exec();
        if (constrain_phi_boundary_air_enabled)
        {
            constrain_phi_boundary_air.exec();
        }
        if (constrain_phi_all_air_enabled)
        {
            constrain_phi_all_air.exec();
        }
        constrain_a_reference_plate.exec();
        constrain_a_reference_coil.exec();
        constrain_a_reference_air.exec();
        apply_time_domain_air_a_boundary_constraint();
    }

    BaseParticles &plate_particles = plate_body.getBaseParticles();
    Real *plate_vol = plate_particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
    Real *plate_temperature = plate_particles.getVariableDataByName<Real>("Temperature");
    Real *plate_electric_potential = plate_particles.getVariableDataByName<Real>("ElectricPotential");
    Real *plate_joule_heat = plate_particles.getVariableDataByName<Real>("JouleHeatSource");
    Vecd *plate_vector_potential = plate_particles.getVariableDataByName<Vecd>("VectorPotential");
    Vecd *plate_vector_potential_dt = plate_particles.getVariableDataByName<Vecd>("VectorPotentialTimeDerivative");
    Vecd *plate_vector_potential_change_rate = plate_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRate");
    AngularVecd *plate_vector_potential_curl = plate_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *plate_electric_potential_gradient = plate_particles.getVariableDataByName<Vecd>("ElectricPotentialGradient");
    Vecd *plate_electric_field = plate_particles.getVariableDataByName<Vecd>("ElectricField");
    Vecd *plate_current_density = plate_particles.getVariableDataByName<Vecd>("CurrentDensity");
    Real *plate_operator_verify_ax_real = plate_particles.getVariableDataByName<Real>("OperatorVerifyAxReal");
    Vecd *plate_operator_verify_ax_real_gradient = plate_particles.getVariableDataByName<Vecd>("OperatorVerifyAxRealGradient");
    VecMatd *plate_operator_verify_ax_real_hessian = plate_particles.getVariableDataByName<VecMatd>("OperatorVerifyAxRealHessian");
    Real *plate_operator_verify_ay_real = plate_particles.getVariableDataByName<Real>("OperatorVerifyAyReal");
    Vecd *plate_operator_verify_ay_real_gradient = plate_particles.getVariableDataByName<Vecd>("OperatorVerifyAyRealGradient");
    VecMatd *plate_operator_verify_ay_real_hessian = plate_particles.getVariableDataByName<VecMatd>("OperatorVerifyAyRealHessian");
    Real *plate_operator_verify_az_real = plate_particles.getVariableDataByName<Real>("OperatorVerifyAzReal");
    Vecd *plate_operator_verify_az_real_gradient = plate_particles.getVariableDataByName<Vecd>("OperatorVerifyAzRealGradient");
    VecMatd *plate_operator_verify_az_real_hessian = plate_particles.getVariableDataByName<VecMatd>("OperatorVerifyAzRealHessian");
    Vecd *plate_operator_verify_curl_nu_b_from_component_hessians_real =
        plate_particles.getVariableDataByName<Vecd>("OperatorVerifyCurlNuBFromComponentHessiansReal");
    Vecd *plate_vector_potential_real = nullptr;
    Vecd *plate_vector_potential_imag = nullptr;
    Vecd *plate_vector_potential_change_rate_real = nullptr;
    Vecd *plate_vector_potential_change_rate_imag = nullptr;
    Vecd *plate_curl_nu_b_real = nullptr;
    Vecd *plate_curl_nu_b_imag = nullptr;
    AngularVecd *plate_vector_potential_curl_real = nullptr;
    AngularVecd *plate_vector_potential_curl_imag = nullptr;
    Vecd *plate_electric_potential_gradient_real = nullptr;
    Vecd *plate_electric_potential_gradient_imag = nullptr;
    Vecd *plate_electric_field_real = nullptr;
    Vecd *plate_electric_field_imag = nullptr;
    Vecd *plate_current_density_real = nullptr;
    Vecd *plate_current_density_imag = nullptr;
    Real *plate_electric_potential_real = nullptr;
    Real *plate_electric_potential_imag = nullptr;

    BaseParticles &coil_particles = coil_body.getBaseParticles();
    Real *coil_vol = coil_particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *coil_particle_positions = coil_particles.getVariableDataByName<Vecd>("Position");
    Vecd *coil_vector_potential = coil_particles.getVariableDataByName<Vecd>("VectorPotential");
    Vecd *coil_vector_potential_dt = coil_particles.getVariableDataByName<Vecd>("VectorPotentialTimeDerivative");
    Vecd *coil_vector_potential_change_rate = coil_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRate");
    Vecd *coil_vector_potential_real = nullptr;
    Vecd *coil_vector_potential_imag = nullptr;
    Vecd *coil_vector_potential_change_rate_real = nullptr;
    Vecd *coil_vector_potential_change_rate_imag = nullptr;
    AngularVecd *coil_vector_potential_curl_real = nullptr;
    AngularVecd *coil_vector_potential_curl_imag = nullptr;
    Vecd *coil_source_current_density = coil_particles.getVariableDataByName<Vecd>("SourceCurrentDensity");
    Vecd *coil_electric_potential_gradient = coil_particles.getVariableDataByName<Vecd>("ElectricPotentialGradient");
    Vecd *coil_curl_nu_b = coil_particles.getVariableDataByName<Vecd>("CurlNuB");
    Vecd *coil_electric_potential_gradient_real = nullptr;
    Vecd *coil_electric_potential_gradient_imag = nullptr;
    Vecd *coil_curl_nu_b_real = nullptr;
    Vecd *coil_curl_nu_b_imag = nullptr;
    Real *coil_electric_potential = coil_particles.getVariableDataByName<Real>("ElectricPotential");
    Vecd *coil_source_current_density_real = nullptr;
    Vecd *coil_source_current_density_imag = nullptr;
    Real *coil_electric_potential_real = nullptr;
    Real *coil_electric_potential_imag = nullptr;

    BaseParticles &air_particles = air_body.getBaseParticles();
    Real *air_vol = air_particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *air_electric_potential = air_particles.getVariableDataByName<Real>("ElectricPotential");
    Vecd *air_vector_potential = air_particles.getVariableDataByName<Vecd>("VectorPotential");
    Vecd *air_vector_potential_dt = air_particles.getVariableDataByName<Vecd>("VectorPotentialTimeDerivative");
    Vecd *air_vector_potential_change_rate = air_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRate");
    Vecd *air_vector_potential_real = nullptr;
    Vecd *air_vector_potential_imag = nullptr;
    Vecd *air_vector_potential_change_rate_real = nullptr;
    Vecd *air_vector_potential_change_rate_imag = nullptr;
    Vecd *air_curl_nu_b_real = nullptr;
    Vecd *air_curl_nu_b_imag = nullptr;
    AngularVecd *air_vector_potential_curl = air_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurl");
    Vecd *air_electric_potential_gradient = air_particles.getVariableDataByName<Vecd>("ElectricPotentialGradient");
    AngularVecd *air_vector_potential_curl_real = nullptr;
    AngularVecd *air_vector_potential_curl_imag = nullptr;
    Real *air_electric_potential_real = nullptr;
    Real *air_electric_potential_imag = nullptr;

    if (use_frequency_aphi)
    {
        plate_vector_potential_real = plate_particles.getVariableDataByName<Vecd>("VectorPotentialReal");
        plate_vector_potential_imag = plate_particles.getVariableDataByName<Vecd>("VectorPotentialImag");
        plate_vector_potential_change_rate_real = plate_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
        plate_vector_potential_change_rate_imag = plate_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
        plate_curl_nu_b_real = plate_particles.getVariableDataByName<Vecd>("CurlNuBReal");
        plate_curl_nu_b_imag = plate_particles.getVariableDataByName<Vecd>("CurlNuBImag");
        plate_vector_potential_curl_real = plate_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
        plate_vector_potential_curl_imag = plate_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
        plate_electric_potential_gradient_real = plate_particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
        plate_electric_potential_gradient_imag = plate_particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
        plate_electric_field_real = plate_particles.getVariableDataByName<Vecd>("ElectricFieldReal");
        plate_electric_field_imag = plate_particles.getVariableDataByName<Vecd>("ElectricFieldImag");
        plate_current_density_real = plate_particles.getVariableDataByName<Vecd>("CurrentDensityReal");
        plate_current_density_imag = plate_particles.getVariableDataByName<Vecd>("CurrentDensityImag");
        plate_electric_potential_real = plate_particles.getVariableDataByName<Real>("ElectricPotentialReal");
        plate_electric_potential_imag = plate_particles.getVariableDataByName<Real>("ElectricPotentialImag");

        coil_vector_potential_real = coil_particles.getVariableDataByName<Vecd>("VectorPotentialReal");
        coil_vector_potential_imag = coil_particles.getVariableDataByName<Vecd>("VectorPotentialImag");
        coil_source_current_density_real = coil_particles.getVariableDataByName<Vecd>("SourceCurrentDensityReal");
        coil_source_current_density_imag = coil_particles.getVariableDataByName<Vecd>("SourceCurrentDensityImag");
        coil_vector_potential_change_rate_real = coil_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
        coil_vector_potential_change_rate_imag = coil_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
        coil_vector_potential_curl_real = coil_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
        coil_vector_potential_curl_imag = coil_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
        coil_electric_potential_gradient_real = coil_particles.getVariableDataByName<Vecd>("ElectricPotentialGradientReal");
        coil_electric_potential_gradient_imag = coil_particles.getVariableDataByName<Vecd>("ElectricPotentialGradientImag");
        coil_curl_nu_b_real = coil_particles.getVariableDataByName<Vecd>("CurlNuBReal");
        coil_curl_nu_b_imag = coil_particles.getVariableDataByName<Vecd>("CurlNuBImag");
        coil_electric_potential_real = coil_particles.getVariableDataByName<Real>("ElectricPotentialReal");
        coil_electric_potential_imag = coil_particles.getVariableDataByName<Real>("ElectricPotentialImag");

        air_vector_potential_real = air_particles.getVariableDataByName<Vecd>("VectorPotentialReal");
        air_vector_potential_imag = air_particles.getVariableDataByName<Vecd>("VectorPotentialImag");
        air_vector_potential_curl_real = air_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlReal");
        air_vector_potential_curl_imag = air_particles.getVariableDataByName<AngularVecd>("VectorPotentialCurlImag");
        air_vector_potential_change_rate_real = air_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateReal");
        air_vector_potential_change_rate_imag = air_particles.getVariableDataByName<Vecd>("VectorPotentialChangeRateImag");
        air_curl_nu_b_real = air_particles.getVariableDataByName<Vecd>("CurlNuBReal");
        air_curl_nu_b_imag = air_particles.getVariableDataByName<Vecd>("CurlNuBImag");
        air_electric_potential_real = air_particles.getVariableDataByName<Real>("ElectricPotentialReal");
        air_electric_potential_imag = air_particles.getVariableDataByName<Real>("ElectricPotentialImag");
    }

    AngularVecd *plate_vector_potential_curl_probe = use_frequency_aphi ? plate_vector_potential_curl_real : plate_vector_potential_curl;
    AngularVecd *air_vector_potential_curl_probe = use_frequency_aphi ? air_vector_potential_curl_real : air_vector_potential_curl;
    const size_t total_plate_particles_count = plate_particles.TotalRealParticles();
    const size_t total_coil_particles_count = coil_particles.TotalRealParticles();
    const size_t total_air_particles_count = air_particles.TotalRealParticles();
    Real *coil_magnetic_reluctivity =
        coil_particles.getVariableDataByName<Real>("MagneticReluctivity");
    Real *air_magnetic_reluctivity =
        air_particles.getVariableDataByName<Real>("MagneticReluctivity");
    std::vector<Real> coil_magnetic_reluctivity_base(
        total_coil_particles_count, magnetic_reluctivity);
    std::vector<Real> air_magnetic_reluctivity_base(
        total_air_particles_count, magnetic_reluctivity);
    Real vacuum_reluctivity_scale_runtime = vacuum_reluctivity_scale_initial;
    bool apply_vacuum_reluctivity_scale_runtime =
        enable_coil_air_vacuum_mode &&
        use_frequency_aphi &&
        (enable_vacuum_reluctivity_feedback ||
         fabs(vacuum_reluctivity_scale_runtime - static_cast<Real>(1.0)) > TinyReal);
    if (coil_magnetic_reluctivity != nullptr)
    {
        for (size_t i = 0; i != total_coil_particles_count; ++i)
        {
            coil_magnetic_reluctivity_base[i] = coil_magnetic_reluctivity[i];
        }
    }
    if (air_magnetic_reluctivity != nullptr)
    {
        for (size_t i = 0; i != total_air_particles_count; ++i)
        {
            air_magnetic_reluctivity_base[i] = air_magnetic_reluctivity[i];
        }
    }
    auto apply_vacuum_reluctivity_scale = [&](Real reluctivity_scale)
    {
        if (!apply_vacuum_reluctivity_scale_runtime)
        {
            return;
        }
        Real clamped_scale =
            SMIN(vacuum_reluctivity_scale_max,
                 SMAX(vacuum_reluctivity_scale_min, reluctivity_scale));
        if (coil_magnetic_reluctivity != nullptr)
        {
            for (size_t i = 0; i != total_coil_particles_count; ++i)
            {
                coil_magnetic_reluctivity[i] =
                    coil_magnetic_reluctivity_base[i] * clamped_scale;
            }
        }
        if (air_magnetic_reluctivity != nullptr)
        {
            for (size_t i = 0; i != total_air_particles_count; ++i)
            {
                air_magnetic_reluctivity[i] =
                    air_magnetic_reluctivity_base[i] * clamped_scale;
            }
        }
    };
    if (apply_vacuum_reluctivity_scale_runtime)
    {
        apply_vacuum_reluctivity_scale(vacuum_reluctivity_scale_runtime);
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-vacuum-reluctivity] apply initial scale="
                  << vacuum_reluctivity_scale_runtime
                  << ", scale_min=" << vacuum_reluctivity_scale_min
                  << ", scale_max=" << vacuum_reluctivity_scale_max
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    }
    std::vector<Vecd> coil_curl_nu_b_real_direct_backup;
    std::vector<Vecd> coil_curl_nu_b_imag_direct_backup;
    std::vector<Vecd> coil_curl_nu_b_real_contact_backup;
    std::vector<Vecd> coil_curl_nu_b_imag_contact_backup;
    std::vector<Vecd> air_curl_nu_b_real_inner_backup;
    std::vector<Vecd> air_curl_nu_b_imag_inner_backup;
    Real air_curl_inner_avg_norm_last = 0.0;
    Real air_curl_contact_avg_norm_last = 0.0;
    Real air_curl_inner_max_norm_last = 0.0;
    Real air_curl_contact_max_norm_last = 0.0;
    if (use_frequency_aphi && use_frequency_coil_component_hessian_inner)
    {
        coil_curl_nu_b_real_direct_backup.resize(total_coil_particles_count);
        coil_curl_nu_b_imag_direct_backup.resize(total_coil_particles_count);
        coil_curl_nu_b_real_contact_backup.resize(total_coil_particles_count);
        coil_curl_nu_b_imag_contact_backup.resize(total_coil_particles_count);
    }
    if (use_frequency_aphi)
    {
        air_curl_nu_b_real_inner_backup.resize(total_air_particles_count);
        air_curl_nu_b_imag_inner_backup.resize(total_air_particles_count);
    }
    std::vector<Vecd> plate_vector_potential_real_backup;
    std::vector<Vecd> plate_vector_potential_imag_backup;
    std::vector<Real> plate_electric_potential_real_backup;
    std::vector<Real> plate_electric_potential_imag_backup;
    std::vector<Vecd> plate_vector_potential_change_rate_real_backup;
    std::vector<Vecd> plate_vector_potential_change_rate_imag_backup;
    if (use_frequency_aphi &&
        (use_frequency_plate_backtracking ||
         (use_frequency_adaptive_plate_block_sweeps && frequency_plate_block_sweeps > 1)))
    {
        plate_vector_potential_real_backup.resize(total_plate_particles_count);
        plate_vector_potential_imag_backup.resize(total_plate_particles_count);
        plate_electric_potential_real_backup.resize(total_plate_particles_count);
        plate_electric_potential_imag_backup.resize(total_plate_particles_count);
        plate_vector_potential_change_rate_real_backup.resize(total_plate_particles_count);
        plate_vector_potential_change_rate_imag_backup.resize(total_plate_particles_count);
    }
    FrequencyEmStateBackup frequency_em_state_backup;
    bool use_frequency_em_state_backup =
        use_frequency_aphi &&
        (use_frequency_global_backtracking || use_frequency_hard_guard);
    if (use_frequency_em_state_backup)
    {
        frequency_em_state_backup.plate_a_real.resize(total_plate_particles_count);
        frequency_em_state_backup.plate_a_imag.resize(total_plate_particles_count);
        frequency_em_state_backup.plate_phi_real.resize(total_plate_particles_count);
        frequency_em_state_backup.plate_phi_imag.resize(total_plate_particles_count);
        frequency_em_state_backup.coil_a_real.resize(coil_particles.TotalRealParticles());
        frequency_em_state_backup.coil_a_imag.resize(coil_particles.TotalRealParticles());
        if (use_frequency_coil_scalar_potential)
        {
            frequency_em_state_backup.coil_phi_real.resize(coil_particles.TotalRealParticles());
            frequency_em_state_backup.coil_phi_imag.resize(coil_particles.TotalRealParticles());
        }
        frequency_em_state_backup.air_a_real.resize(total_air_particles_count);
        frequency_em_state_backup.air_a_imag.resize(total_air_particles_count);
        if (use_frequency_air_scalar_potential)
        {
            frequency_em_state_backup.air_phi_real.resize(total_air_particles_count);
            frequency_em_state_backup.air_phi_imag.resize(total_air_particles_count);
        }
    }

    auto find_nearest_particle_index = [](BaseParticles &particles, Vecd *positions, const Vec3d &target_position) -> size_t
    {
        size_t total_real_particles = particles.TotalRealParticles();
        if (total_real_particles == 0)
        {
            return 0;
        }
        size_t nearest_index = 0;
        Real min_distance_square = std::numeric_limits<Real>::max();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Real distance_square = (positions[i] - target_position).squaredNorm();
            if (distance_square < min_distance_square)
            {
                min_distance_square = distance_square;
                nearest_index = i;
            }
        }
        return nearest_index;
    };

    auto find_nearest_plate_particle_for_j_curve = [&](const Vec3d &target_position) -> size_t
    {
        const size_t total_real_particles = plate_particles.TotalRealParticles();
        if (total_real_particles == 0)
        {
            return 0;
        }
        Real dz_limit = j_curve_plate_dz_band_start_mm >= TinyReal ? j_curve_plate_dz_band_start_mm
                                                                    : dp_plate;
        const Real dz_cap = std::max(j_curve_plate_dz_band_cap_mm, dz_limit);
        while (dz_limit <= dz_cap + TinyReal)
        {
            size_t best_index = 0;
            Real best_xy_sq = std::numeric_limits<Real>::max();
            bool any_in_band = false;
            for (size_t i = 0; i != total_real_particles; ++i)
            {
                const Vecd delta = plate_positions[i] - target_position;
                if (fabs(delta[2]) > dz_limit)
                {
                    continue;
                }
                any_in_band = true;
                const Real xy_sq = delta[0] * delta[0] + delta[1] * delta[1];
                if (xy_sq < best_xy_sq)
                {
                    best_xy_sq = xy_sq;
                    best_index = i;
                }
            }
            if (any_in_band)
            {
                return best_index;
            }
            dz_limit *= 2.0f;
        }
        size_t nearest_index = 0;
        Real min_score = std::numeric_limits<Real>::max();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            const Vecd delta = plate_positions[i] - target_position;
            const Real score = delta[0] * delta[0] + delta[1] * delta[1] +
                               j_curve_plate_z_distance_weight * delta[2] * delta[2];
            if (score < min_score)
            {
                min_score = score;
                nearest_index = i;
            }
        }
        return nearest_index;
    };

    auto build_probe = [&](const std::string &name, BaseParticles &particles, Vecd *positions, const Vec3d &target_position) -> ProbeDescriptor
    {
        ProbeDescriptor probe;
        probe.name = name;
        probe.target_position = target_position;
        probe.particle_index = find_nearest_particle_index(particles, positions, target_position);
        probe.sampled_position = positions[probe.particle_index];
        return probe;
    };

    // NOTE:
    // Default probe coordinates are generated from current geometry bounds.
    // They can be overridden with environment variables:
    // TEAM7_PROBE_B_AIR_ABOVE_X/Y/Z, TEAM7_PROBE_B_AIR_BELOW_X/Y/Z,
    // TEAM7_PROBE_B_PLATE_CENTER_X/Y/Z, TEAM7_PROBE_T_PLATE_CENTER_X/Y/Z,
    // TEAM7_PROBE_T_PLATE_EDGE_XPLUS_X/Y/Z.
    Vec3d target_b_air_above_plate_default =
        plate_center + Vec3d(0.0, 0.0, 0.5 * plate_size[2] + 2.0 * dp_plate);
    Vec3d target_b_air_below_plate_default =
        plate_center - Vec3d(0.0, 0.0, 0.5 * plate_size[2] + 2.0 * dp_plate);
    Vec3d target_b_plate_center_default = plate_center;
    Vec3d target_t_plate_center_default = plate_center;
    Vec3d target_t_plate_edge_xplus_default =
        plate_center + Vec3d(0.4 * plate_size[0], 0.0, 0.0);

    Vec3d target_b_air_above_plate =
        get_env_vec3("TEAM7_PROBE_B_AIR_ABOVE", target_b_air_above_plate_default);
    Vec3d target_b_air_below_plate =
        get_env_vec3("TEAM7_PROBE_B_AIR_BELOW", target_b_air_below_plate_default);
    Vec3d target_b_plate_center =
        get_env_vec3("TEAM7_PROBE_B_PLATE_CENTER", target_b_plate_center_default);
    Vec3d target_t_plate_center =
        get_env_vec3("TEAM7_PROBE_T_PLATE_CENTER", target_t_plate_center_default);
    Vec3d target_t_plate_edge_xplus =
        get_env_vec3("TEAM7_PROBE_T_PLATE_EDGE_XPLUS", target_t_plate_edge_xplus_default);

    ProbeDescriptor probe_b_air_above_plate =
        build_probe("B_air_above_plate", air_particles, air_positions,
                    target_b_air_above_plate);
    ProbeDescriptor probe_b_air_below_plate =
        build_probe("B_air_below_plate", air_particles, air_positions,
                    target_b_air_below_plate);
    ProbeDescriptor probe_b_plate_center =
        build_probe("B_plate_center", plate_particles, plate_positions, target_b_plate_center);
    ProbeDescriptor probe_t_plate_center =
        build_probe("T_plate_center", plate_particles, plate_positions, target_t_plate_center);
    ProbeDescriptor probe_t_plate_edge_xplus =
        build_probe("T_plate_edge_xplus", plate_particles, plate_positions,
                    target_t_plate_edge_xplus);

    std::ofstream probe_locations_file(
        io_environment.OutputFolder() + "/team7_probe_locations.csv",
        std::ios::out | std::ios::trunc);
    probe_locations_file << std::setprecision(12);
    probe_locations_file << "name,body,particle_index,target_x,target_y,target_z,sampled_x,sampled_y,sampled_z\n";
    auto write_probe_location = [&](const ProbeDescriptor &probe, const std::string &body_name)
    {
        probe_locations_file << probe.name << ","
                             << body_name << ","
                             << probe.particle_index << ","
                             << probe.target_position[0] << ","
                             << probe.target_position[1] << ","
                             << probe.target_position[2] << ","
                             << probe.sampled_position[0] << ","
                             << probe.sampled_position[1] << ","
                             << probe.sampled_position[2] << "\n";
    };
    write_probe_location(probe_b_air_above_plate, "air");
    write_probe_location(probe_b_air_below_plate, "air");
    write_probe_location(probe_b_plate_center, "plate");
    write_probe_location(probe_t_plate_center, "plate");
    write_probe_location(probe_t_plate_edge_xplus, "plate");
    probe_locations_file.flush();

    std::vector<Team7CurveReferencePoint> b_curve_reference_points;
    std::vector<Team7CurveReferencePoint> j_curve_reference_points;
    bool b_curve_reference_loaded = load_team7_curve_reference_table(team7_table1_path, b_curve_reference_points);
    bool j_curve_reference_loaded = load_team7_curve_reference_table(team7_table2_path, j_curve_reference_points);

    if (enable_curve_validation)
    {
        std::cout << "[team7-curve-validation] table1_loaded=" << b_curve_reference_loaded
                  << ", table1_points=" << b_curve_reference_points.size()
                  << ", table2_loaded=" << j_curve_reference_loaded
                  << ", table2_points=" << j_curve_reference_points.size()
                  << std::endl;
        if (!b_curve_reference_loaded)
        {
            std::cout << "[team7-curve-validation] warning: cannot read Bz table at "
                      << team7_table1_path << std::endl;
        }
        if (!j_curve_reference_loaded)
        {
            std::cout << "[team7-curve-validation] warning: cannot read Jy table at "
                      << team7_table2_path << std::endl;
        }
    }

    std::vector<Team7CurveSamplePoint> b_curve_sample_points;
    std::vector<Team7CurveSamplePoint> j_curve_sample_points;
    if (enable_curve_validation && b_curve_reference_loaded)
    {
        b_curve_sample_points.reserve(b_curve_reference_points.size());
        for (const Team7CurveReferencePoint &reference_point : b_curve_reference_points)
        {
            Team7CurveSamplePoint sample_point;
            sample_point.quantity_name = "Bz";
            sample_point.x_mm = reference_point.x_mm;
            sample_point.reference_50 = reference_point.ref_50;
            sample_point.reference_200 = reference_point.ref_200;
            sample_point.reference_selected =
                select_reference_value_by_frequency(reference_point, curve_reference_frequency_hz);
            sample_point.target_position = Vec3d(reference_point.x_mm, b_curve_line_y, b_curve_line_z);
            sample_point.sample_index =
                find_nearest_particle_index(air_particles, air_positions, sample_point.target_position);
            sample_point.sampled_position = air_positions[sample_point.sample_index];
            b_curve_sample_points.push_back(sample_point);
        }
    }
    if (enable_curve_validation && j_curve_reference_loaded)
    {
        j_curve_sample_points.reserve(j_curve_reference_points.size());
        for (const Team7CurveReferencePoint &reference_point : j_curve_reference_points)
        {
            Team7CurveSamplePoint sample_point;
            sample_point.quantity_name = "Jy";
            sample_point.x_mm = reference_point.x_mm;
            sample_point.reference_50 = reference_point.ref_50;
            sample_point.reference_200 = reference_point.ref_200;
            sample_point.reference_selected =
                select_reference_value_by_frequency(reference_point, curve_reference_frequency_hz);
            sample_point.target_position = Vec3d(reference_point.x_mm, j_curve_line_y, j_curve_line_z);
            sample_point.sample_index = find_nearest_plate_particle_for_j_curve(sample_point.target_position);
            sample_point.sampled_position = plate_positions[sample_point.sample_index];
            j_curve_sample_points.push_back(sample_point);
        }
    }
    if (enable_curve_validation && (!b_curve_sample_points.empty() || !j_curve_sample_points.empty()))
    {
        std::ofstream curve_probe_locations_file(
            io_environment.OutputFolder() + "/team7_curve_probe_locations.csv",
            std::ios::out | std::ios::trunc);
        curve_probe_locations_file << std::setprecision(12);
        curve_probe_locations_file
            << "quantity,point_id,x_mm,target_x,target_y,target_z,sampled_x,sampled_y,sampled_z,particle_index\n";
        for (size_t i = 0; i != b_curve_sample_points.size(); ++i)
        {
            const Team7CurveSamplePoint &sample = b_curve_sample_points[i];
            curve_probe_locations_file
                << "Bz" << ","
                << i << ","
                << sample.x_mm << ","
                << sample.target_position[0] << ","
                << sample.target_position[1] << ","
                << sample.target_position[2] << ","
                << sample.sampled_position[0] << ","
                << sample.sampled_position[1] << ","
                << sample.sampled_position[2] << ","
                << sample.sample_index << "\n";
        }
        for (size_t i = 0; i != j_curve_sample_points.size(); ++i)
        {
            const Team7CurveSamplePoint &sample = j_curve_sample_points[i];
            curve_probe_locations_file
                << "Jy" << ","
                << i << ","
                << sample.x_mm << ","
                << sample.target_position[0] << ","
                << sample.target_position[1] << ","
                << sample.target_position[2] << ","
                << sample.sampled_position[0] << ","
                << sample.sampled_position[1] << ","
                << sample.sampled_position[2] << ","
                << sample.sample_index << "\n";
        }
        curve_probe_locations_file.flush();
    }

    auto interpolate_complex_curve_component =
        [&](Vecd *positions,
            Real *volumes,
            size_t total_particles,
            Vecd *field_real,
            Vecd *field_imag,
            const Team7CurveSamplePoint &sample,
            int component_index,
            Real interpolation_radius,
            ComplexCurveComponentMode complex_mode) -> Real
    {
        Real safe_radius = SMAX(interpolation_radius, TinyReal);
        Real support_radius_sq = safe_radius * safe_radius;
        Real epsilon_sq = 1.0e-4 * support_radius_sq + TinyReal;
        Real weight_sum = 0.0;
        Real interpolated_real = 0.0;
        Real interpolated_imag = 0.0;

        for (size_t i = 0; i != total_particles; ++i)
        {
            Real distance_sq = (positions[i] - sample.target_position).squaredNorm();
            if (!(distance_sq <= support_radius_sq))
            {
                continue;
            }
            Real particle_volume = volumes != nullptr ? volumes[i] : 1.0;
            Real weight = particle_volume / (distance_sq + epsilon_sq);
            interpolated_real += weight * field_real[i][component_index];
            if (field_imag != nullptr)
            {
                interpolated_imag += weight * field_imag[i][component_index];
            }
            weight_sum += weight;
        }

        if (weight_sum > TinyReal && std::isfinite(weight_sum))
        {
            interpolated_real /= weight_sum;
            interpolated_imag /= weight_sum;
            return project_complex_component(interpolated_real, interpolated_imag, complex_mode);
        }

        Real fallback_real = field_real[sample.sample_index][component_index];
        Real fallback_imag = field_imag != nullptr
                                 ? field_imag[sample.sample_index][component_index]
                                 : 0.0;
        return project_complex_component(fallback_real, fallback_imag, complex_mode);
    };
    auto analyze_curve_support =
        [&](Vecd *positions,
            Real *volumes,
            size_t total_particles,
            const Team7CurveSamplePoint &sample,
            Real interpolation_radius) -> CurveInterpolationSupportInfo
    {
        CurveInterpolationSupportInfo info;
        if (positions == nullptr || total_particles == 0)
        {
            return info;
        }

        Real safe_radius = SMAX(interpolation_radius, TinyReal);
        Real support_radius_sq = safe_radius * safe_radius;
        Real epsilon_sq = 1.0e-4 * support_radius_sq + TinyReal;
        Real weighted_distance_sum = 0.0;
        Real nearest_distance_sq = std::numeric_limits<Real>::max();

        for (size_t i = 0; i != total_particles; ++i)
        {
            Real distance_sq = (positions[i] - sample.target_position).squaredNorm();
            if (!(distance_sq <= support_radius_sq))
            {
                continue;
            }
            Real particle_volume = volumes != nullptr ? volumes[i] : 1.0;
            Real weight = particle_volume / (distance_sq + epsilon_sq);
            info.support_point_count++;
            info.weight_sum += weight;
            weighted_distance_sum += weight * sqrt(distance_sq);
            nearest_distance_sq = SMIN(nearest_distance_sq, distance_sq);
        }

        if (info.support_point_count > 0)
        {
            info.nearest_distance = sqrt(nearest_distance_sq);
        }
        if (info.weight_sum > TinyReal && std::isfinite(info.weight_sum))
        {
            info.weighted_mean_distance = weighted_distance_sum / info.weight_sum;
        }
        return info;
    };
    auto compute_box_boundary_clearance =
        [&](const Vec3d &position, const BoundingBoxd &bounds) -> Vec3d
    {
        return Vec3d(
            SMIN(position[0] - bounds.lower_[0], bounds.upper_[0] - position[0]),
            SMIN(position[1] - bounds.lower_[1], bounds.upper_[1] - position[1]),
            SMIN(position[2] - bounds.lower_[2], bounds.upper_[2] - position[2]));
    };
    auto assign_curve_sample_quality =
        [&](Team7CurveSamplePoint &sample,
            bool is_bz_curve,
            Real nominal_spacing,
            Real interpolation_radius)
    {
        Vecd *positions = is_bz_curve ? air_positions : plate_positions;
        Real *volumes = is_bz_curve ? air_vol : plate_vol;
        size_t total_particles = is_bz_curve ? air_particles.TotalRealParticles()
                                             : plate_particles.TotalRealParticles();
        const BoundingBoxd &bounds = is_bz_curve ? air_bounds : plate_bounds;
        CurveInterpolationSupportInfo support_info =
            analyze_curve_support(positions, volumes, total_particles, sample, interpolation_radius);
        Vec3d clearance = compute_box_boundary_clearance(sample.sampled_position, bounds);
        sample.support_point_count = support_info.support_point_count;
        sample.support_weight_sum = support_info.weight_sum;
        sample.support_nearest_distance = support_info.nearest_distance;
        sample.support_weighted_mean_distance = support_info.weighted_mean_distance;
        sample.boundary_clearance_x = clearance[0];
        sample.boundary_clearance_y = clearance[1];
        sample.boundary_clearance_z = clearance[2];

        bool enough_support = sample.support_point_count >= 3;
        bool near_target = std::isfinite(sample.support_nearest_distance) &&
                           sample.support_nearest_distance <= nominal_spacing + TinyReal;
        bool away_from_curve_end = sample.boundary_clearance_x >= nominal_spacing - TinyReal;
        sample.high_quality = enough_support && near_target && away_from_curve_end;
    };
    auto interpolate_scalar_curve_value =
        [&](Vecd *positions,
            Real *volumes,
            size_t total_particles,
            Real *scalar_field,
            const Team7CurveSamplePoint &sample,
            Real interpolation_radius) -> Real
    {
        if (scalar_field == nullptr || total_particles == 0)
        {
            return 0.0;
        }
        Real safe_radius = SMAX(interpolation_radius, TinyReal);
        Real support_radius_sq = safe_radius * safe_radius;
        Real epsilon_sq = 1.0e-4 * support_radius_sq + TinyReal;
        Real weighted_scalar = 0.0;
        Real weight_sum = 0.0;

        for (size_t i = 0; i != total_particles; ++i)
        {
            Real distance_sq = (positions[i] - sample.target_position).squaredNorm();
            if (!(distance_sq <= support_radius_sq))
            {
                continue;
            }
            Real particle_volume = volumes != nullptr ? volumes[i] : 1.0;
            Real weight = particle_volume / (distance_sq + epsilon_sq);
            weighted_scalar += weight * scalar_field[i];
            weight_sum += weight;
        }
        if (weight_sum > TinyReal && std::isfinite(weight_sum))
        {
            return weighted_scalar / weight_sum;
        }
        return scalar_field[sample.sample_index];
    };

    auto evaluate_j_curve_reference_to_simulated_abs_ratio = [&]() -> Real
    {
        if (!use_frequency_aphi || !enable_curve_validation || !auto_normalize_use_j_curve ||
            j_curve_sample_points.empty() || plate_current_density_real == nullptr)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }

        Real sum_abs_reference = 0.0;
        Real sum_abs_simulated = 0.0;
        size_t valid_points = 0;
        for (size_t i = 0; i != j_curve_sample_points.size(); ++i)
        {
            const Team7CurveSamplePoint &sample = j_curve_sample_points[i];
            if (!std::isfinite(sample.reference_selected))
            {
                continue;
            }
            Real simulated_value =
                interpolate_complex_curve_component(plate_positions,
                                                    plate_vol,
                                                    plate_particles.TotalRealParticles(),
                                                    plate_current_density_real,
                                                    plate_current_density_imag,
                                                    sample,
                                                    j_curve_component_index,
                                                    j_curve_interpolation_radius,
                                                    j_curve_complex_mode) *
                j_curve_unit_scale;
            sum_abs_reference += fabs(sample.reference_selected);
            sum_abs_simulated += fabs(simulated_value);
            valid_points++;
        }
        if (valid_points == 0 || sum_abs_simulated <= TinyReal)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        return (sum_abs_reference + TinyReal) / (sum_abs_simulated + TinyReal);
    };
    auto try_auto_normalize_runtime_source_scale = [&](size_t current_outer_iteration)
    {
        if (!auto_normalize_source || current_outer_iteration % auto_normalize_interval != 0)
        {
            return;
        }
        Real reference_to_simulated_ratio = evaluate_j_curve_reference_to_simulated_abs_ratio();
        if (!std::isfinite(reference_to_simulated_ratio))
        {
            return;
        }
        Real bounded_ratio = SMIN(auto_normalize_max_factor,
                                  SMAX(auto_normalize_min_factor, reference_to_simulated_ratio));
        Real old_scale = runtime_source_scale;
        Real target_scale = old_scale * bounded_ratio;
        runtime_source_scale =
            old_scale + auto_normalize_gain * (target_scale - old_scale);
        runtime_source_scale =
            SMIN(runtime_source_scale_max, SMAX(runtime_source_scale_min, runtime_source_scale));
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-source-normalize] iter=" << current_outer_iteration
                  << ", ratio_ref_to_sim=" << reference_to_simulated_ratio
                  << ", bounded_ratio=" << bounded_ratio
                  << ", source_scale_old=" << old_scale
                  << ", source_scale_new=" << runtime_source_scale
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    };

    auto sample_magnetic_flux_density = [&](AngularVecd *vector_potential_curl, const ProbeDescriptor &probe) -> Vec3d
    {
        const AngularVecd &b = vector_potential_curl[probe.particle_index];
        return Vec3d(b[0], b[1], b[2]);
    };
    auto sample_magnetic_flux_density_magnitude =
        [&](AngularVecd *vector_potential_curl_real_part,
            AngularVecd *vector_potential_curl_imag_part,
            const ProbeDescriptor &probe) -> Real
    {
        Vec3d b_real = sample_magnetic_flux_density(vector_potential_curl_real_part, probe);
        Real magnetic_flux_density_square = b_real.squaredNorm();
        if (use_frequency_aphi && vector_potential_curl_imag_part != nullptr)
        {
            Vec3d b_imag = sample_magnetic_flux_density(vector_potential_curl_imag_part, probe);
            magnetic_flux_density_square += b_imag.squaredNorm();
        }
        return sqrt(magnetic_flux_density_square);
    };
    auto compute_biot_savart_magnetic_flux_density =
        [&](const Vec3d &sample_position) -> std::pair<Vec3d, Vec3d>
    {
        const Vecd *source_current_density_real_part =
            (use_frequency_aphi && coil_source_current_density_real != nullptr)
                ? coil_source_current_density_real
                : coil_source_current_density;
        const Vecd *source_current_density_imag_part =
            (use_frequency_aphi && coil_source_current_density_imag != nullptr)
                ? coil_source_current_density_imag
                : nullptr;

        Vec3d b_real = Vec3d::Zero();
        Vec3d b_imag = Vec3d::Zero();
        if (source_current_density_real_part == nullptr || coil_particle_positions == nullptr || coil_vol == nullptr)
        {
            return {b_real, b_imag};
        }

        const Real mu0_over_4pi = static_cast<Real>(1.0e-7);
        const Real geom_to_si = geom_length_to_m;
        const Real geom_vol_to_si = geom_volume_to_m3;
        const size_t total_coil_particles = coil_particles.TotalRealParticles();

        for (size_t i = 0; i != total_coil_particles; ++i)
        {
            Vec3d r_si = (sample_position - coil_particle_positions[i]) * geom_to_si;
            Real r_sq = r_si.squaredNorm();
            if (r_sq <= TinyReal * TinyReal)
            {
                continue;
            }
            Real inv_r3 = static_cast<Real>(1.0) / (sqrt(r_sq) * r_sq + TinyReal);
            Real dvol_si = coil_vol[i] * geom_vol_to_si;
            Real kernel_factor = dvol_si * inv_r3;

            Vec3d j_real = source_current_density_real_part[i];
            b_real += j_real.cross(r_si) * kernel_factor;
            if (source_current_density_imag_part != nullptr)
            {
                Vec3d j_imag = source_current_density_imag_part[i];
                b_imag += j_imag.cross(r_si) * kernel_factor;
            }
        }

        b_real *= mu0_over_4pi;
        b_imag *= mu0_over_4pi;
        return {b_real, b_imag};
    };
    auto complex_vector_magnitude = [&](const Vec3d &real_part, const Vec3d &imag_part) -> Real
    {
        return sqrt(real_part.squaredNorm() + imag_part.squaredNorm());
    };
    auto evaluate_coil_air_probe_reference_ratio = [&]() -> Real
    {
        Real b_air_above_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_above_plate);
        Real b_air_below_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_below_plate);
        std::pair<Vec3d, Vec3d> b_ref_air_above =
            compute_biot_savart_magnetic_flux_density(probe_b_air_above_plate.sampled_position);
        std::pair<Vec3d, Vec3d> b_ref_air_below =
            compute_biot_savart_magnetic_flux_density(probe_b_air_below_plate.sampled_position);
        Real b_ref_air_above_mag =
            complex_vector_magnitude(b_ref_air_above.first, b_ref_air_above.second);
        Real b_ref_air_below_mag =
            complex_vector_magnitude(b_ref_air_below.first, b_ref_air_below.second);
        auto probe_ratio = [&](Real simulated_mag, Real reference_mag) -> Real
        {
            if (!std::isfinite(simulated_mag) || !std::isfinite(reference_mag) ||
                reference_mag <= TinyReal)
            {
                return std::numeric_limits<Real>::quiet_NaN();
            }
            return simulated_mag / (reference_mag + TinyReal);
        };
        Real above_ratio = probe_ratio(b_air_above_mag, b_ref_air_above_mag);
        Real below_ratio = probe_ratio(b_air_below_mag, b_ref_air_below_mag);
        if (!std::isfinite(above_ratio) && !std::isfinite(below_ratio))
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        if (!std::isfinite(above_ratio))
        {
            return below_ratio;
        }
        if (!std::isfinite(below_ratio))
        {
            return above_ratio;
        }
        return SMAX(above_ratio, below_ratio);
    };
    auto try_update_runtime_source_scale_by_coil_air_probe =
        [&](size_t current_outer_iteration)
    {
        if (!enable_coil_air_probe_feedback || !use_frequency_aphi ||
            !apply_runtime_source_scale ||
            (enable_vacuum_reluctivity_feedback &&
             apply_vacuum_reluctivity_scale_runtime) ||
            current_outer_iteration % coil_air_probe_feedback_interval != 0)
        {
            return;
        }
        Real simulated_to_reference_ratio = evaluate_coil_air_probe_reference_ratio();
        if (!std::isfinite(simulated_to_reference_ratio) ||
            simulated_to_reference_ratio <= TinyReal)
        {
            return;
        }
        Real reference_to_simulated_ratio =
            coil_air_probe_feedback_target_ratio /
            (simulated_to_reference_ratio + TinyReal);
        Real bounded_ratio =
            SMIN(coil_air_probe_feedback_max_factor,
                 SMAX(coil_air_probe_feedback_min_factor,
                      reference_to_simulated_ratio));
        Real old_scale = runtime_source_scale;
        Real target_scale = old_scale * bounded_ratio;
        runtime_source_scale =
            old_scale + coil_air_probe_feedback_gain * (target_scale - old_scale);
        runtime_source_scale =
            SMIN(runtime_source_scale_max,
                 SMAX(runtime_source_scale_min, runtime_source_scale));
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-coil-air-probe-feedback] iter="
                  << current_outer_iteration
                  << ", sim_over_ref=" << simulated_to_reference_ratio
                  << ", target_sim_over_ref=" << coil_air_probe_feedback_target_ratio
                  << ", ratio_ref_to_sim=" << reference_to_simulated_ratio
                  << ", bounded_ratio=" << bounded_ratio
                  << ", source_scale_old=" << old_scale
                  << ", source_scale_new=" << runtime_source_scale
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    };
    auto try_update_vacuum_reluctivity_scale_by_coil_air_probe =
        [&](size_t current_outer_iteration)
    {
        if (!enable_vacuum_reluctivity_feedback || !use_frequency_aphi ||
            !apply_vacuum_reluctivity_scale_runtime ||
            current_outer_iteration % vacuum_reluctivity_feedback_interval != 0)
        {
            return;
        }
        Real simulated_to_reference_ratio = evaluate_coil_air_probe_reference_ratio();
        if (!std::isfinite(simulated_to_reference_ratio) ||
            simulated_to_reference_ratio <= TinyReal)
        {
            return;
        }
        Real sim_to_target_ratio =
            simulated_to_reference_ratio /
            (vacuum_reluctivity_feedback_target_ratio + TinyReal);
        Real bounded_factor =
            SMIN(vacuum_reluctivity_feedback_max_factor,
                 SMAX(vacuum_reluctivity_feedback_min_factor, sim_to_target_ratio));
        Real old_scale = vacuum_reluctivity_scale_runtime;
        Real target_scale = old_scale * bounded_factor;
        vacuum_reluctivity_scale_runtime =
            old_scale + vacuum_reluctivity_feedback_gain * (target_scale - old_scale);
        vacuum_reluctivity_scale_runtime =
            SMIN(vacuum_reluctivity_scale_max,
                 SMAX(vacuum_reluctivity_scale_min, vacuum_reluctivity_scale_runtime));
        apply_vacuum_reluctivity_scale(vacuum_reluctivity_scale_runtime);
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-vacuum-reluctivity-feedback] iter="
                  << current_outer_iteration
                  << ", sim_over_ref=" << simulated_to_reference_ratio
                  << ", target_sim_over_ref=" << vacuum_reluctivity_feedback_target_ratio
                  << ", sim_to_target=" << sim_to_target_ratio
                  << ", bounded_factor=" << bounded_factor
                  << ", reluctivity_scale_old=" << old_scale
                  << ", reluctivity_scale_new=" << vacuum_reluctivity_scale_runtime
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    };
    auto interpolate_complex_air_curl_at_position =
        [&](const Vec3d &target_position, Real interpolation_radius) -> std::pair<Vec3d, Vec3d>
    {
        Vec3d interpolated_real = Vec3d::Zero();
        Vec3d interpolated_imag = Vec3d::Zero();
        if (air_positions == nullptr || air_vol == nullptr || air_vector_potential_curl_probe == nullptr)
        {
            return {interpolated_real, interpolated_imag};
        }

        Real safe_radius = SMAX(interpolation_radius, TinyReal);
        Real support_radius_sq = safe_radius * safe_radius;
        Real epsilon_sq = 1.0e-4 * support_radius_sq + TinyReal;
        Real weight_sum = 0.0;
        size_t total_air_particles = air_particles.TotalRealParticles();

        for (size_t i = 0; i != total_air_particles; ++i)
        {
            Real distance_sq = (air_positions[i] - target_position).squaredNorm();
            if (!(distance_sq <= support_radius_sq))
            {
                continue;
            }
            Real weight = air_vol[i] / (distance_sq + epsilon_sq);
            const AngularVecd &b_real_angular = air_vector_potential_curl_probe[i];
            interpolated_real += weight * Vec3d(b_real_angular[0], b_real_angular[1], b_real_angular[2]);
            if (use_frequency_aphi && air_vector_potential_curl_imag != nullptr)
            {
                const AngularVecd &b_imag_angular = air_vector_potential_curl_imag[i];
                interpolated_imag += weight * Vec3d(b_imag_angular[0], b_imag_angular[1], b_imag_angular[2]);
            }
            weight_sum += weight;
        }

        if (weight_sum > TinyReal && std::isfinite(weight_sum))
        {
            interpolated_real /= weight_sum;
            interpolated_imag /= weight_sum;
            return {interpolated_real, interpolated_imag};
        }

        size_t fallback_index = find_nearest_particle_index(air_particles, air_positions, target_position);
        const AngularVecd &fallback_real_angular = air_vector_potential_curl_probe[fallback_index];
        Vec3d fallback_real = Vec3d(fallback_real_angular[0], fallback_real_angular[1], fallback_real_angular[2]);
        Vec3d fallback_imag = Vec3d::Zero();
        if (use_frequency_aphi && air_vector_potential_curl_imag != nullptr)
        {
            const AngularVecd &fallback_imag_angular = air_vector_potential_curl_imag[fallback_index];
            fallback_imag = Vec3d(fallback_imag_angular[0], fallback_imag_angular[1], fallback_imag_angular[2]);
        }
        return {fallback_real, fallback_imag};
    };

    auto evaluate_plate_diagnostics = [&]() -> PlateDiagnostics
    {
        PlateDiagnostics diagnostics;
        diagnostics.max_joule_heat = -std::numeric_limits<Real>::max();
        diagnostics.max_current_density_magnitude = -std::numeric_limits<Real>::max();
        diagnostics.max_temperature = -std::numeric_limits<Real>::max();

        Real temperature_volume_moment = 0.0;
        Real current_density_magnitude_volume_moment = 0.0;
        size_t total_real_particles = plate_particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Real vol = plate_vol[i] * geom_volume_to_m3;
            Real joule_heat = plate_joule_heat[i];
            Real current_density_magnitude = 0.0;
            Real a_rate_magnitude = 0.0;
            Real a_dot_magnitude = 0.0;
            Real grad_phi_magnitude = 0.0;
            Real electric_field_magnitude = 0.0;
            if (use_frequency_aphi)
            {
                current_density_magnitude =
                    sqrt(plate_current_density_real[i].squaredNorm() + plate_current_density_imag[i].squaredNorm());
                a_rate_magnitude =
                    sqrt(plate_vector_potential_change_rate_real[i].squaredNorm() +
                         plate_vector_potential_change_rate_imag[i].squaredNorm());
                a_dot_magnitude = a_rate_magnitude;
                grad_phi_magnitude =
                    sqrt(plate_electric_potential_gradient_real[i].squaredNorm() +
                         plate_electric_potential_gradient_imag[i].squaredNorm());
                electric_field_magnitude =
                    sqrt(plate_electric_field_real[i].squaredNorm() + plate_electric_field_imag[i].squaredNorm());
            }
            else
            {
                current_density_magnitude = plate_current_density[i].norm();
                a_rate_magnitude = plate_vector_potential_change_rate[i].norm();
                a_dot_magnitude = plate_vector_potential_dt[i].norm();
                grad_phi_magnitude = plate_electric_potential_gradient[i].norm();
                electric_field_magnitude = plate_electric_field[i].norm();
            }
            Real temperature = plate_temperature[i];

            diagnostics.total_volume += vol;
            diagnostics.total_joule_power += joule_heat * vol;
            temperature_volume_moment += temperature * vol;
            current_density_magnitude_volume_moment += current_density_magnitude * vol;
            diagnostics.max_joule_heat = SMAX(diagnostics.max_joule_heat, joule_heat);
            diagnostics.max_current_density_magnitude =
                SMAX(diagnostics.max_current_density_magnitude, current_density_magnitude);
            diagnostics.max_temperature = SMAX(diagnostics.max_temperature, temperature);
            diagnostics.max_a_rate = SMAX(diagnostics.max_a_rate, a_rate_magnitude);
            diagnostics.max_a_dot = SMAX(diagnostics.max_a_dot, a_dot_magnitude);
            diagnostics.max_grad_phi = SMAX(diagnostics.max_grad_phi, grad_phi_magnitude);
            diagnostics.max_electric_field = SMAX(diagnostics.max_electric_field, electric_field_magnitude);
        }

        Real total_volume = diagnostics.total_volume + TinyReal;
        diagnostics.avg_joule_heat = diagnostics.total_joule_power / total_volume;
        diagnostics.avg_current_density_magnitude =
            current_density_magnitude_volume_moment / total_volume;
        diagnostics.avg_temperature = temperature_volume_moment / total_volume;
        return diagnostics;
    };
    Real em_runtime_pseudo_dt_scale = 1.0;
    auto evaluate_max_em_change_rate = [&]() -> Real
    {
        Real max_change_rate = 0.0;
        size_t total_plate_particles = plate_particles.TotalRealParticles();
        size_t total_coil_particles = coil_particles.TotalRealParticles();
        size_t total_air_particles = air_particles.TotalRealParticles();
        auto accumulate_frequency_rate =
            [&](Vecd *change_rate_real,
                Vecd *change_rate_imag,
                size_t total_particles,
                Real body_relax_scaling)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real rate_norm =
                    sqrt(change_rate_real[i].squaredNorm() +
                         change_rate_imag[i].squaredNorm());
                if (em_rate_use_effective_update)
                {
                    rate_norm *= body_relax_scaling;
                    if (em_rate_use_pseudo_dt_scaling)
                    {
                        rate_norm *= em_runtime_pseudo_dt_scale;
                    }
                }
                max_change_rate = SMAX(max_change_rate, rate_norm);
            }
        };
        auto accumulate_harmonic_rate =
            [&](Vecd *change_rate,
                size_t total_particles,
                Real body_relax_scaling)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real rate_norm = change_rate[i].norm();
                if (em_rate_use_effective_update)
                {
                    rate_norm *= body_relax_scaling;
                    if (em_rate_use_pseudo_dt_scaling)
                    {
                        rate_norm *= em_runtime_pseudo_dt_scale;
                    }
                }
                max_change_rate = SMAX(max_change_rate, rate_norm);
            }
        };
        if (use_frequency_aphi)
        {
            if (!use_frequency_operator_coil_air_only)
            {
                accumulate_frequency_rate(plate_vector_potential_change_rate_real,
                                          plate_vector_potential_change_rate_imag,
                                          total_plate_particles,
                                          plate_a_relaxation_scaling);
            }
            accumulate_frequency_rate(coil_vector_potential_change_rate_real,
                                      coil_vector_potential_change_rate_imag,
                                      total_coil_particles,
                                      coil_a_relaxation_scaling);
            if (em_rate_include_air)
            {
                accumulate_frequency_rate(air_vector_potential_change_rate_real,
                                          air_vector_potential_change_rate_imag,
                                          total_air_particles,
                                          air_a_relaxation_scaling);
            }
        }
        else
        {
            if (!use_frequency_operator_coil_air_only)
            {
                accumulate_harmonic_rate(plate_vector_potential_change_rate,
                                         total_plate_particles,
                                         plate_a_relaxation_scaling);
            }
            accumulate_harmonic_rate(coil_vector_potential_change_rate,
                                     total_coil_particles,
                                     coil_a_relaxation_scaling);
            if (em_rate_include_air)
            {
                accumulate_harmonic_rate(air_vector_potential_change_rate,
                                         total_air_particles,
                                         air_a_relaxation_scaling);
            }
        }
        return max_change_rate;
    };
    auto compute_frequency_conductor_residual_terms =
        [&](Vecd vector_potential_real_i,
            Vecd vector_potential_imag_i,
            Vecd electric_potential_gradient_real_i,
            Vecd electric_potential_gradient_imag_i,
            Vecd curl_nu_b_real_i,
            Vecd curl_nu_b_imag_i,
            Vecd source_real_i,
            Vecd source_imag_i,
            Real sigma_i,
            Real &residual_norm,
            Real &reference_norm)
    {
        Vecd rhs_real = source_real_i - curl_nu_b_real_i -
                        sigma_i * electric_potential_gradient_real_i +
                        sigma_i * harmonic_angular_frequency_runtime * vector_potential_imag_i;
        Vecd rhs_imag = source_imag_i - curl_nu_b_imag_i -
                        sigma_i * electric_potential_gradient_imag_i -
                        sigma_i * harmonic_angular_frequency_runtime * vector_potential_real_i;
        residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());
        reference_norm =
            sqrt(source_real_i.squaredNorm() + source_imag_i.squaredNorm()) +
            sqrt(curl_nu_b_real_i.squaredNorm() + curl_nu_b_imag_i.squaredNorm()) +
            sigma_i * sqrt(electric_potential_gradient_real_i.squaredNorm() +
                           electric_potential_gradient_imag_i.squaredNorm()) +
            sigma_i * harmonic_angular_frequency_runtime *
                sqrt(vector_potential_real_i.squaredNorm() +
                     vector_potential_imag_i.squaredNorm());
    };
    auto compute_frequency_magnetic_residual_terms =
        [&](Vecd curl_nu_b_real_i,
            Vecd curl_nu_b_imag_i,
            Vecd source_real_i,
            Vecd source_imag_i,
            Real &residual_norm,
            Real &reference_norm)
    {
        Vecd rhs_real = source_real_i - curl_nu_b_real_i;
        Vecd rhs_imag = source_imag_i - curl_nu_b_imag_i;
        residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());
        reference_norm =
            sqrt(source_real_i.squaredNorm() + source_imag_i.squaredNorm()) +
            sqrt(curl_nu_b_real_i.squaredNorm() + curl_nu_b_imag_i.squaredNorm());
    };
    auto evaluate_em_residual_reference_scale = [&]() -> Real
    {
        if (!use_frequency_aphi)
        {
            return 1.0;
        }

        Real max_reference_scale = 0.0;
        auto update_conductor_reference_scale =
            [&](Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_conductor_residual_terms(
                    vector_potential_real[i], vector_potential_imag[i],
                    electric_potential_gradient_real[i], electric_potential_gradient_imag[i],
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    electrical_conductivity[i], residual_norm, reference_norm);
                max_reference_scale = SMAX(max_reference_scale, reference_norm);
            }
        };
        auto update_magnetic_reference_scale =
            [&](Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                size_t total_particles)
        {
            // Source-free magnetic-only regions such as the surrounding air should be
            // measured against the driven conductor/coil scale, not define it.
            if (source_current_density_real == nullptr &&
                source_current_density_imag == nullptr)
            {
                return;
            }
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_magnetic_residual_terms(
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    residual_norm, reference_norm);
                max_reference_scale = SMAX(max_reference_scale, reference_norm);
            }
        };

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        if (!use_frequency_operator_coil_air_only)
        {
            update_conductor_reference_scale(plate_vector_potential_real,
                                             plate_vector_potential_imag,
                                             plate_electric_potential_gradient_real,
                                             plate_electric_potential_gradient_imag,
                                             plate_curl_nu_b_real,
                                             plate_curl_nu_b_imag,
                                             nullptr,
                                             nullptr,
                                             plate_sigma,
                                             plate_particles.TotalRealParticles());
        }
        if (use_frequency_coil_magnetic_only)
        {
            update_magnetic_reference_scale(coil_curl_nu_b_real,
                                            coil_curl_nu_b_imag,
                                            coil_source_current_density_real,
                                            coil_source_current_density_imag,
                                            coil_particles.TotalRealParticles());
        }
        else
        {
            update_conductor_reference_scale(coil_vector_potential_real,
                                             coil_vector_potential_imag,
                                             coil_electric_potential_gradient_real,
                                             coil_electric_potential_gradient_imag,
                                             coil_curl_nu_b_real,
                                             coil_curl_nu_b_imag,
                                             coil_source_current_density_real,
                                             coil_source_current_density_imag,
                                             coil_sigma,
                                             coil_particles.TotalRealParticles());
        }
        if (em_rate_include_air)
        {
            update_magnetic_reference_scale(air_curl_nu_b_real,
                                            air_curl_nu_b_imag,
                                            nullptr,
                                            nullptr,
                                            air_particles.TotalRealParticles());
        }
        return SMAX(max_reference_scale, TinyReal);
    };
    auto evaluate_max_em_residual_ratio = [&]() -> Real
    {
        if (!use_frequency_aphi)
        {
            return 0.0;
        }

        Real max_residual_ratio = 0.0;
        Real global_reference_scale = evaluate_em_residual_reference_scale();
        auto update_conductor_residual_ratio =
            [&](Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_conductor_residual_terms(
                    vector_potential_real[i], vector_potential_imag[i],
                    electric_potential_gradient_real[i], electric_potential_gradient_imag[i],
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    electrical_conductivity[i], residual_norm, reference_norm);
                max_residual_ratio =
                    SMAX(max_residual_ratio, residual_norm / (global_reference_scale + TinyReal));
            }
        };
        auto update_magnetic_only_residual_ratio =
            [&](Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_magnetic_residual_terms(
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    residual_norm, reference_norm);
                max_residual_ratio =
                    SMAX(max_residual_ratio, residual_norm / (global_reference_scale + TinyReal));
            }
        };

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        size_t total_plate_particles = plate_particles.TotalRealParticles();
        size_t total_coil_particles = coil_particles.TotalRealParticles();
        size_t total_air_particles = air_particles.TotalRealParticles();

        if (!use_frequency_operator_coil_air_only)
        {
            update_conductor_residual_ratio(plate_vector_potential_real,
                                            plate_vector_potential_imag,
                                            plate_electric_potential_gradient_real,
                                            plate_electric_potential_gradient_imag,
                                            plate_curl_nu_b_real,
                                            plate_curl_nu_b_imag,
                                            nullptr,
                                            nullptr,
                                            plate_sigma,
                                            total_plate_particles);
        }
        if (use_frequency_coil_magnetic_only)
        {
            update_magnetic_only_residual_ratio(coil_curl_nu_b_real,
                                                coil_curl_nu_b_imag,
                                                coil_source_current_density_real,
                                                coil_source_current_density_imag,
                                                total_coil_particles);
        }
        else
        {
            update_conductor_residual_ratio(coil_vector_potential_real,
                                            coil_vector_potential_imag,
                                            coil_electric_potential_gradient_real,
                                            coil_electric_potential_gradient_imag,
                                            coil_curl_nu_b_real,
                                            coil_curl_nu_b_imag,
                                            coil_source_current_density_real,
                                            coil_source_current_density_imag,
                                            coil_sigma,
                                            total_coil_particles);
        }
        if (em_rate_include_air)
        {
            update_magnetic_only_residual_ratio(air_curl_nu_b_real,
                                                air_curl_nu_b_imag,
                                                nullptr,
                                                nullptr,
                                                total_air_particles);
        }
        return max_residual_ratio;
    };
    auto evaluate_em_rate_breakdown = [&]() -> EmRateBreakdown
    {
        EmRateBreakdown breakdown;
        size_t total_plate_particles = plate_particles.TotalRealParticles();
        size_t total_coil_particles = coil_particles.TotalRealParticles();
        size_t total_air_particles = air_particles.TotalRealParticles();
        auto update_frequency_body_rate =
            [&](Real &body_max_rate,
                Vecd *change_rate_real,
                Vecd *change_rate_imag,
                size_t total_particles,
                Real body_relax_scaling)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real rate_norm =
                    sqrt(change_rate_real[i].squaredNorm() +
                         change_rate_imag[i].squaredNorm());
                if (em_rate_use_effective_update)
                {
                    rate_norm *= body_relax_scaling;
                    if (em_rate_use_pseudo_dt_scaling)
                    {
                        rate_norm *= em_runtime_pseudo_dt_scale;
                    }
                }
                body_max_rate = SMAX(body_max_rate, rate_norm);
            }
        };
        auto update_harmonic_body_rate =
            [&](Real &body_max_rate,
                Vecd *change_rate,
                size_t total_particles,
                Real body_relax_scaling)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real rate_norm = change_rate[i].norm();
                if (em_rate_use_effective_update)
                {
                    rate_norm *= body_relax_scaling;
                    if (em_rate_use_pseudo_dt_scaling)
                    {
                        rate_norm *= em_runtime_pseudo_dt_scale;
                    }
                }
                body_max_rate = SMAX(body_max_rate, rate_norm);
            }
        };
        if (use_frequency_aphi)
        {
            if (!use_frequency_operator_coil_air_only)
            {
                update_frequency_body_rate(breakdown.plate,
                                           plate_vector_potential_change_rate_real,
                                           plate_vector_potential_change_rate_imag,
                                           total_plate_particles,
                                           plate_a_relaxation_scaling);
            }
            update_frequency_body_rate(breakdown.coil,
                                       coil_vector_potential_change_rate_real,
                                       coil_vector_potential_change_rate_imag,
                                       total_coil_particles,
                                       coil_a_relaxation_scaling);
            update_frequency_body_rate(breakdown.air,
                                       air_vector_potential_change_rate_real,
                                       air_vector_potential_change_rate_imag,
                                       total_air_particles,
                                       air_a_relaxation_scaling);
        }
        else
        {
            if (!use_frequency_operator_coil_air_only)
            {
                update_harmonic_body_rate(breakdown.plate,
                                          plate_vector_potential_change_rate,
                                          total_plate_particles,
                                          plate_a_relaxation_scaling);
            }
            update_harmonic_body_rate(breakdown.coil,
                                      coil_vector_potential_change_rate,
                                      total_coil_particles,
                                      coil_a_relaxation_scaling);
            update_harmonic_body_rate(breakdown.air,
                                      air_vector_potential_change_rate,
                                      total_air_particles,
                                      air_a_relaxation_scaling);
        }
        return breakdown;
    };
    auto evaluate_em_residual_breakdown = [&]() -> EmResidualBreakdown
    {
        EmResidualBreakdown breakdown;
        if (!use_frequency_aphi)
        {
            return breakdown;
        }

        Real global_reference_scale = evaluate_em_residual_reference_scale();
        auto update_conductor_body_residual =
            [&](Real &body_max_residual,
                Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_conductor_residual_terms(
                    vector_potential_real[i], vector_potential_imag[i],
                    electric_potential_gradient_real[i], electric_potential_gradient_imag[i],
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    electrical_conductivity[i], residual_norm, reference_norm);
                body_max_residual =
                    SMAX(body_max_residual, residual_norm / (global_reference_scale + TinyReal));
            }
        };
        auto update_magnetic_only_body_residual =
            [&](Real &body_max_residual,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_magnetic_residual_terms(
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    residual_norm, reference_norm);
                body_max_residual =
                    SMAX(body_max_residual, residual_norm / (global_reference_scale + TinyReal));
            }
        };

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        if (!use_frequency_operator_coil_air_only)
        {
            update_conductor_body_residual(breakdown.plate,
                                           plate_vector_potential_real,
                                           plate_vector_potential_imag,
                                           plate_electric_potential_gradient_real,
                                           plate_electric_potential_gradient_imag,
                                           plate_curl_nu_b_real,
                                           plate_curl_nu_b_imag,
                                           nullptr,
                                           nullptr,
                                           plate_sigma,
                                           plate_particles.TotalRealParticles());
        }
        if (use_frequency_coil_magnetic_only)
        {
            update_magnetic_only_body_residual(breakdown.coil,
                                               coil_curl_nu_b_real,
                                               coil_curl_nu_b_imag,
                                               coil_source_current_density_real,
                                               coil_source_current_density_imag,
                                               coil_particles.TotalRealParticles());
        }
        else
        {
            update_conductor_body_residual(breakdown.coil,
                                           coil_vector_potential_real,
                                           coil_vector_potential_imag,
                                           coil_electric_potential_gradient_real,
                                           coil_electric_potential_gradient_imag,
                                           coil_curl_nu_b_real,
                                           coil_curl_nu_b_imag,
                                           coil_source_current_density_real,
                                           coil_source_current_density_imag,
                                           coil_sigma,
                                           coil_particles.TotalRealParticles());
        }
        if (em_rate_include_air)
        {
            update_magnetic_only_body_residual(breakdown.air,
                                               air_curl_nu_b_real,
                                               air_curl_nu_b_imag,
                                               nullptr,
                                               nullptr,
                                               air_particles.TotalRealParticles());
        }
        return breakdown;
    };
    auto coil_source_tangential_direction = [&](const Vecd &position) -> Vecd
    {
        if (!use_circular_coil_source)
        {
            return coil_source_direction_used;
        }
        return extra_electromagnetics::TangentialDirectionAroundAxis(position, coil_center,
                                                                     coil_source_axis_used);
    };
    auto finalize_body_term_diagnostics = [&](EmBodyTermDiagnostics &body)
    {
        Real inv_total_volume = 1.0 / (body.total_volume + TinyReal);
        body.avg_source *= inv_total_volume;
        body.avg_curl_nu_b *= inv_total_volume;
        body.avg_sigma_grad_phi *= inv_total_volume;
        body.avg_omega_sigma_a *= inv_total_volume;
        body.avg_residual *= inv_total_volume;
    };
    auto point_box_distance = [](const Vecd &position, const Vec3d &center, const Vec3d &halfsize) -> Real
    {
        Vec3d distance = (position - center).cwiseAbs() - halfsize;
        Vec3d outside = distance.cwiseMax(static_cast<Real>(0.0));
        return outside.norm();
    };
    auto is_in_inner_box_shell = [&](const Vecd &position, const Vec3d &center,
                                     const Vec3d &halfsize, Real thickness) -> bool
    {
        Vec3d margin = halfsize - (position - center).cwiseAbs();
        if (margin.minCoeff() < static_cast<Real>(0.0))
        {
            return false;
        }
        return margin.minCoeff() <= thickness;
    };
    auto is_in_outer_box_shell = [&](const Vecd &position, const Vec3d &center,
                                     const Vec3d &halfsize, Real thickness) -> bool
    {
        Vec3d distance = (position - center).cwiseAbs() - halfsize;
        if (distance.maxCoeff() <= static_cast<Real>(0.0))
        {
            return false;
        }
        return point_box_distance(position, center, halfsize) <= thickness;
    };
    auto evaluate_coil_source_diagnostics = [&]() -> CoilSourceDiagnostics
    {
        CoilSourceDiagnostics diagnostics;
        diagnostics.path_length_si = coil_projection_length_si;
        Real configured_source_scale =
            apply_runtime_source_scale ? runtime_source_scale : static_cast<Real>(1.0);
        diagnostics.configured_ampere_turns_real =
            configured_source_scale * frequency_source_real_magnitude * coil_effective_area_used_si;
        diagnostics.configured_ampere_turns_imag =
            configured_source_scale * frequency_source_imag_magnitude * coil_effective_area_used_si;
        if (!use_frequency_aphi ||
            coil_source_current_density_real == nullptr ||
            coil_source_current_density_imag == nullptr)
        {
            return diagnostics;
        }

        Real source_density_volume_moment = 0.0;
        Real integrated_tangential_source_real = 0.0;
        Real integrated_tangential_source_imag = 0.0;
        size_t total_real_particles = coil_particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Real vol_si = coil_vol[i] * geom_volume_to_m3;
            Vecd tangent = coil_source_tangential_direction(coil_particle_positions[i]);
            Real projected_source_real = coil_source_current_density_real[i].dot(tangent);
            Real projected_source_imag = coil_source_current_density_imag[i].dot(tangent);
            Real source_density_magnitude =
                sqrt(coil_source_current_density_real[i].squaredNorm() +
                     coil_source_current_density_imag[i].squaredNorm());

            diagnostics.total_volume_si += vol_si;
            source_density_volume_moment += source_density_magnitude * vol_si;
            diagnostics.max_source_density_magnitude =
                SMAX(diagnostics.max_source_density_magnitude, source_density_magnitude);
            integrated_tangential_source_real += projected_source_real * vol_si;
            integrated_tangential_source_imag += projected_source_imag * vol_si;
        }

        diagnostics.avg_source_density_magnitude =
            source_density_volume_moment / (diagnostics.total_volume_si + TinyReal);
        if (diagnostics.path_length_si > TinyReal)
        {
            diagnostics.equivalent_ampere_turns_real =
                integrated_tangential_source_real / diagnostics.path_length_si;
            diagnostics.equivalent_ampere_turns_imag =
                integrated_tangential_source_imag / diagnostics.path_length_si;
            diagnostics.equivalent_ampere_turns_magnitude =
                sqrt(diagnostics.equivalent_ampere_turns_real *
                         diagnostics.equivalent_ampere_turns_real +
                     diagnostics.equivalent_ampere_turns_imag *
                         diagnostics.equivalent_ampere_turns_imag);
        }
        return diagnostics;
    };
    auto try_update_runtime_source_scale_by_circuit =
        [&](size_t current_outer_iteration, const CoilSourceDiagnostics &coil_source_diagnostics)
    {
        if (!enable_circuit_ni_feedback ||
            current_outer_iteration % circuit_ni_feedback_interval != 0)
        {
            return;
        }
        if (target_coil_ampere_turns_magnitude <= TinyReal)
        {
            return;
        }
        Real equivalent_ampere_turns =
            coil_source_diagnostics.equivalent_ampere_turns_magnitude;
        if (!(equivalent_ampere_turns > TinyReal) ||
            !std::isfinite(equivalent_ampere_turns))
        {
            return;
        }
        Real reference_to_simulated_ratio =
            target_coil_ampere_turns_magnitude /
            (equivalent_ampere_turns + TinyReal);
        if (!std::isfinite(reference_to_simulated_ratio))
        {
            return;
        }
        Real bounded_ratio =
            SMIN(circuit_ni_feedback_max_factor,
                 SMAX(circuit_ni_feedback_min_factor, reference_to_simulated_ratio));
        Real old_scale = runtime_source_scale;
        Real target_scale = old_scale * bounded_ratio;
        runtime_source_scale =
            old_scale + circuit_ni_feedback_gain * (target_scale - old_scale);
        runtime_source_scale =
            SMIN(runtime_source_scale_max,
                 SMAX(runtime_source_scale_min, runtime_source_scale));
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-circuit-feedback] iter=" << current_outer_iteration
                  << ", NI_target=" << target_coil_ampere_turns_magnitude
                  << ", NI_equivalent=" << equivalent_ampere_turns
                  << ", ratio_target_to_equivalent=" << reference_to_simulated_ratio
                  << ", bounded_ratio=" << bounded_ratio
                  << ", source_scale_old=" << old_scale
                  << ", source_scale_new=" << runtime_source_scale
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    };
    auto evaluate_em_term_breakdown = [&]() -> EmTermBreakdown
    {
        EmTermBreakdown breakdown;
        if (!use_frequency_aphi)
        {
            return breakdown;
        }

        auto update_conductor_body_terms =
            [&](EmBodyTermDiagnostics &body,
                Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                Real *volumetric_measure,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real weight = volumetric_measure[i];
                Real sigma_i = electrical_conductivity[i];
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Vecd rhs_real = source_real - curl_nu_b_real[i] -
                                sigma_i * electric_potential_gradient_real[i] +
                                sigma_i * harmonic_angular_frequency_runtime * vector_potential_imag[i];
                Vecd rhs_imag = source_imag - curl_nu_b_imag[i] -
                                sigma_i * electric_potential_gradient_imag[i] -
                                sigma_i * harmonic_angular_frequency_runtime * vector_potential_real[i];
                Real source_norm = sqrt(source_real.squaredNorm() + source_imag.squaredNorm());
                Real curl_norm = sqrt(curl_nu_b_real[i].squaredNorm() + curl_nu_b_imag[i].squaredNorm());
                Real grad_norm =
                    sigma_i * sqrt(electric_potential_gradient_real[i].squaredNorm() +
                                   electric_potential_gradient_imag[i].squaredNorm());
                Real omega_norm =
                    sigma_i * harmonic_angular_frequency_runtime *
                    sqrt(vector_potential_real[i].squaredNorm() +
                         vector_potential_imag[i].squaredNorm());
                Real residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());

                body.total_volume += weight;
                body.avg_source += source_norm * weight;
                body.avg_curl_nu_b += curl_norm * weight;
                body.avg_sigma_grad_phi += grad_norm * weight;
                body.avg_omega_sigma_a += omega_norm * weight;
                body.avg_residual += residual_norm * weight;
                body.max_source = SMAX(body.max_source, source_norm);
                body.max_curl_nu_b = SMAX(body.max_curl_nu_b, curl_norm);
                body.max_sigma_grad_phi = SMAX(body.max_sigma_grad_phi, grad_norm);
                body.max_omega_sigma_a = SMAX(body.max_omega_sigma_a, omega_norm);
                body.max_residual = SMAX(body.max_residual, residual_norm);
            }
            finalize_body_term_diagnostics(body);
        };
        auto update_magnetic_only_body_terms =
            [&](EmBodyTermDiagnostics &body,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *volumetric_measure,
                size_t total_particles)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                Real weight = volumetric_measure[i];
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Vecd rhs_real = source_real - curl_nu_b_real[i];
                Vecd rhs_imag = source_imag - curl_nu_b_imag[i];
                Real source_norm = sqrt(source_real.squaredNorm() + source_imag.squaredNorm());
                Real curl_norm = sqrt(curl_nu_b_real[i].squaredNorm() + curl_nu_b_imag[i].squaredNorm());
                Real residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());
                body.total_volume += weight;
                body.avg_source += source_norm * weight;
                body.avg_curl_nu_b += curl_norm * weight;
                body.avg_residual += residual_norm * weight;
                body.max_source = SMAX(body.max_source, source_norm);
                body.max_curl_nu_b = SMAX(body.max_curl_nu_b, curl_norm);
                body.max_residual = SMAX(body.max_residual, residual_norm);
            }
            finalize_body_term_diagnostics(body);
        };

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        if (!use_frequency_operator_coil_air_only)
        {
            update_conductor_body_terms(breakdown.plate,
                                        plate_vector_potential_real,
                                        plate_vector_potential_imag,
                                        plate_electric_potential_gradient_real,
                                        plate_electric_potential_gradient_imag,
                                        plate_curl_nu_b_real,
                                        plate_curl_nu_b_imag,
                                        nullptr,
                                        nullptr,
                                        plate_sigma,
                                        plate_vol,
                                        plate_particles.TotalRealParticles());
        }
        if (use_frequency_coil_magnetic_only)
        {
            update_magnetic_only_body_terms(breakdown.coil,
                                            coil_curl_nu_b_real,
                                            coil_curl_nu_b_imag,
                                            coil_source_current_density_real,
                                            coil_source_current_density_imag,
                                            coil_vol,
                                            coil_particles.TotalRealParticles());
        }
        else
        {
            update_conductor_body_terms(breakdown.coil,
                                        coil_vector_potential_real,
                                        coil_vector_potential_imag,
                                        coil_electric_potential_gradient_real,
                                        coil_electric_potential_gradient_imag,
                                        coil_curl_nu_b_real,
                                        coil_curl_nu_b_imag,
                                        coil_source_current_density_real,
                                        coil_source_current_density_imag,
                                        coil_sigma,
                                        coil_vol,
                                        coil_particles.TotalRealParticles());
        }
        if (em_rate_include_air)
        {
            update_magnetic_only_body_terms(breakdown.air,
                                            air_curl_nu_b_real,
                                            air_curl_nu_b_imag,
                                            nullptr,
                                            nullptr,
                                            air_vol,
                                            air_particles.TotalRealParticles());
        }
        return breakdown;
    };

    auto evaluate_em_interface_shell_breakdown = [&]() -> EmInterfaceShellBreakdown
    {
        EmInterfaceShellBreakdown breakdown;
        if (!use_frequency_aphi)
        {
            return breakdown;
        }

        auto update_conductor_shell_terms =
            [&](EmBodyTermDiagnostics &body,
                Vecd *positions,
                Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                Real *volumetric_measure,
                size_t total_particles,
                const auto &selector)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }

                Real weight = volumetric_measure[i];
                Real sigma_i = electrical_conductivity[i];
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Vecd rhs_real = source_real - curl_nu_b_real[i] -
                                sigma_i * electric_potential_gradient_real[i] +
                                sigma_i * harmonic_angular_frequency_runtime * vector_potential_imag[i];
                Vecd rhs_imag = source_imag - curl_nu_b_imag[i] -
                                sigma_i * electric_potential_gradient_imag[i] -
                                sigma_i * harmonic_angular_frequency_runtime * vector_potential_real[i];
                Real source_norm = sqrt(source_real.squaredNorm() + source_imag.squaredNorm());
                Real curl_norm = sqrt(curl_nu_b_real[i].squaredNorm() + curl_nu_b_imag[i].squaredNorm());
                Real grad_norm =
                    sigma_i * sqrt(electric_potential_gradient_real[i].squaredNorm() +
                                   electric_potential_gradient_imag[i].squaredNorm());
                Real omega_norm =
                    sigma_i * harmonic_angular_frequency_runtime *
                    sqrt(vector_potential_real[i].squaredNorm() +
                         vector_potential_imag[i].squaredNorm());
                Real residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());

                body.total_volume += weight;
                body.avg_source += source_norm * weight;
                body.avg_curl_nu_b += curl_norm * weight;
                body.avg_sigma_grad_phi += grad_norm * weight;
                body.avg_omega_sigma_a += omega_norm * weight;
                body.avg_residual += residual_norm * weight;
                body.max_source = SMAX(body.max_source, source_norm);
                body.max_curl_nu_b = SMAX(body.max_curl_nu_b, curl_norm);
                body.max_sigma_grad_phi = SMAX(body.max_sigma_grad_phi, grad_norm);
                body.max_omega_sigma_a = SMAX(body.max_omega_sigma_a, omega_norm);
                body.max_residual = SMAX(body.max_residual, residual_norm);
            }
            finalize_body_term_diagnostics(body);
        };
        auto update_magnetic_shell_terms =
            [&](EmBodyTermDiagnostics &body,
                Vecd *positions,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *volumetric_measure,
                size_t total_particles,
                const auto &selector)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }

                Real weight = volumetric_measure[i];
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Vecd rhs_real = source_real - curl_nu_b_real[i];
                Vecd rhs_imag = source_imag - curl_nu_b_imag[i];
                Real source_norm = sqrt(source_real.squaredNorm() + source_imag.squaredNorm());
                Real curl_norm = sqrt(curl_nu_b_real[i].squaredNorm() + curl_nu_b_imag[i].squaredNorm());
                Real residual_norm = sqrt(rhs_real.squaredNorm() + rhs_imag.squaredNorm());

                body.total_volume += weight;
                body.avg_source += source_norm * weight;
                body.avg_curl_nu_b += curl_norm * weight;
                body.avg_residual += residual_norm * weight;
                body.max_source = SMAX(body.max_source, source_norm);
                body.max_curl_nu_b = SMAX(body.max_curl_nu_b, curl_norm);
                body.max_residual = SMAX(body.max_residual, residual_norm);
            }
            finalize_body_term_diagnostics(body);
        };

        auto in_plate_inner_shell = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_inner_shell_thickness);
        };
        auto in_coil_inner_shell = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, coil_center, coil_halfsize,
                                         coil_interface_inner_shell_thickness);
        };
        auto in_plate_air_shell = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_air_shell_thickness);
        };
        auto in_coil_air_shell = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, coil_center, coil_halfsize,
                                         coil_interface_air_shell_thickness);
        };

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        if (!use_frequency_operator_coil_air_only)
        {
            update_conductor_shell_terms(breakdown.plate_inner,
                                         plate_positions,
                                         plate_vector_potential_real,
                                         plate_vector_potential_imag,
                                         plate_electric_potential_gradient_real,
                                         plate_electric_potential_gradient_imag,
                                         plate_curl_nu_b_real,
                                         plate_curl_nu_b_imag,
                                         nullptr,
                                         nullptr,
                                         plate_sigma,
                                         plate_vol,
                                         plate_particles.TotalRealParticles(),
                                         in_plate_inner_shell);
            update_magnetic_shell_terms(breakdown.plate_air,
                                        air_positions,
                                        air_curl_nu_b_real,
                                        air_curl_nu_b_imag,
                                        nullptr,
                                        nullptr,
                                        air_vol,
                                        air_particles.TotalRealParticles(),
                                        in_plate_air_shell);
        }
        if (use_frequency_coil_magnetic_only)
        {
            update_magnetic_shell_terms(breakdown.coil_inner,
                                        coil_positions,
                                        coil_curl_nu_b_real,
                                        coil_curl_nu_b_imag,
                                        coil_source_current_density_real,
                                        coil_source_current_density_imag,
                                        coil_vol,
                                        coil_particles.TotalRealParticles(),
                                        in_coil_inner_shell);
        }
        else
        {
            update_conductor_shell_terms(breakdown.coil_inner,
                                         coil_positions,
                                         coil_vector_potential_real,
                                         coil_vector_potential_imag,
                                         coil_electric_potential_gradient_real,
                                         coil_electric_potential_gradient_imag,
                                         coil_curl_nu_b_real,
                                         coil_curl_nu_b_imag,
                                         coil_source_current_density_real,
                                         coil_source_current_density_imag,
                                         coil_sigma,
                                         coil_vol,
                                         coil_particles.TotalRealParticles(),
                                         in_coil_inner_shell);
        }
        update_magnetic_shell_terms(breakdown.coil_air,
                                    air_positions,
                                    air_curl_nu_b_real,
                                    air_curl_nu_b_imag,
                                    nullptr,
                                    nullptr,
                                    air_vol,
                                    air_particles.TotalRealParticles(),
                                    in_coil_air_shell);
        return breakdown;
    };
    auto evaluate_em_coil_air_interface_ab_diagnostics = [&]() -> EmCoilAirInterfaceShellFieldDiagnostics
    {
        EmCoilAirInterfaceShellFieldDiagnostics out;
        if (!use_frequency_aphi || !enable_coil_air_interface_b_diag)
        {
            return out;
        }
        if (coil_vector_potential_curl_real == nullptr || air_vector_potential_curl_real == nullptr ||
            coil_vector_potential_real == nullptr || air_vector_potential_real == nullptr)
        {
            return out;
        }

        auto in_coil_inner_shell_for_ab = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, coil_center, coil_halfsize,
                                         coil_interface_inner_shell_thickness);
        };
        auto in_coil_air_shell_for_ab = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, coil_center, coil_halfsize,
                                          coil_interface_air_shell_thickness);
        };

        auto accumulate_ab_shell = [&](Vecd *positions,
                                       Real *volumetric_measure,
                                       size_t total_particles,
                                       const auto &selector,
                                       Vecd *vector_potential_real,
                                       Vecd *vector_potential_imag,
                                       AngularVecd *curl_real,
                                       AngularVecd *curl_imag,
                                       Real &shell_volume,
                                       Real &sum_a_mag_weighted,
                                       Real &sum_b_mag_weighted)
        {
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }
                const Real weight = volumetric_measure[i];
                shell_volume += weight;

                Real a_mag_sq = vector_potential_real[i].squaredNorm();
                if (vector_potential_imag != nullptr)
                {
                    a_mag_sq += vector_potential_imag[i].squaredNorm();
                }
                sum_a_mag_weighted += sqrt(a_mag_sq) * weight;

                const Vec3d b_real(curl_real[i][0], curl_real[i][1], curl_real[i][2]);
                Real b_mag_sq = b_real.squaredNorm();
                if (curl_imag != nullptr)
                {
                    const Vec3d b_imag(curl_imag[i][0], curl_imag[i][1], curl_imag[i][2]);
                    b_mag_sq += b_imag.squaredNorm();
                }
                sum_b_mag_weighted += sqrt(b_mag_sq) * weight;
            }
        };

        Real sum_a_coil_inner = 0.0;
        Real sum_b_coil_inner = 0.0;
        Real sum_a_air = 0.0;
        Real sum_b_air = 0.0;
        accumulate_ab_shell(coil_positions,
                            coil_vol,
                            coil_particles.TotalRealParticles(),
                            in_coil_inner_shell_for_ab,
                            coil_vector_potential_real,
                            coil_vector_potential_imag,
                            coil_vector_potential_curl_real,
                            coil_vector_potential_curl_imag,
                            out.coil_inner_volume,
                            sum_a_coil_inner,
                            sum_b_coil_inner);
        accumulate_ab_shell(air_positions,
                            air_vol,
                            air_particles.TotalRealParticles(),
                            in_coil_air_shell_for_ab,
                            air_vector_potential_real,
                            air_vector_potential_imag,
                            air_vector_potential_curl_real,
                            air_vector_potential_curl_imag,
                            out.coil_air_volume,
                            sum_a_air,
                            sum_b_air);

        if (out.coil_inner_volume > TinyReal)
        {
            const Real inv_v = 1.0 / out.coil_inner_volume;
            out.coil_inner_avg_a_mag = sum_a_coil_inner * inv_v;
            out.coil_inner_avg_b_mag = sum_b_coil_inner * inv_v;
        }
        if (out.coil_air_volume > TinyReal)
        {
            const Real inv_v = 1.0 / out.coil_air_volume;
            out.coil_air_avg_a_mag = sum_a_air * inv_v;
            out.coil_air_avg_b_mag = sum_b_air * inv_v;
        }
        return out;
    };
    auto evaluate_em_magnetic_diagonal_breakdown = [&]() -> EmMagneticDiagonalBreakdown
    {
        EmMagneticDiagonalBreakdown breakdown;
        if (!use_frequency_aphi || !enable_contact_magnetic_diagonal_diagnostics)
        {
            return breakdown;
        }

        auto finalize_diagonal_diagnostics = [&](EmMagneticDiagonalBodyDiagnostics &body)
        {
            if (body.total_volume <= TinyReal)
            {
                return;
            }
            body.avg_local_diagonal /= body.total_volume;
            body.avg_contact_diagonal /= body.total_volume;
            body.avg_jacobi_diagonal /= body.total_volume;
            body.avg_conservative_diagonal /= body.total_volume;
            body.avg_balanced_diagonal /= body.total_volume;
            body.avg_contact_ratio /= body.total_volume;
            body.avg_conservative_to_jacobi /= body.total_volume;
            body.avg_balanced_to_jacobi /= body.total_volume;
        };
        auto collect_contact_material_fields =
            [&](BaseContactRelation &contact_relation,
                StdVec<Real *> &contact_volumes,
                StdVec<Real *> &contact_reluctivities)
        {
            contact_volumes.clear();
            contact_reluctivities.clear();
            StdVec<BaseParticles *> contact_particles = contact_relation.getContactParticles();
            contact_volumes.reserve(contact_particles.size());
            contact_reluctivities.reserve(contact_particles.size());
            for (BaseParticles *contact_particle : contact_particles)
            {
                contact_volumes.push_back(
                    contact_particle->getVariableDataByName<Real>("VolumetricMeasure"));
                contact_reluctivities.push_back(
                    contact_particle->getVariableDataByName<Real>("MagneticReluctivity"));
            }
        };
        auto update_body_diagonal_terms =
            [&](EmMagneticDiagonalBodyDiagnostics &body,
                Vecd *positions,
                Real *volumetric_measure,
                Real *magnetic_reluctivity,
                size_t total_particles,
                BaseInnerRelation &inner_relation,
                BaseContactRelation &contact_relation,
                const auto &selector)
        {
            StdVec<Real *> contact_volumes;
            StdVec<Real *> contact_reluctivities;
            collect_contact_material_fields(contact_relation, contact_volumes, contact_reluctivities);
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }
                Real nu_i = magnetic_reluctivity[i];
                Real local_diagonal = 0.0;
                Real contact_diagonal = 0.0;
                Real local_row_sum = 0.0;
                Real contact_row_sum = 0.0;
                Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (inner_neighborhood.r_ij_[n] < TinyReal)
                    {
                        continue;
                    }
                    Vecd gradW_ijV_j =
                        inner_neighborhood.dW_ij_[n] * volumetric_measure[index_j] *
                        inner_neighborhood.e_ij_[n];
                    Real grad_norm_sq = gradW_ijV_j.squaredNorm();
                    if (!std::isfinite(grad_norm_sq))
                    {
                        continue;
                    }
                    Real nu_ij = 0.5 * (nu_i + magnetic_reluctivity[index_j]);
                    local_diagonal +=
                        frequency_magnetic_diagonal_scaling * nu_ij * grad_norm_sq;
                    local_row_sum +=
                        sqrt(SMAX(static_cast<Real>(0.0), nu_ij)) * sqrt(grad_norm_sq);
                }
                for (size_t k = 0; k != contact_relation.contact_configuration_.size(); ++k)
                {
                    Real *contact_vol_k = contact_volumes[k];
                    Real *contact_nu_k = contact_reluctivities[k];
                    Neighborhood &contact_neighborhood = contact_relation.contact_configuration_[k][i];
                    for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = contact_neighborhood.j_[n];
                        if (contact_neighborhood.r_ij_[n] < TinyReal)
                        {
                            continue;
                        }
                        Vecd gradW_ijV_j =
                            contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j] *
                            contact_neighborhood.e_ij_[n];
                        Real grad_norm_sq = gradW_ijV_j.squaredNorm();
                        if (!std::isfinite(grad_norm_sq))
                        {
                            continue;
                        }
                        Real nu_ij = 0.5 * (nu_i + contact_nu_k[index_j]);
                        contact_diagonal +=
                            frequency_magnetic_diagonal_scaling * nu_ij * grad_norm_sq;
                        contact_row_sum +=
                            sqrt(SMAX(static_cast<Real>(0.0), nu_ij)) * sqrt(grad_norm_sq);
                    }
                }
                Real jacobi_diagonal = local_diagonal + contact_diagonal;
                Real conservative_diagonal =
                    frequency_magnetic_diagonal_scaling *
                    (local_row_sum + contact_row_sum) * (local_row_sum + contact_row_sum);
                Real safe_jacobi = SMAX(jacobi_diagonal, TinyReal);
                Real safe_conservative = SMAX(conservative_diagonal, safe_jacobi);
                Real balanced_diagonal = SMAX(
                    safe_jacobi, static_cast<Real>(sqrt(safe_jacobi * safe_conservative)));
                Real contact_ratio =
                    jacobi_diagonal > TinyReal ? contact_diagonal / (jacobi_diagonal + TinyReal) : 0.0;
                Real weight = volumetric_measure[i];
                body.total_volume += weight;
                body.avg_local_diagonal += local_diagonal * weight;
                body.avg_contact_diagonal += contact_diagonal * weight;
                body.avg_jacobi_diagonal += jacobi_diagonal * weight;
                body.avg_conservative_diagonal += conservative_diagonal * weight;
                body.avg_balanced_diagonal += balanced_diagonal * weight;
                body.avg_contact_ratio += contact_ratio * weight;
                body.avg_conservative_to_jacobi +=
                    (conservative_diagonal / (safe_jacobi + TinyReal)) * weight;
                body.avg_balanced_to_jacobi +=
                    (balanced_diagonal / (safe_jacobi + TinyReal)) * weight;
                body.max_contact_ratio = SMAX(body.max_contact_ratio, contact_ratio);
            }
            finalize_diagonal_diagnostics(body);
        };

        auto in_plate_inner_shell = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_inner_shell_thickness);
        };
        auto in_coil_inner_shell = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, coil_center, coil_halfsize,
                                         coil_interface_inner_shell_thickness);
        };
        auto in_plate_air_shell = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_air_shell_thickness);
        };
        auto in_coil_air_shell = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, coil_center, coil_halfsize,
                                         coil_interface_air_shell_thickness);
        };
        auto all_particles = [&](const Vecd &) -> bool
        {
            return true;
        };

        if (!use_frequency_operator_coil_air_only)
        {
            update_body_diagonal_terms(breakdown.plate,
                                       plate_positions,
                                       plate_vol,
                                       plate_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                       plate_particles.TotalRealParticles(),
                                       plate_inner,
                                       plate_contact,
                                       all_particles);
            update_body_diagonal_terms(breakdown.plate_inner,
                                       plate_positions,
                                       plate_vol,
                                       plate_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                       plate_particles.TotalRealParticles(),
                                       plate_inner,
                                       plate_contact,
                                       in_plate_inner_shell);
            update_body_diagonal_terms(breakdown.plate_air,
                                       air_positions,
                                       air_vol,
                                       air_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                       air_particles.TotalRealParticles(),
                                       air_inner,
                                       air_contact,
                                       in_plate_air_shell);
        }

        update_body_diagonal_terms(breakdown.coil,
                                   coil_positions,
                                   coil_vol,
                                   coil_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                   coil_particles.TotalRealParticles(),
                                   coil_inner,
                                   coil_contact,
                                   all_particles);
        update_body_diagonal_terms(breakdown.coil_inner,
                                   coil_positions,
                                   coil_vol,
                                   coil_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                   coil_particles.TotalRealParticles(),
                                   coil_inner,
                                   coil_contact,
                                   in_coil_inner_shell);
        update_body_diagonal_terms(breakdown.coil_air,
                                   air_positions,
                                   air_vol,
                                   air_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                   air_particles.TotalRealParticles(),
                                   air_inner,
                                   air_contact,
                                   in_coil_air_shell);
        update_body_diagonal_terms(breakdown.air,
                                   air_positions,
                                   air_vol,
                                   air_particles.getVariableDataByName<Real>("MagneticReluctivity"),
                                   air_particles.TotalRealParticles(),
                                   air_inner,
                                   air_contact,
                                   all_particles);
        return breakdown;
    };
    auto write_coil_air_transfer_diagnostics = [&]()
    {
        if (!use_frequency_aphi || !enable_coil_air_transfer_diagnostics)
        {
            return;
        }

        struct TransferBucket
        {
            std::string name;
            size_t particle_count = 0;
            Real total_volume_si = 0.0;
            Real sum_distance = 0.0;
            Real sum_source = 0.0;
            Real sum_curl_nu_b = 0.0;
            Real sum_residual = 0.0;
            Real sum_a_magnitude = 0.0;
            Real sum_b_magnitude = 0.0;
            Real max_residual = 0.0;
        };
        auto accumulate_bucket =
            [&](TransferBucket &bucket,
                Real distance,
                Real volume_si,
                Real source_norm,
                Real curl_norm,
                Real residual_norm,
                Real a_magnitude,
                Real b_magnitude)
        {
            if (!(volume_si > TinyReal))
            {
                return;
            }
            bucket.particle_count++;
            bucket.total_volume_si += volume_si;
            bucket.sum_distance += distance * volume_si;
            bucket.sum_source += source_norm * volume_si;
            bucket.sum_curl_nu_b += curl_norm * volume_si;
            bucket.sum_residual += residual_norm * volume_si;
            bucket.sum_a_magnitude += a_magnitude * volume_si;
            bucket.sum_b_magnitude += b_magnitude * volume_si;
            bucket.max_residual = SMAX(bucket.max_residual, residual_norm);
        };
        auto avg_from_bucket = [&](Real weighted_sum, const TransferBucket &bucket) -> Real
        {
            return weighted_sum / (bucket.total_volume_si + TinyReal);
        };
        auto signed_distance_to_coil_box = [&](const Vecd &position) -> Real
        {
            Vec3d distance = (position - coil_center).cwiseAbs() - coil_halfsize;
            if (distance.maxCoeff() <= static_cast<Real>(0.0))
            {
                Vec3d margin = coil_halfsize - (position - coil_center).cwiseAbs();
                return -margin.minCoeff();
            }
            return point_box_distance(position, coil_center, coil_halfsize);
        };
        auto vector_magnitude = [&](const Vecd &real_part, const Vecd &imag_part) -> Real
        {
            return sqrt(real_part.squaredNorm() + imag_part.squaredNorm());
        };

        TransferBucket coil_core;
        coil_core.name = "coil_core";
        TransferBucket coil_inner_shell;
        coil_inner_shell.name = "coil_inner_shell";
        TransferBucket air_near_shell;
        air_near_shell.name = "air_near_shell";
        TransferBucket air_far;
        air_far.name = "air_far";

        Real *coil_sigma = coil_particles.getVariableDataByName<Real>("ElectricalConductivity");
        size_t total_coil_particles = coil_particles.TotalRealParticles();
        for (size_t i = 0; i != total_coil_particles; ++i)
        {
            Real distance = signed_distance_to_coil_box(coil_positions[i]);
            TransferBucket *bucket = nullptr;
            if (distance <= -coil_interface_inner_shell_thickness)
            {
                bucket = &coil_core;
            }
            else
            {
                bucket = &coil_inner_shell;
            }

            Vecd source_real = coil_source_current_density_real != nullptr
                                   ? coil_source_current_density_real[i]
                                   : ZeroData<Vecd>::value;
            Vecd source_imag = coil_source_current_density_imag != nullptr
                                   ? coil_source_current_density_imag[i]
                                   : ZeroData<Vecd>::value;
            Vecd curl_real = coil_curl_nu_b_real != nullptr
                                 ? coil_curl_nu_b_real[i]
                                 : ZeroData<Vecd>::value;
            Vecd curl_imag = coil_curl_nu_b_imag != nullptr
                                 ? coil_curl_nu_b_imag[i]
                                 : ZeroData<Vecd>::value;
            Vecd rhs_real = source_real - curl_real;
            Vecd rhs_imag = source_imag - curl_imag;
            if (!use_frequency_coil_magnetic_only &&
                coil_vector_potential_real != nullptr &&
                coil_vector_potential_imag != nullptr &&
                coil_electric_potential_gradient_real != nullptr &&
                coil_electric_potential_gradient_imag != nullptr &&
                coil_sigma != nullptr)
            {
                Real sigma_i = coil_sigma[i];
                rhs_real -= sigma_i * coil_electric_potential_gradient_real[i];
                rhs_real += sigma_i * harmonic_angular_frequency_runtime *
                            coil_vector_potential_imag[i];
                rhs_imag -= sigma_i * coil_electric_potential_gradient_imag[i];
                rhs_imag -= sigma_i * harmonic_angular_frequency_runtime *
                            coil_vector_potential_real[i];
            }

            Vecd a_real = coil_vector_potential_real != nullptr
                              ? coil_vector_potential_real[i]
                              : ZeroData<Vecd>::value;
            Vecd a_imag = coil_vector_potential_imag != nullptr
                              ? coil_vector_potential_imag[i]
                              : ZeroData<Vecd>::value;
            Vecd b_real = coil_vector_potential_curl_real != nullptr
                              ? Vecd(coil_vector_potential_curl_real[i][0],
                                     coil_vector_potential_curl_real[i][1],
                                     coil_vector_potential_curl_real[i][2])
                              : ZeroData<Vecd>::value;
            Vecd b_imag = coil_vector_potential_curl_imag != nullptr
                              ? Vecd(coil_vector_potential_curl_imag[i][0],
                                     coil_vector_potential_curl_imag[i][1],
                                     coil_vector_potential_curl_imag[i][2])
                              : ZeroData<Vecd>::value;

            Real source_norm = vector_magnitude(source_real, source_imag);
            Real curl_norm = vector_magnitude(curl_real, curl_imag);
            Real residual_norm = vector_magnitude(rhs_real, rhs_imag);
            Real a_magnitude = vector_magnitude(a_real, a_imag);
            Real b_magnitude = vector_magnitude(b_real, b_imag);
            Real volume_si = coil_vol[i] * geom_volume_to_m3;
            accumulate_bucket(*bucket, distance, volume_si, source_norm, curl_norm,
                              residual_norm, a_magnitude, b_magnitude);
        }

        size_t total_air_particles = air_particles.TotalRealParticles();
        for (size_t i = 0; i != total_air_particles; ++i)
        {
            Real distance = signed_distance_to_coil_box(air_positions[i]);
            if (distance <= static_cast<Real>(0.0))
            {
                continue;
            }
            TransferBucket *bucket =
                (distance <= coil_interface_air_shell_thickness)
                    ? &air_near_shell
                    : &air_far;

            Vecd curl_real = air_curl_nu_b_real != nullptr
                                 ? air_curl_nu_b_real[i]
                                 : ZeroData<Vecd>::value;
            Vecd curl_imag = air_curl_nu_b_imag != nullptr
                                 ? air_curl_nu_b_imag[i]
                                 : ZeroData<Vecd>::value;
            Vecd a_real = air_vector_potential_real != nullptr
                              ? air_vector_potential_real[i]
                              : ZeroData<Vecd>::value;
            Vecd a_imag = air_vector_potential_imag != nullptr
                              ? air_vector_potential_imag[i]
                              : ZeroData<Vecd>::value;
            Vecd b_real = air_vector_potential_curl_real != nullptr
                              ? Vecd(air_vector_potential_curl_real[i][0],
                                     air_vector_potential_curl_real[i][1],
                                     air_vector_potential_curl_real[i][2])
                              : ZeroData<Vecd>::value;
            Vecd b_imag = air_vector_potential_curl_imag != nullptr
                              ? Vecd(air_vector_potential_curl_imag[i][0],
                                     air_vector_potential_curl_imag[i][1],
                                     air_vector_potential_curl_imag[i][2])
                              : ZeroData<Vecd>::value;

            Real source_norm = 0.0;
            Real curl_norm = vector_magnitude(curl_real, curl_imag);
            Real residual_norm = curl_norm;
            Real a_magnitude = vector_magnitude(a_real, a_imag);
            Real b_magnitude = vector_magnitude(b_real, b_imag);
            Real volume_si = air_vol[i] * geom_volume_to_m3;
            accumulate_bucket(*bucket, distance, volume_si, source_norm, curl_norm,
                              residual_norm, a_magnitude, b_magnitude);
        }

        std::ofstream transfer_summary_file(
            io_environment.OutputFolder() + "/team7_coil_air_transfer_summary.csv",
            std::ios::out | std::ios::trunc);
        transfer_summary_file << std::setprecision(12);
        transfer_summary_file
            << "region,particle_count,total_volume_si,avg_signed_distance,"
            << "avg_source,avg_curl_nu_b,avg_residual,max_residual,avg_a_mag,avg_b_mag\n";
        auto write_summary_row = [&](const TransferBucket &bucket)
        {
            transfer_summary_file << bucket.name << ","
                                  << bucket.particle_count << ","
                                  << bucket.total_volume_si << ","
                                  << avg_from_bucket(bucket.sum_distance, bucket) << ","
                                  << avg_from_bucket(bucket.sum_source, bucket) << ","
                                  << avg_from_bucket(bucket.sum_curl_nu_b, bucket) << ","
                                  << avg_from_bucket(bucket.sum_residual, bucket) << ","
                                  << bucket.max_residual << ","
                                  << avg_from_bucket(bucket.sum_a_magnitude, bucket) << ","
                                  << avg_from_bucket(bucket.sum_b_magnitude, bucket) << "\n";
        };
        write_summary_row(coil_core);
        write_summary_row(coil_inner_shell);
        write_summary_row(air_near_shell);
        write_summary_row(air_far);
        transfer_summary_file.flush();

        std::vector<TransferBucket> profile_bins(coil_air_transfer_profile_bins);
        Real profile_distance_min = -coil_interface_inner_shell_thickness;
        Real profile_distance_max = coil_interface_air_shell_thickness;
        Real profile_distance_span =
            SMAX(profile_distance_max - profile_distance_min, TinyReal);
        for (size_t i = 0; i != coil_air_transfer_profile_bins; ++i)
        {
            profile_bins[i].name = "bin_" + std::to_string(i);
        }
        auto accumulate_profile_bin =
            [&](Real distance,
                Real volume_si,
                Real source_norm,
                Real curl_norm,
                Real residual_norm,
                Real a_magnitude,
                Real b_magnitude)
        {
            if (distance < profile_distance_min || distance > profile_distance_max)
            {
                return;
            }
            Real normalized = (distance - profile_distance_min) / profile_distance_span;
            size_t bin_id = static_cast<size_t>(normalized *
                                                static_cast<Real>(coil_air_transfer_profile_bins));
            if (bin_id >= coil_air_transfer_profile_bins)
            {
                bin_id = coil_air_transfer_profile_bins - 1;
            }
            accumulate_bucket(profile_bins[bin_id], distance, volume_si, source_norm,
                              curl_norm, residual_norm, a_magnitude, b_magnitude);
        };

        for (size_t i = 0; i != total_coil_particles; ++i)
        {
            Real distance = signed_distance_to_coil_box(coil_positions[i]);
            Vecd source_real = coil_source_current_density_real != nullptr
                                   ? coil_source_current_density_real[i]
                                   : ZeroData<Vecd>::value;
            Vecd source_imag = coil_source_current_density_imag != nullptr
                                   ? coil_source_current_density_imag[i]
                                   : ZeroData<Vecd>::value;
            Vecd curl_real = coil_curl_nu_b_real != nullptr
                                 ? coil_curl_nu_b_real[i]
                                 : ZeroData<Vecd>::value;
            Vecd curl_imag = coil_curl_nu_b_imag != nullptr
                                 ? coil_curl_nu_b_imag[i]
                                 : ZeroData<Vecd>::value;
            Vecd rhs_real = source_real - curl_real;
            Vecd rhs_imag = source_imag - curl_imag;
            if (!use_frequency_coil_magnetic_only &&
                coil_vector_potential_real != nullptr &&
                coil_vector_potential_imag != nullptr &&
                coil_electric_potential_gradient_real != nullptr &&
                coil_electric_potential_gradient_imag != nullptr &&
                coil_sigma != nullptr)
            {
                Real sigma_i = coil_sigma[i];
                rhs_real -= sigma_i * coil_electric_potential_gradient_real[i];
                rhs_real += sigma_i * harmonic_angular_frequency_runtime *
                            coil_vector_potential_imag[i];
                rhs_imag -= sigma_i * coil_electric_potential_gradient_imag[i];
                rhs_imag -= sigma_i * harmonic_angular_frequency_runtime *
                            coil_vector_potential_real[i];
            }
            Vecd a_real = coil_vector_potential_real != nullptr
                              ? coil_vector_potential_real[i]
                              : ZeroData<Vecd>::value;
            Vecd a_imag = coil_vector_potential_imag != nullptr
                              ? coil_vector_potential_imag[i]
                              : ZeroData<Vecd>::value;
            Vecd b_real = coil_vector_potential_curl_real != nullptr
                              ? Vecd(coil_vector_potential_curl_real[i][0],
                                     coil_vector_potential_curl_real[i][1],
                                     coil_vector_potential_curl_real[i][2])
                              : ZeroData<Vecd>::value;
            Vecd b_imag = coil_vector_potential_curl_imag != nullptr
                              ? Vecd(coil_vector_potential_curl_imag[i][0],
                                     coil_vector_potential_curl_imag[i][1],
                                     coil_vector_potential_curl_imag[i][2])
                              : ZeroData<Vecd>::value;
            accumulate_profile_bin(distance,
                                   coil_vol[i] * geom_volume_to_m3,
                                   vector_magnitude(source_real, source_imag),
                                   vector_magnitude(curl_real, curl_imag),
                                   vector_magnitude(rhs_real, rhs_imag),
                                   vector_magnitude(a_real, a_imag),
                                   vector_magnitude(b_real, b_imag));
        }
        for (size_t i = 0; i != total_air_particles; ++i)
        {
            Real distance = signed_distance_to_coil_box(air_positions[i]);
            Vecd curl_real = air_curl_nu_b_real != nullptr
                                 ? air_curl_nu_b_real[i]
                                 : ZeroData<Vecd>::value;
            Vecd curl_imag = air_curl_nu_b_imag != nullptr
                                 ? air_curl_nu_b_imag[i]
                                 : ZeroData<Vecd>::value;
            Vecd a_real = air_vector_potential_real != nullptr
                              ? air_vector_potential_real[i]
                              : ZeroData<Vecd>::value;
            Vecd a_imag = air_vector_potential_imag != nullptr
                              ? air_vector_potential_imag[i]
                              : ZeroData<Vecd>::value;
            Vecd b_real = air_vector_potential_curl_real != nullptr
                              ? Vecd(air_vector_potential_curl_real[i][0],
                                     air_vector_potential_curl_real[i][1],
                                     air_vector_potential_curl_real[i][2])
                              : ZeroData<Vecd>::value;
            Vecd b_imag = air_vector_potential_curl_imag != nullptr
                              ? Vecd(air_vector_potential_curl_imag[i][0],
                                     air_vector_potential_curl_imag[i][1],
                                     air_vector_potential_curl_imag[i][2])
                              : ZeroData<Vecd>::value;
            Real curl_norm = vector_magnitude(curl_real, curl_imag);
            accumulate_profile_bin(distance,
                                   air_vol[i] * geom_volume_to_m3,
                                   0.0,
                                   curl_norm,
                                   curl_norm,
                                   vector_magnitude(a_real, a_imag),
                                   vector_magnitude(b_real, b_imag));
        }

        std::ofstream transfer_profile_file(
            io_environment.OutputFolder() + "/team7_coil_air_transfer_profile.csv",
            std::ios::out | std::ios::trunc);
        transfer_profile_file << std::setprecision(12);
        transfer_profile_file
            << "bin_id,distance_min,distance_max,distance_center,particle_count,total_volume_si,"
            << "avg_signed_distance,avg_source,avg_curl_nu_b,avg_residual,max_residual,avg_a_mag,avg_b_mag\n";
        for (size_t i = 0; i != profile_bins.size(); ++i)
        {
            Real bin_min = profile_distance_min +
                           static_cast<Real>(i) * profile_distance_span /
                               static_cast<Real>(profile_bins.size());
            Real bin_max = profile_distance_min +
                           static_cast<Real>(i + 1) * profile_distance_span /
                               static_cast<Real>(profile_bins.size());
            Real bin_center = 0.5 * (bin_min + bin_max);
            const TransferBucket &bucket = profile_bins[i];
            transfer_profile_file << i << ","
                                  << bin_min << ","
                                  << bin_max << ","
                                  << bin_center << ","
                                  << bucket.particle_count << ","
                                  << bucket.total_volume_si << ","
                                  << avg_from_bucket(bucket.sum_distance, bucket) << ","
                                  << avg_from_bucket(bucket.sum_source, bucket) << ","
                                  << avg_from_bucket(bucket.sum_curl_nu_b, bucket) << ","
                                  << avg_from_bucket(bucket.sum_residual, bucket) << ","
                                  << bucket.max_residual << ","
                                  << avg_from_bucket(bucket.sum_a_magnitude, bucket) << ","
                                  << avg_from_bucket(bucket.sum_b_magnitude, bucket) << "\n";
        }
        transfer_profile_file.flush();

        Real coil_inner_b =
            avg_from_bucket(coil_inner_shell.sum_b_magnitude, coil_inner_shell);
        Real air_near_b =
            avg_from_bucket(air_near_shell.sum_b_magnitude, air_near_shell);
        Real transfer_ratio = air_near_b / (coil_inner_b + TinyReal);
        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-coil-air-transfer] avg_B_coil_inner=" << coil_inner_b
                  << ", avg_B_air_near=" << air_near_b
                  << ", transfer_ratio_air_over_coil=" << transfer_ratio
                  << ", inner_shell_thickness=" << coil_interface_inner_shell_thickness
                  << ", air_shell_thickness=" << coil_interface_air_shell_thickness
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    };
    bool plate_component_hessian_inner_ck_prepared = false;
    auto ensure_plate_component_hessian_ck_prepared = [&]()
    {
        if (plate_component_hessian_inner_ck_prepared)
        {
            return;
        }
        update_plate_cell_linked_list_ck.exec();
        update_plate_inner_ck.exec();
        plate_verify_linear_correction_ck.exec();
        plate_verify_displacement_matrix_gradient_ck.exec();
        plate_verify_hessian_correction_ck.exec();
        plate_component_hessian_inner_ck_prepared = true;
    };
    bool plate_component_hessian_contact_ck_prepared = false;
    auto ensure_plate_component_hessian_contact_ck_prepared = [&]()
    {
        if (plate_component_hessian_contact_ck_prepared)
        {
            return;
        }
        update_plate_cell_linked_list_ck.exec();
        update_plate_inner_ck.exec();
        update_plate_contact_ck.exec();
        plate_component_linear_correction_ck.exec();
        plate_component_displacement_matrix_gradient_ck.exec();
        plate_component_hessian_correction_ck.exec();
        plate_component_hessian_contact_ck_prepared = true;
    };
    bool coil_component_hessian_inner_ck_prepared = false;
    auto ensure_coil_component_hessian_ck_prepared = [&]()
    {
        if (coil_component_hessian_inner_ck_prepared)
        {
            return;
        }
        update_coil_cell_linked_list_ck.exec();
        update_coil_inner_ck.exec();
        coil_verify_linear_correction_ck.exec();
        coil_verify_displacement_matrix_gradient_ck.exec();
        coil_verify_hessian_correction_ck.exec();
        coil_component_hessian_inner_ck_prepared = true;
    };
    bool coil_component_hessian_contact_ck_prepared = false;
    auto ensure_coil_component_hessian_contact_ck_prepared = [&]()
    {
        if (coil_component_hessian_contact_ck_prepared)
        {
            return;
        }
        update_coil_cell_linked_list_ck.exec();
        update_coil_inner_ck.exec();
        update_coil_contact_ck.exec();
        coil_component_linear_correction_ck.exec();
        coil_component_displacement_matrix_gradient_ck.exec();
        coil_component_hessian_correction_ck.exec();
        coil_component_hessian_contact_ck_prepared = true;
    };
    bool coil_component_hessian_runtime_warning_printed = false;
    size_t coil_component_hessian_runtime_guard_activations = 0;
    size_t coil_component_hessian_runtime_shell_blend_activations = 0;
    size_t coil_component_hessian_runtime_core_gate_activations = 0;
    auto refresh_frequency_air_scalar_gradients = [&]()
    {
        if (!use_frequency_air_scalar_potential)
        {
            return;
        }
        electric_potential_gradient_real_air_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_real_air_contact.exec();
        }
        electric_potential_gradient_imag_air_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_imag_air_contact.exec();
        }
    };
    auto refresh_frequency_coil_scalar_gradients = [&]()
    {
        if (!use_frequency_coil_scalar_potential)
        {
            return;
        }
        electric_potential_gradient_real_coil_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_real_coil_contact.exec();
        }
        electric_potential_gradient_imag_coil_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_imag_coil_contact.exec();
        }
    };
    auto refresh_frequency_coil_magnetic_inner_by_component_hessian = [&]()
    {
        ensure_coil_component_hessian_ck_prepared();
        copy_coil_a_real_components.exec();
        coil_a_real_gradient_x_ck.exec();
        coil_a_real_hessian_x_ck.exec();
        coil_a_real_gradient_y_ck.exec();
        coil_a_real_hessian_y_ck.exec();
        coil_a_real_gradient_z_ck.exec();
        coil_a_real_hessian_z_ck.exec();
        reconstruct_coil_curl_nu_b_real_from_components.exec();

        copy_coil_a_imag_components.exec();
        coil_a_imag_gradient_x_ck.exec();
        coil_a_imag_hessian_x_ck.exec();
        coil_a_imag_gradient_y_ck.exec();
        coil_a_imag_hessian_y_ck.exec();
        coil_a_imag_gradient_z_ck.exec();
        coil_a_imag_hessian_z_ck.exec();
        reconstruct_coil_curl_nu_b_imag_from_components.exec();
    };
    auto refresh_frequency_coil_magnetic_contact_by_component_hessian = [&]()
    {
        ensure_coil_component_hessian_contact_ck_prepared();

        copy_coil_a_real_components.exec();
        copy_plate_a_real_components_for_coil.exec();
        copy_air_a_real_components_for_coil.exec();
        coil_a_real_gradient_x_contact_ck.exec();
        coil_a_real_hessian_x_contact_ck.exec();
        coil_a_real_gradient_y_contact_ck.exec();
        coil_a_real_hessian_y_contact_ck.exec();
        coil_a_real_gradient_z_contact_ck.exec();
        coil_a_real_hessian_z_contact_ck.exec();
        reconstruct_coil_curl_nu_b_real_from_components.exec();

        copy_coil_a_imag_components.exec();
        copy_plate_a_imag_components_for_coil.exec();
        copy_air_a_imag_components_for_coil.exec();
        coil_a_imag_gradient_x_contact_ck.exec();
        coil_a_imag_hessian_x_contact_ck.exec();
        coil_a_imag_gradient_y_contact_ck.exec();
        coil_a_imag_hessian_y_contact_ck.exec();
        coil_a_imag_gradient_z_contact_ck.exec();
        coil_a_imag_hessian_z_contact_ck.exec();
        reconstruct_coil_curl_nu_b_imag_from_components.exec();
    };
    auto refresh_frequency_coil_magnetic_terms = [&]()
    {
        vector_potential_curl_real_coil_inner.exec();
        if (use_frequency_operator_coil_air_only)
        {
            vector_potential_curl_real_coil_contact_air_only.exec();
        }
        else
        {
            vector_potential_curl_real_coil_contact.exec();
        }
        curl_nu_b_real_coil_inner.exec();
        if (use_frequency_operator_coil_air_only)
        {
            curl_nu_b_real_coil_contact_air_only.exec();
        }
        else
        {
            curl_nu_b_real_coil_contact.exec();
        }
        vector_potential_curl_imag_coil_inner.exec();
        if (use_frequency_operator_coil_air_only)
        {
            vector_potential_curl_imag_coil_contact_air_only.exec();
        }
        else
        {
            vector_potential_curl_imag_coil_contact.exec();
        }
        curl_nu_b_imag_coil_inner.exec();
        if (use_frequency_operator_coil_air_only)
        {
            curl_nu_b_imag_coil_contact_air_only.exec();
        }
        else
        {
            curl_nu_b_imag_coil_contact.exec();
        }
        if (use_frequency_coil_component_hessian_inner)
        {
            for (size_t i = 0; i != total_coil_particles_count; ++i)
            {
                coil_curl_nu_b_real_direct_backup[i] = coil_curl_nu_b_real[i];
                coil_curl_nu_b_imag_direct_backup[i] = coil_curl_nu_b_imag[i];
            }
            if (run_operator_verification)
            {
                // Verification path: exercise the full component-hessian (inner + contact) rebuild.
                refresh_frequency_coil_magnetic_contact_by_component_hessian();
            }
            else
            {
                // Runtime path: use a stable hybrid
                // (component-hessian inner) + (direct contact contribution).
                for (size_t i = 0; i != total_coil_particles_count; ++i)
                {
                    coil_curl_nu_b_real[i] = ZeroData<Vecd>::value;
                    coil_curl_nu_b_imag[i] = ZeroData<Vecd>::value;
                }
                if (use_frequency_operator_coil_air_only)
                {
                    curl_nu_b_real_coil_contact_air_only.exec();
                    curl_nu_b_imag_coil_contact_air_only.exec();
                }
                else
                {
                    curl_nu_b_real_coil_contact.exec();
                    curl_nu_b_imag_coil_contact.exec();
                }
                for (size_t i = 0; i != total_coil_particles_count; ++i)
                {
                    coil_curl_nu_b_real_contact_backup[i] = coil_curl_nu_b_real[i];
                    coil_curl_nu_b_imag_contact_backup[i] = coil_curl_nu_b_imag[i];
                }
                refresh_frequency_coil_magnetic_inner_by_component_hessian();
                for (size_t i = 0; i != total_coil_particles_count; ++i)
                {
                    coil_curl_nu_b_real[i] += coil_curl_nu_b_real_contact_backup[i];
                    coil_curl_nu_b_imag[i] += coil_curl_nu_b_imag_contact_backup[i];
                }
                // Keep direct operator values in a thin inner shell where CK hessian
                // remains most sensitive to neighbor anisotropy near boundaries.
                size_t shell_blend_particles = 0;
                for (size_t i = 0; i != total_coil_particles_count; ++i)
                {
                    if (!is_in_inner_box_shell(coil_positions[i], coil_center, coil_halfsize,
                                               coil_interface_inner_shell_thickness))
                    {
                        continue;
                    }
                    coil_curl_nu_b_real[i] = coil_curl_nu_b_real_direct_backup[i];
                    coil_curl_nu_b_imag[i] = coil_curl_nu_b_imag_direct_backup[i];
                    shell_blend_particles++;
                }
                if (shell_blend_particles > 0)
                {
                    coil_component_hessian_runtime_shell_blend_activations++;
                    if (enable_coil_component_hessian_diagnostics &&
                        (coil_component_hessian_runtime_shell_blend_activations <= 3 ||
                         coil_component_hessian_runtime_shell_blend_activations % 100 == 0))
                    {
                        std::cout << "[team7-em] coil component-hessian shell blend kept direct curlNuB on "
                                  << shell_blend_particles
                                  << " inner-shell particles (activation="
                                  << coil_component_hessian_runtime_shell_blend_activations
                                  << ")."
                                  << std::endl;
                    }
                }

                if (enable_coil_component_hessian_core_quality_gate)
                {
                    size_t core_gate_candidates = 0;
                    size_t core_gate_rejected = 0;
                    for (size_t i = 0; i != total_coil_particles_count; ++i)
                    {
                        bool in_inner_shell = is_in_inner_box_shell(
                            coil_positions[i], coil_center, coil_halfsize,
                            coil_interface_inner_shell_thickness);
                        bool in_outer_shell = is_in_outer_box_shell(
                            coil_positions[i], coil_center, coil_halfsize,
                            coil_interface_air_shell_thickness);
                        if (in_inner_shell || in_outer_shell)
                        {
                            continue;
                        }

                        core_gate_candidates++;
                        Real diff_real =
                            (coil_curl_nu_b_real[i] - coil_curl_nu_b_real_direct_backup[i]).norm();
                        Real diff_imag =
                            (coil_curl_nu_b_imag[i] - coil_curl_nu_b_imag_direct_backup[i]).norm();
                        Real ref_real =
                            coil_curl_nu_b_real_direct_backup[i].norm() +
                            coil_component_hessian_core_gate_abs_ref;
                        Real ref_imag =
                            coil_curl_nu_b_imag_direct_backup[i].norm() +
                            coil_component_hessian_core_gate_abs_ref;
                        bool quality_ok =
                            diff_real <= coil_component_hessian_core_gate_rel_tol * ref_real &&
                            diff_imag <= coil_component_hessian_core_gate_rel_tol * ref_imag;
                        if (quality_ok)
                        {
                            continue;
                        }

                        coil_curl_nu_b_real[i] = coil_curl_nu_b_real_direct_backup[i];
                        coil_curl_nu_b_imag[i] = coil_curl_nu_b_imag_direct_backup[i];
                        core_gate_rejected++;
                    }

                    if (core_gate_rejected > 0)
                    {
                        coil_component_hessian_runtime_core_gate_activations++;
                        if (enable_coil_component_hessian_diagnostics &&
                            (coil_component_hessian_runtime_core_gate_activations <= 3 ||
                             coil_component_hessian_runtime_core_gate_activations % 100 == 0))
                        {
                            std::cout << "[team7-em] coil component-hessian core gate rejected "
                                      << core_gate_rejected << "/" << core_gate_candidates
                                      << " core particles (activation="
                                      << coil_component_hessian_runtime_core_gate_activations
                                      << ")."
                                      << std::endl;
                        }
                    }
                }
            }

            size_t fallback_particles = 0;
            size_t fallback_nonfinite_particles = 0;
            size_t fallback_outlier_particles = 0;
            size_t fallback_inner_shell_particles = 0;
            size_t fallback_core_particles = 0;
            size_t fallback_outer_shell_particles = 0;
            Real fallback_max_ratio_real = 0.0;
            Real fallback_max_ratio_imag = 0.0;
            for (size_t i = 0; i != total_coil_particles_count; ++i)
            {
                bool real_finite = std::isfinite(coil_curl_nu_b_real[i][0]) &&
                                   std::isfinite(coil_curl_nu_b_real[i][1]) &&
                                   std::isfinite(coil_curl_nu_b_real[i][2]);
                bool imag_finite = std::isfinite(coil_curl_nu_b_imag[i][0]) &&
                                   std::isfinite(coil_curl_nu_b_imag[i][1]) &&
                                   std::isfinite(coil_curl_nu_b_imag[i][2]);
                Real rebuilt_real_max_abs = SMAX(fabs(coil_curl_nu_b_real[i][0]),
                                                 SMAX(fabs(coil_curl_nu_b_real[i][1]),
                                                      fabs(coil_curl_nu_b_real[i][2])));
                Real rebuilt_imag_max_abs = SMAX(fabs(coil_curl_nu_b_imag[i][0]),
                                                 SMAX(fabs(coil_curl_nu_b_imag[i][1]),
                                                      fabs(coil_curl_nu_b_imag[i][2])));
                Real direct_real_max_abs = SMAX(fabs(coil_curl_nu_b_real_direct_backup[i][0]),
                                                SMAX(fabs(coil_curl_nu_b_real_direct_backup[i][1]),
                                                     fabs(coil_curl_nu_b_real_direct_backup[i][2])));
                Real direct_imag_max_abs = SMAX(fabs(coil_curl_nu_b_imag_direct_backup[i][0]),
                                                SMAX(fabs(coil_curl_nu_b_imag_direct_backup[i][1]),
                                                     fabs(coil_curl_nu_b_imag_direct_backup[i][2])));
                Real ratio_real = rebuilt_real_max_abs / (direct_real_max_abs + 1.0);
                Real ratio_imag = rebuilt_imag_max_abs / (direct_imag_max_abs + 1.0);
                bool is_extreme_outlier =
                    (ratio_real > coil_component_hessian_outlier_ratio_threshold &&
                     rebuilt_real_max_abs > coil_component_hessian_outlier_abs_threshold) ||
                    (ratio_imag > coil_component_hessian_outlier_ratio_threshold &&
                     rebuilt_imag_max_abs > coil_component_hessian_outlier_abs_threshold);
                if (real_finite && imag_finite && !is_extreme_outlier)
                {
                    continue;
                }
                coil_curl_nu_b_real[i] = coil_curl_nu_b_real_direct_backup[i];
                coil_curl_nu_b_imag[i] = coil_curl_nu_b_imag_direct_backup[i];
                fallback_particles++;
                fallback_max_ratio_real = SMAX(fallback_max_ratio_real, ratio_real);
                fallback_max_ratio_imag = SMAX(fallback_max_ratio_imag, ratio_imag);
                if (!real_finite || !imag_finite)
                {
                    fallback_nonfinite_particles++;
                }
                else
                {
                    fallback_outlier_particles++;
                }
                if (is_in_inner_box_shell(coil_positions[i], coil_center, coil_halfsize,
                                          coil_interface_inner_shell_thickness))
                {
                    fallback_inner_shell_particles++;
                }
                else if (is_in_outer_box_shell(coil_positions[i], coil_center, coil_halfsize,
                                               coil_interface_air_shell_thickness))
                {
                    fallback_outer_shell_particles++;
                }
                else
                {
                    fallback_core_particles++;
                }
            }
            if (fallback_particles > 0)
            {
                if (run_operator_verification)
                {
                    std::cout << "[team7-em] coil component-hessian fallback restored "
                              << fallback_particles << " particles to direct curlNuB."
                              << std::endl;
                }
                else if (!coil_component_hessian_runtime_warning_printed)
                {
                    std::cout << "[team7-em] coil component-hessian fallback active "
                              << "(first sweep restored " << fallback_particles
                              << " particles)." << std::endl;
                    coil_component_hessian_runtime_warning_printed = true;
                }
                if (enable_coil_component_hessian_diagnostics)
                {
                    std::cout << "[team7-em] coil fallback details: nonfinite="
                              << fallback_nonfinite_particles
                              << ", outlier=" << fallback_outlier_particles
                              << ", inner_shell=" << fallback_inner_shell_particles
                              << ", core=" << fallback_core_particles
                              << ", outer_shell=" << fallback_outer_shell_particles
                              << ", max_ratio_real=" << fallback_max_ratio_real
                              << ", max_ratio_imag=" << fallback_max_ratio_imag
                              << std::endl;
                }
            }
            if (!run_operator_verification &&
                fallback_particles * 100 > total_coil_particles_count)
            {
                for (size_t i = 0; i != total_coil_particles_count; ++i)
                {
                    coil_curl_nu_b_real[i] = coil_curl_nu_b_real_direct_backup[i];
                    coil_curl_nu_b_imag[i] = coil_curl_nu_b_imag_direct_backup[i];
                }
                coil_component_hessian_runtime_guard_activations++;
                if (coil_component_hessian_runtime_guard_activations <= 3 ||
                    coil_component_hessian_runtime_guard_activations % 100 == 0)
                {
                    std::cout << "[team7-em] coil component-hessian guard switched this sweep "
                              << "to direct curlNuB (fallback ratio > 1%, activation="
                              << coil_component_hessian_runtime_guard_activations << ")."
                              << std::endl;
                }
            }
        }
    };
    auto refresh_frequency_air_magnetic_terms = [&]()
    {
        air_curl_inner_avg_norm_last = 0.0;
        air_curl_contact_avg_norm_last = 0.0;
        air_curl_inner_max_norm_last = 0.0;
        air_curl_contact_max_norm_last = 0.0;

        if (!use_frequency_air_inner_term && !use_frequency_air_contact_term)
        {
            for (size_t i = 0; i != total_air_particles_count; ++i)
            {
                air_vector_potential_curl_real[i] = ZeroData<AngularVecd>::value;
                air_vector_potential_curl_imag[i] = ZeroData<AngularVecd>::value;
                air_curl_nu_b_real[i] = ZeroData<Vecd>::value;
                air_curl_nu_b_imag[i] = ZeroData<Vecd>::value;
                air_curl_nu_b_real_inner_backup[i] = ZeroData<Vecd>::value;
                air_curl_nu_b_imag_inner_backup[i] = ZeroData<Vecd>::value;
            }
            return;
        }

        if (use_frequency_air_inner_term)
        {
            vector_potential_curl_real_air_inner.exec();
            vector_potential_curl_imag_air_inner.exec();
        }
        else
        {
            for (size_t i = 0; i != total_air_particles_count; ++i)
            {
                air_vector_potential_curl_real[i] = ZeroData<AngularVecd>::value;
                air_vector_potential_curl_imag[i] = ZeroData<AngularVecd>::value;
            }
        }

        if (use_frequency_air_contact_term)
        {
            if (use_frequency_operator_coil_air_only)
            {
                vector_potential_curl_real_air_contact_coil_only.exec();
                vector_potential_curl_imag_air_contact_coil_only.exec();
            }
            else
            {
                vector_potential_curl_real_air_contact.exec();
                vector_potential_curl_imag_air_contact.exec();
            }
        }

        if (use_frequency_air_inner_term)
        {
            curl_nu_b_real_air_inner.exec();
            curl_nu_b_imag_air_inner.exec();
        }
        else
        {
            for (size_t i = 0; i != total_air_particles_count; ++i)
            {
                air_curl_nu_b_real[i] = ZeroData<Vecd>::value;
                air_curl_nu_b_imag[i] = ZeroData<Vecd>::value;
            }
        }

        for (size_t i = 0; i != total_air_particles_count; ++i)
        {
            air_curl_nu_b_real_inner_backup[i] = air_curl_nu_b_real[i];
            air_curl_nu_b_imag_inner_backup[i] = air_curl_nu_b_imag[i];
        }

        if (use_frequency_air_contact_term)
        {
            if (use_frequency_operator_coil_air_only)
            {
                curl_nu_b_real_air_contact_coil_only.exec();
                curl_nu_b_imag_air_contact_coil_only.exec();
            }
            else
            {
                curl_nu_b_real_air_contact.exec();
                curl_nu_b_imag_air_contact.exec();
            }
        }

        Real weighted_inner = 0.0;
        Real weighted_contact = 0.0;
        Real total_weight = 0.0;
        for (size_t i = 0; i != total_air_particles_count; ++i)
        {
            Vecd inner_real = air_curl_nu_b_real_inner_backup[i];
            Vecd inner_imag = air_curl_nu_b_imag_inner_backup[i];
            Vecd contact_real = air_curl_nu_b_real[i] - inner_real;
            Vecd contact_imag = air_curl_nu_b_imag[i] - inner_imag;

            Real inner_norm =
                sqrt(inner_real.squaredNorm() + inner_imag.squaredNorm());
            Real contact_norm =
                sqrt(contact_real.squaredNorm() + contact_imag.squaredNorm());
            Real weight = air_vol[i];

            weighted_inner += inner_norm * weight;
            weighted_contact += contact_norm * weight;
            total_weight += weight;
            air_curl_inner_max_norm_last =
                SMAX(air_curl_inner_max_norm_last, inner_norm);
            air_curl_contact_max_norm_last =
                SMAX(air_curl_contact_max_norm_last, contact_norm);
        }
        air_curl_inner_avg_norm_last =
            weighted_inner / (total_weight + TinyReal);
        air_curl_contact_avg_norm_last =
            weighted_contact / (total_weight + TinyReal);
    };
    auto refresh_frequency_plate_scalar_gradients = [&]()
    {
        electric_potential_gradient_real_plate_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_real_plate_contact.exec();
        }
        electric_potential_gradient_imag_plate_inner.exec();
        if (use_frequency_scalar_contact_coupling)
        {
            electric_potential_gradient_imag_plate_contact.exec();
        }
    };
    auto refresh_frequency_plate_magnetic_inner_by_component_hessian = [&]()
    {
        ensure_plate_component_hessian_ck_prepared();
        copy_plate_a_real_components.exec();
        plate_a_real_gradient_x_ck.exec();
        plate_a_real_hessian_x_ck.exec();
        plate_a_real_gradient_y_ck.exec();
        plate_a_real_hessian_y_ck.exec();
        plate_a_real_gradient_z_ck.exec();
        plate_a_real_hessian_z_ck.exec();
        reconstruct_plate_curl_nu_b_real_from_components.exec();

        copy_plate_a_imag_components.exec();
        plate_a_imag_gradient_x_ck.exec();
        plate_a_imag_hessian_x_ck.exec();
        plate_a_imag_gradient_y_ck.exec();
        plate_a_imag_hessian_y_ck.exec();
        plate_a_imag_gradient_z_ck.exec();
        plate_a_imag_hessian_z_ck.exec();
        reconstruct_plate_curl_nu_b_imag_from_components.exec();
    };
    auto refresh_frequency_plate_magnetic_contact_by_component_hessian = [&]()
    {
        ensure_plate_component_hessian_contact_ck_prepared();

        copy_plate_a_real_components.exec();
        copy_coil_a_real_components_for_plate.exec();
        copy_air_a_real_components_for_plate.exec();
        plate_a_real_gradient_x_contact_ck.exec();
        plate_a_real_hessian_x_contact_ck.exec();
        plate_a_real_gradient_y_contact_ck.exec();
        plate_a_real_hessian_y_contact_ck.exec();
        plate_a_real_gradient_z_contact_ck.exec();
        plate_a_real_hessian_z_contact_ck.exec();
        reconstruct_plate_curl_nu_b_real_from_components.exec();

        copy_plate_a_imag_components.exec();
        copy_coil_a_imag_components_for_plate.exec();
        copy_air_a_imag_components_for_plate.exec();
        plate_a_imag_gradient_x_contact_ck.exec();
        plate_a_imag_hessian_x_contact_ck.exec();
        plate_a_imag_gradient_y_contact_ck.exec();
        plate_a_imag_hessian_y_contact_ck.exec();
        plate_a_imag_gradient_z_contact_ck.exec();
        plate_a_imag_hessian_z_contact_ck.exec();
        reconstruct_plate_curl_nu_b_imag_from_components.exec();
    };
    auto refresh_frequency_plate_magnetic_terms = [&]()
    {
        vector_potential_curl_real_plate_inner.exec();
        vector_potential_curl_real_plate_contact.exec();
        if (use_frequency_plate_component_hessian_inner)
        {
            refresh_frequency_plate_magnetic_contact_by_component_hessian();
        }
        else
        {
            curl_nu_b_real_plate_inner.exec();
            curl_nu_b_real_plate_contact.exec();
        }
        vector_potential_curl_imag_plate_inner.exec();
        vector_potential_curl_imag_plate_contact.exec();
        if (!use_frequency_plate_component_hessian_inner)
        {
            curl_nu_b_imag_plate_inner.exec();
            curl_nu_b_imag_plate_contact.exec();
        }
    };
    auto apply_frequency_plate_a_constraints = [&]()
    {
        if (constrain_a_reference_plate_enabled)
        {
            constrain_a_real_reference_plate.exec();
            constrain_a_imag_reference_plate.exec();
        }
    };
    auto apply_frequency_plate_phi_constraints = [&]()
    {
        constrain_phi_real_reference_plate.exec();
        if (constrain_phi_boundary_plate_enabled)
        {
            constrain_phi_real_boundary_plate.exec();
        }
        constrain_phi_imag_reference_plate.exec();
        if (constrain_phi_boundary_plate_enabled)
        {
            constrain_phi_imag_boundary_plate.exec();
        }
    };
    auto apply_frequency_phi_constraints = [&]()
    {
        apply_frequency_plate_phi_constraints();

        if (use_frequency_coil_scalar_potential)
        {
            constrain_phi_real_reference_coil.exec();
            if (constrain_phi_boundary_coil_enabled)
            {
                constrain_phi_real_boundary_coil.exec();
            }
            constrain_phi_imag_reference_coil.exec();
            if (constrain_phi_boundary_coil_enabled)
            {
                constrain_phi_imag_boundary_coil.exec();
            }
        }

        if (use_frequency_air_scalar_potential)
        {
            constrain_phi_real_reference_air.exec();
            if (constrain_phi_boundary_air_enabled)
            {
                constrain_phi_real_boundary_air.exec();
            }
            if (constrain_phi_all_air_enabled)
            {
                constrain_phi_real_all_air.exec();
            }
            constrain_phi_imag_reference_air.exec();
            if (constrain_phi_boundary_air_enabled)
            {
                constrain_phi_imag_boundary_air.exec();
            }
            if (constrain_phi_all_air_enabled)
            {
                constrain_phi_imag_all_air.exec();
            }
        }
    };
    auto apply_frequency_a_constraints = [&]()
    {
        if (constrain_a_reference_coil_enabled)
        {
            constrain_a_real_reference_coil.exec();
            constrain_a_imag_reference_coil.exec();
        }
        if (constrain_a_reference_air_enabled)
        {
            constrain_a_real_reference_air.exec();
            constrain_a_imag_reference_air.exec();
        }
        if (constrain_a_boundary_air_enabled)
        {
            apply_frequency_air_a_boundary_constraint();
        }
        apply_frequency_plate_a_constraints();
    };
    auto execute_frequency_air_magnetic_sweep = [&](Real dt_air_step)
    {
        if (!(dt_air_step > TinyReal))
        {
            return;
        }
        size_t air_substeps = SMAX(static_cast<size_t>(1), frequency_air_block_sweeps);
        Real dt_air_substep = dt_air_step / static_cast<Real>(air_substeps);
        for (size_t air_substep = 0; air_substep < air_substeps; ++air_substep)
        {
            refresh_frequency_air_scalar_gradients();
            refresh_frequency_air_magnetic_terms();
            if (use_frequency_operator_coil_air_only)
            {
                vector_potential_equation_real_air_coil_only.exec(dt_air_substep);
                vector_potential_equation_imag_air_coil_only.exec(dt_air_substep);
            }
            else
            {
                vector_potential_equation_real_air.exec(dt_air_substep);
                vector_potential_equation_imag_air.exec(dt_air_substep);
            }
            if (constrain_a_reference_air_enabled)
            {
                constrain_a_real_reference_air.exec();
                constrain_a_imag_reference_air.exec();
            }
            if (constrain_a_boundary_air_enabled)
            {
                apply_frequency_air_a_boundary_constraint();
            }
        }
    };
    auto backup_frequency_plate_local_state = [&]()
    {
        if (!use_frequency_aphi ||
            plate_vector_potential_real_backup.empty() ||
            plate_vector_potential_imag_backup.empty() ||
            plate_electric_potential_real_backup.empty() ||
            plate_electric_potential_imag_backup.empty())
        {
            return;
        }
        std::copy(plate_vector_potential_real,
                  plate_vector_potential_real + total_plate_particles_count,
                  plate_vector_potential_real_backup.begin());
        std::copy(plate_vector_potential_imag,
                  plate_vector_potential_imag + total_plate_particles_count,
                  plate_vector_potential_imag_backup.begin());
        std::copy(plate_electric_potential_real,
                  plate_electric_potential_real + total_plate_particles_count,
                  plate_electric_potential_real_backup.begin());
        std::copy(plate_electric_potential_imag,
                  plate_electric_potential_imag + total_plate_particles_count,
                  plate_electric_potential_imag_backup.begin());
        std::copy(plate_vector_potential_change_rate_real,
                  plate_vector_potential_change_rate_real + total_plate_particles_count,
                  plate_vector_potential_change_rate_real_backup.begin());
        std::copy(plate_vector_potential_change_rate_imag,
                  plate_vector_potential_change_rate_imag + total_plate_particles_count,
                  plate_vector_potential_change_rate_imag_backup.begin());
    };
    auto restore_frequency_plate_local_state = [&]()
    {
        if (!use_frequency_aphi ||
            plate_vector_potential_real_backup.empty() ||
            plate_vector_potential_imag_backup.empty() ||
            plate_electric_potential_real_backup.empty() ||
            plate_electric_potential_imag_backup.empty())
        {
            return;
        }
        std::copy(plate_vector_potential_real_backup.begin(),
                  plate_vector_potential_real_backup.end(),
                  plate_vector_potential_real);
        std::copy(plate_vector_potential_imag_backup.begin(),
                  plate_vector_potential_imag_backup.end(),
                  plate_vector_potential_imag);
        std::copy(plate_electric_potential_real_backup.begin(),
                  plate_electric_potential_real_backup.end(),
                  plate_electric_potential_real);
        std::copy(plate_electric_potential_imag_backup.begin(),
                  plate_electric_potential_imag_backup.end(),
                  plate_electric_potential_imag);
        std::copy(plate_vector_potential_change_rate_real_backup.begin(),
                  plate_vector_potential_change_rate_real_backup.end(),
                  plate_vector_potential_change_rate_real);
        std::copy(plate_vector_potential_change_rate_imag_backup.begin(),
                  plate_vector_potential_change_rate_imag_backup.end(),
                  plate_vector_potential_change_rate_imag);
        apply_frequency_plate_phi_constraints();
        apply_frequency_plate_a_constraints();
    };
    auto zero_frequency_plate_change_rate = [&]()
    {
        if (!use_frequency_aphi)
        {
            return;
        }
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            plate_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
            plate_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
        }
    };
    auto backup_frequency_em_state = [&]()
    {
        if (!use_frequency_em_state_backup)
        {
            return;
        }
        std::copy(plate_vector_potential_real,
                  plate_vector_potential_real + total_plate_particles_count,
                  frequency_em_state_backup.plate_a_real.begin());
        std::copy(plate_vector_potential_imag,
                  plate_vector_potential_imag + total_plate_particles_count,
                  frequency_em_state_backup.plate_a_imag.begin());
        std::copy(plate_electric_potential_real,
                  plate_electric_potential_real + total_plate_particles_count,
                  frequency_em_state_backup.plate_phi_real.begin());
        std::copy(plate_electric_potential_imag,
                  plate_electric_potential_imag + total_plate_particles_count,
                  frequency_em_state_backup.plate_phi_imag.begin());

        size_t total_coil_particles_count = coil_particles.TotalRealParticles();
        std::copy(coil_vector_potential_real,
                  coil_vector_potential_real + total_coil_particles_count,
                  frequency_em_state_backup.coil_a_real.begin());
        std::copy(coil_vector_potential_imag,
                  coil_vector_potential_imag + total_coil_particles_count,
                  frequency_em_state_backup.coil_a_imag.begin());
        if (use_frequency_coil_scalar_potential)
        {
            std::copy(coil_electric_potential_real,
                      coil_electric_potential_real + total_coil_particles_count,
                      frequency_em_state_backup.coil_phi_real.begin());
            std::copy(coil_electric_potential_imag,
                      coil_electric_potential_imag + total_coil_particles_count,
                      frequency_em_state_backup.coil_phi_imag.begin());
        }

        std::copy(air_vector_potential_real,
                  air_vector_potential_real + total_air_particles_count,
                  frequency_em_state_backup.air_a_real.begin());
        std::copy(air_vector_potential_imag,
                  air_vector_potential_imag + total_air_particles_count,
                  frequency_em_state_backup.air_a_imag.begin());
        if (use_frequency_air_scalar_potential)
        {
            std::copy(air_electric_potential_real,
                      air_electric_potential_real + total_air_particles_count,
                      frequency_em_state_backup.air_phi_real.begin());
            std::copy(air_electric_potential_imag,
                      air_electric_potential_imag + total_air_particles_count,
                      frequency_em_state_backup.air_phi_imag.begin());
        }
    };
    auto restore_frequency_em_state = [&]()
    {
        if (!use_frequency_em_state_backup)
        {
            return;
        }
        std::copy(frequency_em_state_backup.plate_a_real.begin(),
                  frequency_em_state_backup.plate_a_real.end(),
                  plate_vector_potential_real);
        std::copy(frequency_em_state_backup.plate_a_imag.begin(),
                  frequency_em_state_backup.plate_a_imag.end(),
                  plate_vector_potential_imag);
        std::copy(frequency_em_state_backup.plate_phi_real.begin(),
                  frequency_em_state_backup.plate_phi_real.end(),
                  plate_electric_potential_real);
        std::copy(frequency_em_state_backup.plate_phi_imag.begin(),
                  frequency_em_state_backup.plate_phi_imag.end(),
                  plate_electric_potential_imag);

        size_t total_coil_particles_count = coil_particles.TotalRealParticles();
        std::copy(frequency_em_state_backup.coil_a_real.begin(),
                  frequency_em_state_backup.coil_a_real.end(),
                  coil_vector_potential_real);
        std::copy(frequency_em_state_backup.coil_a_imag.begin(),
                  frequency_em_state_backup.coil_a_imag.end(),
                  coil_vector_potential_imag);
        if (use_frequency_coil_scalar_potential)
        {
            std::copy(frequency_em_state_backup.coil_phi_real.begin(),
                      frequency_em_state_backup.coil_phi_real.end(),
                      coil_electric_potential_real);
            std::copy(frequency_em_state_backup.coil_phi_imag.begin(),
                      frequency_em_state_backup.coil_phi_imag.end(),
                      coil_electric_potential_imag);
        }

        std::copy(frequency_em_state_backup.air_a_real.begin(),
                  frequency_em_state_backup.air_a_real.end(),
                  air_vector_potential_real);
        std::copy(frequency_em_state_backup.air_a_imag.begin(),
                  frequency_em_state_backup.air_a_imag.end(),
                  air_vector_potential_imag);
        if (use_frequency_air_scalar_potential)
        {
            std::copy(frequency_em_state_backup.air_phi_real.begin(),
                      frequency_em_state_backup.air_phi_real.end(),
                      air_electric_potential_real);
            std::copy(frequency_em_state_backup.air_phi_imag.begin(),
                      frequency_em_state_backup.air_phi_imag.end(),
                      air_electric_potential_imag);
        }

        apply_frequency_phi_constraints();
        apply_frequency_a_constraints();
        refresh_frequency_coil_scalar_gradients();
        refresh_frequency_coil_magnetic_terms();
        refresh_frequency_air_scalar_gradients();
        refresh_frequency_air_magnetic_terms();
        refresh_frequency_plate_scalar_gradients();
        refresh_frequency_plate_magnetic_terms();
    };
    auto zero_frequency_em_change_rate = [&]()
    {
        if (!use_frequency_aphi)
        {
            return;
        }
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            plate_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
            plate_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
        }
        for (size_t i = 0; i != total_coil_particles_count; ++i)
        {
            coil_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
            coil_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
        }
        for (size_t i = 0; i != total_air_particles_count; ++i)
        {
            air_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
            air_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
        }
    };
    auto evaluate_plate_acceptance_diagnostics = [&]() -> FrequencyPlateAcceptanceDiagnostics
    {
        FrequencyPlateAcceptanceDiagnostics diagnostics;
        if (!use_frequency_aphi)
        {
            return diagnostics;
        }

        Real global_reference_scale = evaluate_em_residual_reference_scale();
        auto accumulate_conductor_average_ratio =
            [&](Real &ratio_accumulator,
                Vecd *positions,
                Vecd *vector_potential_real,
                Vecd *vector_potential_imag,
                Vecd *electric_potential_gradient_real,
                Vecd *electric_potential_gradient_imag,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *electrical_conductivity,
                Real *volumetric_measure,
                size_t total_particles,
                const auto &selector)
        {
            Real weighted_sum = 0.0;
            Real total_weight = 0.0;
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_conductor_residual_terms(
                    vector_potential_real[i], vector_potential_imag[i],
                    electric_potential_gradient_real[i], electric_potential_gradient_imag[i],
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    electrical_conductivity[i], residual_norm, reference_norm);
                Real weight = volumetric_measure[i];
                weighted_sum +=
                    (residual_norm / (global_reference_scale + TinyReal)) * weight;
                total_weight += weight;
            }
            ratio_accumulator = weighted_sum / (total_weight + TinyReal);
        };
        auto accumulate_magnetic_average_ratio =
            [&](Real &ratio_accumulator,
                Vecd *positions,
                Vecd *curl_nu_b_real,
                Vecd *curl_nu_b_imag,
                Vecd *source_current_density_real,
                Vecd *source_current_density_imag,
                Real *volumetric_measure,
                size_t total_particles,
                const auto &selector)
        {
            Real weighted_sum = 0.0;
            Real total_weight = 0.0;
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector(positions[i]))
                {
                    continue;
                }
                Vecd source_real = source_current_density_real != nullptr
                                       ? source_current_density_real[i]
                                       : ZeroData<Vecd>::value;
                Vecd source_imag = source_current_density_imag != nullptr
                                       ? source_current_density_imag[i]
                                       : ZeroData<Vecd>::value;
                Real residual_norm = 0.0;
                Real reference_norm = 0.0;
                compute_frequency_magnetic_residual_terms(
                    curl_nu_b_real[i], curl_nu_b_imag[i], source_real, source_imag,
                    residual_norm, reference_norm);
                Real weight = volumetric_measure[i];
                weighted_sum +=
                    (residual_norm / (global_reference_scale + TinyReal)) * weight;
                total_weight += weight;
            }
            ratio_accumulator = weighted_sum / (total_weight + TinyReal);
        };

        auto select_all_plate = [](const Vecd &) -> bool
        { return true; };
        auto select_plate_inner_shell = [&](const Vecd &position) -> bool
        {
            return is_in_inner_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_inner_shell_thickness);
        };
        auto select_plate_air_shell = [&](const Vecd &position) -> bool
        {
            return is_in_outer_box_shell(position, plate_center, plate_halfsize,
                                         plate_interface_air_shell_thickness);
        };
        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        accumulate_conductor_average_ratio(diagnostics.body_avg_residual_ratio,
                                           plate_positions,
                                           plate_vector_potential_real,
                                           plate_vector_potential_imag,
                                           plate_electric_potential_gradient_real,
                                           plate_electric_potential_gradient_imag,
                                           plate_curl_nu_b_real,
                                           plate_curl_nu_b_imag,
                                           nullptr,
                                           nullptr,
                                           plate_sigma,
                                           plate_vol,
                                           total_plate_particles_count,
                                           select_all_plate);
        accumulate_conductor_average_ratio(diagnostics.inner_shell_avg_residual_ratio,
                                           plate_positions,
                                           plate_vector_potential_real,
                                           plate_vector_potential_imag,
                                           plate_electric_potential_gradient_real,
                                           plate_electric_potential_gradient_imag,
                                           plate_curl_nu_b_real,
                                           plate_curl_nu_b_imag,
                                           nullptr,
                                           nullptr,
                                           plate_sigma,
                                           plate_vol,
                                           total_plate_particles_count,
                                           select_plate_inner_shell);
        accumulate_magnetic_average_ratio(diagnostics.air_shell_avg_residual_ratio,
                                          air_positions,
                                          air_curl_nu_b_real,
                                          air_curl_nu_b_imag,
                                          nullptr,
                                          nullptr,
                                          air_vol,
                                          total_air_particles_count,
                                          select_plate_air_shell);
        return diagnostics;
    };
    auto evaluate_plate_acceptance_metric =
        [&](const FrequencyPlateAcceptanceDiagnostics &diagnostics) -> Real
    {
        return SMAX(diagnostics.body_avg_residual_ratio,
                    SMAX(diagnostics.inner_shell_avg_residual_ratio,
                         diagnostics.air_shell_avg_residual_ratio));
    };
    auto evaluate_frequency_em_backtracking_diagnostics =
        [&]() -> FrequencyEmAcceptanceDiagnostics
    {
        FrequencyEmAcceptanceDiagnostics diagnostics;
        if (!use_frequency_aphi)
        {
            return diagnostics;
        }

        EmResidualBreakdown residual_breakdown = evaluate_em_residual_breakdown();
        if (!use_frequency_operator_coil_air_only)
        {
            FrequencyPlateAcceptanceDiagnostics plate_diagnostics =
                evaluate_plate_acceptance_diagnostics();
            diagnostics.plate_shell_metric =
                evaluate_plate_acceptance_metric(plate_diagnostics);
            diagnostics.conductor_metric =
                SMAX(diagnostics.plate_shell_metric,
                     SMAX(residual_breakdown.plate, residual_breakdown.coil));
        }
        else
        {
            diagnostics.plate_shell_metric = 0.0;
            diagnostics.conductor_metric = residual_breakdown.coil;
        }
        diagnostics.air_metric = residual_breakdown.air;
        diagnostics.weighted_metric =
            diagnostics.conductor_metric +
            frequency_global_backtracking_air_penalty * diagnostics.air_metric;
        return diagnostics;
    };
    auto evaluate_frequency_em_hard_guard_diagnostics =
        [&]() -> FrequencyEmHardGuardDiagnostics
    {
        FrequencyEmHardGuardDiagnostics diagnostics;
        if (!use_frequency_aphi)
        {
            return diagnostics;
        }
        FrequencyEmAcceptanceDiagnostics acceptance_diagnostics =
            evaluate_frequency_em_backtracking_diagnostics();
        diagnostics.weighted_metric = acceptance_diagnostics.weighted_metric;
        Real b_air_above_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_above_plate);
        Real b_air_below_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_below_plate);
        diagnostics.probe_b_metric = SMAX(b_air_above_mag, b_air_below_mag);
        std::pair<Vec3d, Vec3d> b_ref_air_above = compute_biot_savart_magnetic_flux_density(
            probe_b_air_above_plate.sampled_position);
        std::pair<Vec3d, Vec3d> b_ref_air_below = compute_biot_savart_magnetic_flux_density(
            probe_b_air_below_plate.sampled_position);
        Real b_ref_air_above_mag =
            complex_vector_magnitude(b_ref_air_above.first, b_ref_air_above.second);
        Real b_ref_air_below_mag =
            complex_vector_magnitude(b_ref_air_below.first, b_ref_air_below.second);
        auto compute_probe_reference_ratio = [&](Real simulated_mag, Real reference_mag) -> Real
        {
            if (!std::isfinite(simulated_mag) || !std::isfinite(reference_mag) ||
                reference_mag <= TinyReal)
            {
                return std::numeric_limits<Real>::infinity();
            }
            return simulated_mag / (reference_mag + TinyReal);
        };
        Real probe_above_ref_ratio =
            compute_probe_reference_ratio(b_air_above_mag, b_ref_air_above_mag);
        Real probe_below_ref_ratio =
            compute_probe_reference_ratio(b_air_below_mag, b_ref_air_below_mag);
        diagnostics.probe_ref_ratio_metric =
            SMAX(probe_above_ref_ratio, probe_below_ref_ratio);
        return diagnostics;
    };
    auto report_coil_nan = [&](const std::string &tag)
    {
        size_t phi_nan = 0;
        size_t a_nan = 0;
        size_t a_dot_nan = 0;
        size_t a_change_rate_nan = 0;
        size_t total_real_particles = coil_particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            if (!std::isfinite(coil_electric_potential[i]))
            {
                phi_nan++;
            }
            for (int d = 0; d != Dimensions; ++d)
            {
                if (!std::isfinite(coil_vector_potential[i][d]))
                {
                    a_nan++;
                }
                if (!std::isfinite(coil_vector_potential_dt[i][d]))
                {
                    a_dot_nan++;
                }
                if (!std::isfinite(coil_vector_potential_change_rate[i][d]))
                {
                    a_change_rate_nan++;
                }
            }
        }
        if (phi_nan + a_nan + a_dot_nan + a_change_rate_nan > 0)
        {
            std::cout << "[debug] coil nan detected at " << tag
                      << ", phi_nan=" << phi_nan
                      << ", A_nan_components=" << a_nan
                      << ", A_dot_nan_components=" << a_dot_nan
                      << ", A_rate_nan_components=" << a_change_rate_nan
                      << std::endl;
        }
    };
    auto report_coil_rhs = [&](const std::string &tag)
    {
        Real max_source = 0.0;
        Real max_grad_phi = 0.0;
        Real max_curl_nu_b = 0.0;
        Real max_a_rate = 0.0;
        size_t total_real_particles = coil_particles.TotalRealParticles();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            if (use_frequency_aphi)
            {
                Real source_mag = static_cast<Real>(
                    sqrt(coil_source_current_density_real[i].squaredNorm() +
                         coil_source_current_density_imag[i].squaredNorm()));
                max_source = SMAX(max_source, source_mag);
                max_grad_phi = SMAX(max_grad_phi, Real(0));
                max_curl_nu_b = SMAX(max_curl_nu_b, Real(0));
                max_a_rate = SMAX(max_a_rate, Real(0));
            }
            else
            {
                max_source = SMAX(max_source, coil_source_current_density[i].norm());
                max_grad_phi = SMAX(max_grad_phi, coil_electric_potential_gradient[i].norm());
                max_curl_nu_b = SMAX(max_curl_nu_b, coil_curl_nu_b[i].norm());
                max_a_rate = SMAX(max_a_rate, coil_vector_potential_change_rate[i].norm());
            }
        }
        std::cout << "[debug] " << tag
                  << ", max|Js|=" << max_source
                  << ", max|gradPhi|=" << max_grad_phi
                  << ", max|curlNuB|=" << max_curl_nu_b
                  << ", max|A_rate|=" << max_a_rate
                  << std::endl;
    };
    auto report_plate_air_nan = [&](const std::string &tag)
    {
        size_t plate_phi_nan = 0;
        size_t plate_a_nan = 0;
        size_t plate_a_dot_nan = 0;
        size_t plate_grad_phi_nan = 0;
        size_t plate_e_nan = 0;
        size_t plate_j_nan = 0;
        size_t plate_q_nan = 0;
        size_t total_plate_particles = plate_particles.TotalRealParticles();
        for (size_t i = 0; i != total_plate_particles; ++i)
        {
            if (!std::isfinite(plate_electric_potential[i]))
            {
                plate_phi_nan++;
            }
            if (!std::isfinite(plate_joule_heat[i]))
            {
                plate_q_nan++;
            }
            for (int d = 0; d != Dimensions; ++d)
            {
                if (!std::isfinite(plate_vector_potential[i][d]))
                {
                    plate_a_nan++;
                }
                if (!std::isfinite(plate_vector_potential_dt[i][d]))
                {
                    plate_a_dot_nan++;
                }
                if (!std::isfinite(plate_electric_potential_gradient[i][d]))
                {
                    plate_grad_phi_nan++;
                }
                if (!std::isfinite(plate_electric_field[i][d]))
                {
                    plate_e_nan++;
                }
                if (!std::isfinite(plate_current_density[i][d]))
                {
                    plate_j_nan++;
                }
            }
        }

        size_t air_phi_nan = 0;
        size_t air_a_nan = 0;
        size_t air_a_dot_nan = 0;
        size_t air_grad_phi_nan = 0;
        size_t total_air_particles = air_particles.TotalRealParticles();
        for (size_t i = 0; i != total_air_particles; ++i)
        {
            if (!std::isfinite(air_electric_potential[i]))
            {
                air_phi_nan++;
            }
            for (int d = 0; d != Dimensions; ++d)
            {
                if (!std::isfinite(air_vector_potential[i][d]))
                {
                    air_a_nan++;
                }
                if (!std::isfinite(air_vector_potential_dt[i][d]))
                {
                    air_a_dot_nan++;
                }
                if (!std::isfinite(air_electric_potential_gradient[i][d]))
                {
                    air_grad_phi_nan++;
                }
            }
        }

        size_t total_nan_count =
            plate_phi_nan + plate_a_nan + plate_a_dot_nan + plate_grad_phi_nan +
            plate_e_nan + plate_j_nan + plate_q_nan +
            air_phi_nan + air_a_nan + air_a_dot_nan + air_grad_phi_nan;
        if (total_nan_count > 0)
        {
            std::cout << "[debug] plate/air nan detected at " << tag
                      << ", plate_phi_nan=" << plate_phi_nan
                      << ", plate_A_nan_components=" << plate_a_nan
                      << ", plate_A_dot_nan_components=" << plate_a_dot_nan
                      << ", plate_gradPhi_nan_components=" << plate_grad_phi_nan
                      << ", plate_E_nan_components=" << plate_e_nan
                      << ", plate_J_nan_components=" << plate_j_nan
                      << ", plate_Q_nan=" << plate_q_nan
                      << ", air_phi_nan=" << air_phi_nan
                      << ", air_A_nan_components=" << air_a_nan
                      << ", air_A_dot_nan_components=" << air_a_dot_nan
                      << ", air_gradPhi_nan_components=" << air_grad_phi_nan
                      << std::endl;
        }
    };
    auto run_frequency_operator_verification = [&]() -> int
    {
        if (!run_operator_verification)
        {
            return -1;
        }
        if (!use_frequency_aphi)
        {
            std::cerr << "[team7-verify] TEAM7_RUN_OPERATOR_VERIFICATION requires TEAM7_USE_FREQ_APHI=1"
                      << std::endl;
            return 1;
        }

        Real *plate_sigma = plate_particles.getVariableDataByName<Real>("ElectricalConductivity");
        Real *plate_source_real = plate_particles.getVariableDataByName<Real>("ElectricPotentialSourceReal");
        Real *plate_source_imag = plate_particles.getVariableDataByName<Real>("ElectricPotentialSourceImag");
        Real *plate_temperature_change_rate_by_joule =
            plate_particles.getVariableDataByName<Real>("TemperatureChangeRateByJoule");
        std::ofstream verification_file(
            io_environment.OutputFolder() + "/team7_operator_verification.csv",
            std::ios::out | std::ios::trunc);
        verification_file << std::setprecision(12);
        verification_file
            << "test,quantity,checked_particles,invalid_particles,mean_abs_error,max_abs_error,mean_rel_error,mean_projection_ratio,mean_alignment\n";

        Real verification_margin =
            SMIN(operator_verification_margin_dp * dp_plate,
                 static_cast<Real>(0.45) * plate_halfsize.minCoeff());
        Vecd plate_verification_halfsize =
            plate_halfsize - Vecd::Ones() * verification_margin;
        Real coil_verification_margin =
            SMIN(operator_verification_margin_dp * dp_coil,
                 static_cast<Real>(0.45) * coil_halfsize.minCoeff());
        Vecd coil_verification_halfsize =
            coil_halfsize - Vecd::Ones() * coil_verification_margin;
        auto select_plate_core = [&](const Vecd &position) -> bool
        {
            if (plate_verification_halfsize.minCoeff() <= TinyReal)
            {
                return true;
            }
            Vecd local_abs = (position - plate_center).cwiseAbs();
            return (local_abs.array() <= plate_verification_halfsize.array()).all();
        };
        auto select_plate_contact_band = [&](const Vecd &position) -> bool
        {
            if (plate_verification_halfsize.minCoeff() <= TinyReal)
            {
                return true;
            }
            Vecd local_abs = (position - plate_center).cwiseAbs();
            return (local_abs.array() > plate_verification_halfsize.array()).any();
        };
        auto select_coil_core = [&](const Vecd &position) -> bool
        {
            if (coil_verification_halfsize.minCoeff() <= TinyReal)
            {
                return true;
            }
            Vecd local_abs = (position - coil_center).cwiseAbs();
            return (local_abs.array() <= coil_verification_halfsize.array()).all();
        };
        auto select_coil_contact_band = [&](const Vecd &position) -> bool
        {
            if (coil_verification_halfsize.minCoeff() <= TinyReal)
            {
                return true;
            }
            Vecd local_abs = (position - coil_center).cwiseAbs();
            return (local_abs.array() > coil_verification_halfsize.array()).any();
        };
        auto plate_local_position_si = [&](const Vecd &position) -> Vecd
        {
            return (position - plate_center) * geom_length_to_m;
        };
        auto coil_local_position_si = [&](const Vecd &position) -> Vecd
        {
            return (position - coil_center) * geom_length_to_m;
        };
        auto is_finite_scalar = [&](Real value) -> bool
        {
            return std::isfinite(value);
        };
        auto is_finite_vector = [&](const Vecd &value) -> bool
        {
            for (int axis = 0; axis != static_cast<int>(value.size()); ++axis)
            {
                if (!std::isfinite(value[axis]))
                {
                    return false;
                }
            }
            return true;
        };
        auto reset_plate_frequency_fields = [&]()
        {
            for (size_t i = 0; i != total_plate_particles_count; ++i)
            {
                plate_vector_potential_real[i] = ZeroData<Vecd>::value;
                plate_vector_potential_imag[i] = ZeroData<Vecd>::value;
                plate_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
                plate_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
                plate_curl_nu_b_real[i] = ZeroData<Vecd>::value;
                plate_curl_nu_b_imag[i] = ZeroData<Vecd>::value;
                plate_electric_potential_gradient_real[i] = ZeroData<Vecd>::value;
                plate_electric_potential_gradient_imag[i] = ZeroData<Vecd>::value;
                plate_electric_field_real[i] = ZeroData<Vecd>::value;
                plate_electric_field_imag[i] = ZeroData<Vecd>::value;
                plate_current_density_real[i] = ZeroData<Vecd>::value;
                plate_current_density_imag[i] = ZeroData<Vecd>::value;
                plate_vector_potential_curl_real[i] = ZeroData<AngularVecd>::value;
                plate_vector_potential_curl_imag[i] = ZeroData<AngularVecd>::value;
                plate_electric_potential_real[i] = 0.0;
                plate_electric_potential_imag[i] = 0.0;
                plate_source_real[i] = 0.0;
                plate_source_imag[i] = 0.0;
                plate_joule_heat[i] = 0.0;
                plate_temperature_change_rate_by_joule[i] = 0.0;
            }
        };
        auto record_scalar_verification =
            [&](const std::string &test_name,
                const std::string &quantity_name,
                Real *sample_values,
                const auto &reference_fn)
        {
            OperatorVerificationRow row;
            row.test_name = test_name;
            row.quantity_name = quantity_name;
            bool has_first_invalid = false;
            Real first_invalid_sample = 0.0;
            Real first_invalid_reference = 0.0;
            for (size_t i = 0; i != total_plate_particles_count; ++i)
            {
                if (!select_plate_core(plate_positions[i]))
                {
                    continue;
                }
                Real reference_value = reference_fn(i);
                Real sample_value = sample_values[i];
                if (!is_finite_scalar(reference_value) || !is_finite_scalar(sample_value))
                {
                    if (!has_first_invalid)
                    {
                        has_first_invalid = true;
                        first_invalid_sample = sample_value;
                        first_invalid_reference = reference_value;
                    }
                    row.invalid_particles++;
                    continue;
                }
                Real abs_error = fabs(sample_value - reference_value);
                Real rel_error = abs_error / (fabs(reference_value) + TinyReal);
                row.checked_particles++;
                row.mean_abs_error += abs_error;
                row.max_abs_error = SMAX(row.max_abs_error, abs_error);
                row.mean_rel_error += rel_error;
            }
            if (row.checked_particles > 0)
            {
                Real inv_count = 1.0 / static_cast<Real>(row.checked_particles);
                row.mean_abs_error *= inv_count;
                row.mean_rel_error *= inv_count;
            }
            verification_file << row.test_name << ","
                              << row.quantity_name << ","
                              << row.checked_particles << ","
                              << row.invalid_particles << ","
                              << row.mean_abs_error << ","
                              << row.max_abs_error << ","
                              << row.mean_rel_error << ","
                              << row.mean_projection_ratio << ","
                              << row.mean_alignment << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-verify] " << row.test_name
                      << " / " << row.quantity_name
                      << ": count=" << row.checked_particles
                      << ", invalid=" << row.invalid_particles
                      << ", mean_abs=" << row.mean_abs_error
                      << ", max_abs=" << row.max_abs_error
                      << ", mean_rel=" << row.mean_rel_error
                      << ", proj_ratio=" << row.mean_projection_ratio
                      << ", alignment=" << row.mean_alignment
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (has_first_invalid)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-verify] first invalid scalar sample/ref for "
                          << row.test_name << " / " << row.quantity_name
                          << ": sample=" << first_invalid_sample
                          << ", reference=" << first_invalid_reference
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
        };
        auto record_scalar_function_verification =
            [&](const std::string &test_name,
                const std::string &quantity_name,
                const auto &sample_fn,
                const auto &reference_fn)
        {
            OperatorVerificationRow row;
            row.test_name = test_name;
            row.quantity_name = quantity_name;
            bool has_first_invalid = false;
            Real first_invalid_sample = 0.0;
            Real first_invalid_reference = 0.0;
            for (size_t i = 0; i != total_plate_particles_count; ++i)
            {
                if (!select_plate_core(plate_positions[i]))
                {
                    continue;
                }
                Real reference_value = reference_fn(i);
                Real sample_value = sample_fn(i);
                if (!is_finite_scalar(reference_value) || !is_finite_scalar(sample_value))
                {
                    if (!has_first_invalid)
                    {
                        has_first_invalid = true;
                        first_invalid_sample = sample_value;
                        first_invalid_reference = reference_value;
                    }
                    row.invalid_particles++;
                    continue;
                }
                Real abs_error = fabs(sample_value - reference_value);
                Real rel_error = abs_error / (fabs(reference_value) + TinyReal);
                row.checked_particles++;
                row.mean_abs_error += abs_error;
                row.max_abs_error = SMAX(row.max_abs_error, abs_error);
                row.mean_rel_error += rel_error;
            }
            if (row.checked_particles > 0)
            {
                Real inv_count = 1.0 / static_cast<Real>(row.checked_particles);
                row.mean_abs_error *= inv_count;
                row.mean_rel_error *= inv_count;
            }
            verification_file << row.test_name << ","
                              << row.quantity_name << ","
                              << row.checked_particles << ","
                              << row.invalid_particles << ","
                              << row.mean_abs_error << ","
                              << row.max_abs_error << ","
                              << row.mean_rel_error << ","
                              << row.mean_projection_ratio << ","
                              << row.mean_alignment << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-verify] " << row.test_name
                      << " / " << row.quantity_name
                      << ": count=" << row.checked_particles
                      << ", invalid=" << row.invalid_particles
                      << ", mean_abs=" << row.mean_abs_error
                      << ", max_abs=" << row.max_abs_error
                      << ", mean_rel=" << row.mean_rel_error
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (has_first_invalid)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-verify] first invalid scalar sample/ref for "
                          << row.test_name << " / " << row.quantity_name
                          << ": sample=" << first_invalid_sample
                          << ", reference=" << first_invalid_reference
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
        };
        auto record_vector_verification =
            [&](const std::string &test_name,
                const std::string &quantity_name,
                const auto &sample_fn,
                const auto &reference_fn)
        {
            OperatorVerificationRow row;
            row.test_name = test_name;
            row.quantity_name = quantity_name;
            row.mean_projection_ratio = 0.0;
            row.mean_alignment = 0.0;
            bool has_first_invalid = false;
            Vecd first_invalid_sample = ZeroData<Vecd>::value;
            Vecd first_invalid_reference = ZeroData<Vecd>::value;
            for (size_t i = 0; i != total_plate_particles_count; ++i)
            {
                if (!select_plate_core(plate_positions[i]))
                {
                    continue;
                }
                Vecd reference_value = reference_fn(i);
                Vecd sample_value = sample_fn(i);
                if (!is_finite_vector(reference_value) || !is_finite_vector(sample_value))
                {
                    if (!has_first_invalid)
                    {
                        has_first_invalid = true;
                        first_invalid_sample = sample_value;
                        first_invalid_reference = reference_value;
                    }
                    row.invalid_particles++;
                    continue;
                }
                Real abs_error = (sample_value - reference_value).norm();
                Real reference_norm = reference_value.norm();
                Real sample_norm = sample_value.norm();
                Real rel_error = abs_error / (reference_norm + TinyReal);
                if (reference_norm > TinyReal)
                {
                    row.mean_projection_ratio +=
                        sample_value.dot(reference_value) /
                        (reference_norm * reference_norm);
                    row.mean_alignment +=
                        sample_value.dot(reference_value) /
                        ((sample_norm + TinyReal) * reference_norm);
                }
                row.checked_particles++;
                row.mean_abs_error += abs_error;
                row.max_abs_error = SMAX(row.max_abs_error, abs_error);
                row.mean_rel_error += rel_error;
            }
            if (row.checked_particles > 0)
            {
                Real inv_count = 1.0 / static_cast<Real>(row.checked_particles);
                row.mean_abs_error *= inv_count;
                row.mean_rel_error *= inv_count;
                row.mean_projection_ratio *= inv_count;
                row.mean_alignment *= inv_count;
            }
            verification_file << row.test_name << ","
                              << row.quantity_name << ","
                              << row.checked_particles << ","
                              << row.invalid_particles << ","
                              << row.mean_abs_error << ","
                              << row.max_abs_error << ","
                              << row.mean_rel_error << ","
                              << row.mean_projection_ratio << ","
                              << row.mean_alignment << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-verify] " << row.test_name
                      << " / " << row.quantity_name
                      << ": count=" << row.checked_particles
                      << ", invalid=" << row.invalid_particles
                      << ", mean_abs=" << row.mean_abs_error
                      << ", max_abs=" << row.max_abs_error
                      << ", mean_rel=" << row.mean_rel_error
                      << ", proj_ratio=" << row.mean_projection_ratio
                      << ", alignment=" << row.mean_alignment
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (has_first_invalid)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-verify] first invalid vector sample/ref for "
                          << row.test_name << " / " << row.quantity_name
                          << ": sample=(" << first_invalid_sample[0] << ","
                          << first_invalid_sample[1] << ","
                          << first_invalid_sample[2] << ")"
                          << ", reference=(" << first_invalid_reference[0] << ","
                          << first_invalid_reference[1] << ","
                          << first_invalid_reference[2] << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
        };
        auto record_vector_verification_with_selector =
            [&](const auto &selector_fn,
                const std::string &test_name,
                const std::string &quantity_name,
                const auto &sample_fn,
                const auto &reference_fn)
        {
            OperatorVerificationRow row;
            row.test_name = test_name;
            row.quantity_name = quantity_name;
            row.mean_projection_ratio = 0.0;
            row.mean_alignment = 0.0;
            bool has_first_invalid = false;
            Vecd first_invalid_sample = ZeroData<Vecd>::value;
            Vecd first_invalid_reference = ZeroData<Vecd>::value;
            for (size_t i = 0; i != total_plate_particles_count; ++i)
            {
                if (!selector_fn(plate_positions[i]))
                {
                    continue;
                }
                Vecd reference_value = reference_fn(i);
                Vecd sample_value = sample_fn(i);
                if (!is_finite_vector(reference_value) || !is_finite_vector(sample_value))
                {
                    if (!has_first_invalid)
                    {
                        has_first_invalid = true;
                        first_invalid_sample = sample_value;
                        first_invalid_reference = reference_value;
                    }
                    row.invalid_particles++;
                    continue;
                }
                Real abs_error = (sample_value - reference_value).norm();
                Real reference_norm = reference_value.norm();
                Real sample_norm = sample_value.norm();
                Real rel_error = abs_error / (reference_norm + TinyReal);
                if (reference_norm > TinyReal)
                {
                    row.mean_projection_ratio +=
                        sample_value.dot(reference_value) /
                        (reference_norm * reference_norm);
                    row.mean_alignment +=
                        sample_value.dot(reference_value) /
                        ((sample_norm + TinyReal) * reference_norm);
                }
                row.checked_particles++;
                row.mean_abs_error += abs_error;
                row.max_abs_error = SMAX(row.max_abs_error, abs_error);
                row.mean_rel_error += rel_error;
            }
            if (row.checked_particles > 0)
            {
                Real inv_count = 1.0 / static_cast<Real>(row.checked_particles);
                row.mean_abs_error *= inv_count;
                row.mean_rel_error *= inv_count;
                row.mean_projection_ratio *= inv_count;
                row.mean_alignment *= inv_count;
            }
            verification_file << row.test_name << ","
                              << row.quantity_name << ","
                              << row.checked_particles << ","
                              << row.invalid_particles << ","
                              << row.mean_abs_error << ","
                              << row.max_abs_error << ","
                              << row.mean_rel_error << ","
                              << row.mean_projection_ratio << ","
                              << row.mean_alignment << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-verify] " << row.test_name
                      << " / " << row.quantity_name
                      << ": count=" << row.checked_particles
                      << ", invalid=" << row.invalid_particles
                      << ", mean_abs=" << row.mean_abs_error
                      << ", max_abs=" << row.max_abs_error
                      << ", mean_rel=" << row.mean_rel_error
                      << ", proj_ratio=" << row.mean_projection_ratio
                      << ", alignment=" << row.mean_alignment
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (has_first_invalid)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-verify] first invalid vector sample/ref for "
                          << row.test_name << " / " << row.quantity_name
                          << ": sample=(" << first_invalid_sample[0] << ","
                          << first_invalid_sample[1] << ","
                          << first_invalid_sample[2] << ")"
                          << ", reference=(" << first_invalid_reference[0] << ","
                          << first_invalid_reference[1] << ","
                          << first_invalid_reference[2] << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
        };
        auto record_vector_verification_with_selector_on_body =
            [&](size_t total_particles,
                Vecd *positions,
                const auto &selector_fn,
                const std::string &test_name,
                const std::string &quantity_name,
                const auto &sample_fn,
                const auto &reference_fn)
        {
            OperatorVerificationRow row;
            row.test_name = test_name;
            row.quantity_name = quantity_name;
            row.mean_projection_ratio = 0.0;
            row.mean_alignment = 0.0;
            bool has_first_invalid = false;
            Vecd first_invalid_sample = ZeroData<Vecd>::value;
            Vecd first_invalid_reference = ZeroData<Vecd>::value;
            for (size_t i = 0; i != total_particles; ++i)
            {
                if (!selector_fn(positions[i]))
                {
                    continue;
                }
                Vecd reference_value = reference_fn(i);
                Vecd sample_value = sample_fn(i);
                if (!is_finite_vector(reference_value) || !is_finite_vector(sample_value))
                {
                    if (!has_first_invalid)
                    {
                        has_first_invalid = true;
                        first_invalid_sample = sample_value;
                        first_invalid_reference = reference_value;
                    }
                    row.invalid_particles++;
                    continue;
                }
                Real abs_error = (sample_value - reference_value).norm();
                Real reference_norm = reference_value.norm();
                Real sample_norm = sample_value.norm();
                Real rel_error = abs_error / (reference_norm + TinyReal);
                if (reference_norm > TinyReal)
                {
                    row.mean_projection_ratio +=
                        sample_value.dot(reference_value) /
                        (reference_norm * reference_norm);
                    row.mean_alignment +=
                        sample_value.dot(reference_value) /
                        ((sample_norm + TinyReal) * reference_norm);
                }
                row.checked_particles++;
                row.mean_abs_error += abs_error;
                row.max_abs_error = SMAX(row.max_abs_error, abs_error);
                row.mean_rel_error += rel_error;
            }
            if (row.checked_particles > 0)
            {
                Real inv_count = 1.0 / static_cast<Real>(row.checked_particles);
                row.mean_abs_error *= inv_count;
                row.mean_rel_error *= inv_count;
                row.mean_projection_ratio *= inv_count;
                row.mean_alignment *= inv_count;
            }
            verification_file << row.test_name << ","
                              << row.quantity_name << ","
                              << row.checked_particles << ","
                              << row.invalid_particles << ","
                              << row.mean_abs_error << ","
                              << row.max_abs_error << ","
                              << row.mean_rel_error << ","
                              << row.mean_projection_ratio << ","
                              << row.mean_alignment << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-verify] " << row.test_name
                      << " / " << row.quantity_name
                      << ": count=" << row.checked_particles
                      << ", invalid=" << row.invalid_particles
                      << ", mean_abs=" << row.mean_abs_error
                      << ", max_abs=" << row.max_abs_error
                      << ", mean_rel=" << row.mean_rel_error
                      << ", proj_ratio=" << row.mean_projection_ratio
                      << ", alignment=" << row.mean_alignment
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (has_first_invalid)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-verify] first invalid vector sample/ref for "
                          << row.test_name << " / " << row.quantity_name
                          << ": sample=(" << first_invalid_sample[0] << ","
                          << first_invalid_sample[1] << ","
                          << first_invalid_sample[2] << ")"
                          << ", reference=(" << first_invalid_reference[0] << ","
                          << first_invalid_reference[1] << ","
                          << first_invalid_reference[2] << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
        };

        // Test 1: curl(A) and curl(nu curl(A)) with quadratic Az(x,y).
        reset_plate_frequency_fields();
        auto reset_plate_contact_component_sources = [&]()
        {
            for (size_t i = 0; i != coil_particles.TotalRealParticles(); ++i)
            {
                coil_vector_potential_real[i] = ZeroData<Vecd>::value;
                coil_vector_potential_imag[i] = ZeroData<Vecd>::value;
            }
            for (size_t i = 0; i != air_particles.TotalRealParticles(); ++i)
            {
                air_vector_potential_real[i] = ZeroData<Vecd>::value;
                air_vector_potential_imag[i] = ZeroData<Vecd>::value;
            }
        };
        reset_plate_contact_component_sources();
        const Real curl_real_coeff = 1.0e-3;
        const Real curl_imag_coeff = -7.5e-4;
        ensure_plate_component_hessian_ck_prepared();
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            Vecd x_si = plate_local_position_si(plate_positions[i]);
            plate_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            plate_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        copy_operator_verify_a_real_components.exec();
        plate_verify_scalar_gradient_ax_ck.exec();
        plate_verify_scalar_hessian_ax_ck.exec();
        plate_verify_scalar_gradient_ay_ck.exec();
        plate_verify_scalar_hessian_ay_ck.exec();
        plate_verify_scalar_gradient_ck.exec();
        plate_verify_scalar_hessian_ck.exec();
        reconstruct_operator_verify_curl_nu_b_real_from_components.exec();
        auto sample_az_gradient_si = [&](size_t i, int axis) -> Real
        {
            return differential_operator_scaling *
                   plate_operator_verify_az_real_gradient[i][axis];
        };
        auto sample_az_hessian_si = [&](size_t i, int component) -> Real
        {
            return second_order_operator_scaling *
                   plate_operator_verify_az_real_hessian[i][component];
        };
        auto sample_az_laplacian_si = [&](size_t i) -> Real
        {
            const VecMatd &sample_hessian = plate_operator_verify_az_real_hessian[i];
            return second_order_operator_scaling *
                   (sample_hessian[0] + sample_hessian[1] + sample_hessian[2]);
        };
        record_scalar_function_verification(
            "scalar_hessian_operator", "gradAz_dx_real",
            [&](size_t i)
            { return sample_az_gradient_si(i, 0); },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return 2.0 * curl_real_coeff * x_si[0];
            });
        record_scalar_function_verification(
            "scalar_hessian_operator", "gradAz_dy_real",
            [&](size_t i)
            { return sample_az_gradient_si(i, 1); },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return 2.0 * curl_real_coeff * x_si[1];
            });
        record_scalar_function_verification(
            "scalar_hessian_operator", "hessAz_xx_real",
            [&](size_t i)
            { return sample_az_hessian_si(i, 0); },
            [&](size_t i)
            {
                (void)i;
                return 2.0 * curl_real_coeff;
            });
        record_scalar_function_verification(
            "scalar_hessian_operator", "hessAz_yy_real",
            [&](size_t i)
            { return sample_az_hessian_si(i, 1); },
            [&](size_t i)
            {
                (void)i;
                return 2.0 * curl_real_coeff;
            });
        record_scalar_function_verification(
            "scalar_hessian_operator", "laplacianAz_real",
            [&](size_t i)
            { return sample_az_laplacian_si(i); },
            [&](size_t i)
            {
                (void)i;
                return 4.0 * curl_real_coeff;
            });
        record_scalar_function_verification(
            "scalar_hessian_operator", "curlNuBz_from_hessian_real",
            [&](size_t i)
            { return -magnetic_reluctivity * sample_az_laplacian_si(i); },
            [&](size_t i)
            {
                (void)i;
                return -4.0 * curl_real_coeff * magnetic_reluctivity;
            });
        record_vector_verification(
            "component_hessian_identity", "curlNuB_real_from_component_hessians",
            [&](size_t i)
            { return plate_operator_verify_curl_nu_b_from_component_hessians_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
            });
        reset_plate_frequency_fields();
        reset_plate_contact_component_sources();
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            Vecd x_si = plate_local_position_si(plate_positions[i]);
            plate_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            plate_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        for (size_t i = 0; i != coil_particles.TotalRealParticles(); ++i)
        {
            Vecd x_si = (coil_positions[i] - plate_center) * geom_length_to_m;
            coil_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            coil_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        for (size_t i = 0; i != air_particles.TotalRealParticles(); ++i)
        {
            Vecd x_si = (air_positions[i] - plate_center) * geom_length_to_m;
            air_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            air_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        refresh_frequency_plate_magnetic_contact_by_component_hessian();
        record_vector_verification_with_selector(
            select_plate_contact_band,
            "component_hessian_contact_identity", "curlNuB_real_contact_band",
            [&](size_t i)
            { return plate_curl_nu_b_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
            });
        vector_potential_curl_real_plate_inner.exec();
        vector_potential_curl_imag_plate_inner.exec();
        refresh_frequency_plate_magnetic_inner_by_component_hessian();
        record_vector_verification(
            "curl_operator", "curlA_real",
            [&](size_t i)
            {
                const AngularVecd &sample = plate_vector_potential_curl_real[i];
                return Vecd(sample[0], sample[1], sample[2]);
            },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return Vecd(2.0 * curl_real_coeff * x_si[1],
                            -2.0 * curl_real_coeff * x_si[0],
                            0.0);
            });
        record_vector_verification(
            "curl_operator", "curlA_imag",
            [&](size_t i)
            {
                const AngularVecd &sample = plate_vector_potential_curl_imag[i];
                return Vecd(sample[0], sample[1], sample[2]);
            },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return Vecd(2.0 * curl_imag_coeff * x_si[1],
                            -2.0 * curl_imag_coeff * x_si[0],
                            0.0);
            });
        record_vector_verification(
            "curl_nu_b_operator", "curlNuB_real",
            [&](size_t i)
            { return plate_curl_nu_b_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
            });
        record_vector_verification(
            "curl_nu_b_operator", "curlNuB_imag",
            [&](size_t i)
            { return plate_curl_nu_b_imag[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_imag_coeff * magnetic_reluctivity);
            });
        record_scalar_function_verification(
            "curl_nu_b_crosscheck", "curlNuBz_real_vs_hessian",
            [&](size_t i)
            { return plate_curl_nu_b_real[i][2]; },
            [&](size_t i)
            { return -magnetic_reluctivity * sample_az_laplacian_si(i); });
        record_vector_verification(
            "curl_nu_b_crosscheck", "curlNuB_real_vs_component_hessians",
            [&](size_t i)
            { return plate_curl_nu_b_real[i]; },
            [&](size_t i)
            { return plate_operator_verify_curl_nu_b_from_component_hessians_real[i]; });

        // Test 1b: coil curl(A) and curl(nu curl(A)) under the same manufactured field.
        reset_plate_frequency_fields();
        reset_plate_contact_component_sources();
        size_t total_coil_particles_count = coil_particles.TotalRealParticles();
        for (size_t i = 0; i != total_coil_particles_count; ++i)
        {
            coil_vector_potential_change_rate_real[i] = ZeroData<Vecd>::value;
            coil_vector_potential_change_rate_imag[i] = ZeroData<Vecd>::value;
            coil_curl_nu_b_real[i] = ZeroData<Vecd>::value;
            coil_curl_nu_b_imag[i] = ZeroData<Vecd>::value;
            coil_vector_potential_curl_real[i] = ZeroData<AngularVecd>::value;
            coil_vector_potential_curl_imag[i] = ZeroData<AngularVecd>::value;
            Vecd x_si = coil_local_position_si(coil_positions[i]);
            coil_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            coil_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            Vecd x_si = coil_local_position_si(plate_positions[i]);
            plate_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            plate_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        for (size_t i = 0; i != total_air_particles_count; ++i)
        {
            Vecd x_si = coil_local_position_si(air_positions[i]);
            air_vector_potential_real[i] =
                Vecd(0.0, 0.0, curl_real_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
            air_vector_potential_imag[i] =
                Vecd(0.0, 0.0, curl_imag_coeff * (x_si[0] * x_si[0] + x_si[1] * x_si[1]));
        }
        refresh_frequency_plate_magnetic_terms();
        refresh_frequency_air_magnetic_terms();
        refresh_frequency_coil_magnetic_terms();
        record_vector_verification_with_selector_on_body(
            total_coil_particles_count, coil_positions, select_coil_core,
            "coil_curl_operator", "curlA_real",
            [&](size_t i)
            {
                const AngularVecd &sample = coil_vector_potential_curl_real[i];
                return Vecd(sample[0], sample[1], sample[2]);
            },
            [&](size_t i)
            {
                Vecd x_si = coil_local_position_si(coil_positions[i]);
                return Vecd(2.0 * curl_real_coeff * x_si[1],
                            -2.0 * curl_real_coeff * x_si[0], 0.0);
            });
        record_vector_verification_with_selector_on_body(
            total_coil_particles_count, coil_positions, select_coil_core,
            "coil_curl_nu_b_operator", "curlNuB_real",
            [&](size_t i)
            { return coil_curl_nu_b_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
            });
        record_vector_verification_with_selector_on_body(
            total_coil_particles_count, coil_positions, select_coil_contact_band,
            "coil_curl_nu_b_operator", "curlNuB_real_contact_band",
            [&](size_t i)
            { return coil_curl_nu_b_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
            });
        if (use_frequency_coil_component_hessian_inner)
        {
            refresh_frequency_coil_magnetic_inner_by_component_hessian();
            record_vector_verification_with_selector_on_body(
                total_coil_particles_count, coil_positions, select_coil_core,
                "coil_component_hessian_inner_identity", "curlNuB_real",
                [&](size_t i)
                { return coil_curl_nu_b_real[i]; },
                [&](size_t i)
                {
                    (void)i;
                    return Vecd(0.0, 0.0, -4.0 * curl_real_coeff * magnetic_reluctivity);
                });
        }

        // Test 2: scalar source and gradient with linear phi and quadratic Ax.
        reset_plate_frequency_fields();
        const Vecd grad_phi_real_ref(0.35, -0.22, 0.18);
        const Vecd grad_phi_imag_ref(-0.14, 0.27, -0.09);
        const Real source_real_coeff = 3.0e-4;
        const Real source_imag_coeff = -2.5e-4;
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            Vecd x_si = plate_local_position_si(plate_positions[i]);
            plate_electric_potential_real[i] = grad_phi_real_ref.dot(x_si);
            plate_electric_potential_imag[i] = grad_phi_imag_ref.dot(x_si);
            plate_vector_potential_imag[i] =
                Vecd(source_real_coeff * x_si[0] * x_si[0], 0.0, 0.0);
            plate_vector_potential_real[i] =
                Vecd(source_imag_coeff * x_si[0] * x_si[0], 0.0, 0.0);
        }
        electric_potential_source_real_plate.exec();
        electric_potential_source_imag_plate.exec();
        electric_potential_gradient_real_plate_inner.exec();
        electric_potential_gradient_imag_plate_inner.exec();
        record_vector_verification(
            "phi_gradient_operator", "gradPhi_real",
            [&](size_t i)
            { return plate_electric_potential_gradient_real[i]; },
            [&](size_t i)
            {
                (void)i;
                return grad_phi_real_ref;
            });
        record_vector_verification(
            "phi_gradient_operator", "gradPhi_imag",
            [&](size_t i)
            { return plate_electric_potential_gradient_imag[i]; },
            [&](size_t i)
            {
                (void)i;
                return grad_phi_imag_ref;
            });
        record_scalar_verification(
            "phi_source_operator", "source_real",
            plate_source_real,
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return 2.0 * plate_sigma[i] * harmonic_angular_frequency_runtime *
                       source_real_coeff * x_si[0];
            });
        record_scalar_verification(
            "phi_source_operator", "source_imag",
            plate_source_imag,
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                return -2.0 * plate_sigma[i] * harmonic_angular_frequency_runtime *
                       source_imag_coeff * x_si[0];
            });

        // Test 3: constitutive closure E/J/Q from prescribed A and phi.
        reset_plate_frequency_fields();
        const Real ejq_a_real_coeff = 1.8e-4;
        const Real ejq_a_imag_coeff = -2.1e-4;
        const Vecd ejq_grad_phi_real(0.16, -0.11, 0.07);
        const Vecd ejq_grad_phi_imag(-0.09, 0.13, -0.05);
        for (size_t i = 0; i != total_plate_particles_count; ++i)
        {
            Vecd x_si = plate_local_position_si(plate_positions[i]);
            plate_electric_potential_real[i] = ejq_grad_phi_real.dot(x_si);
            plate_electric_potential_imag[i] = ejq_grad_phi_imag.dot(x_si);
            plate_vector_potential_real[i] = Vecd(0.0, 0.0, ejq_a_real_coeff * x_si[0]);
            plate_vector_potential_imag[i] = Vecd(0.0, 0.0, ejq_a_imag_coeff * x_si[1]);
        }
        electric_potential_gradient_real_plate_inner.exec();
        electric_potential_gradient_imag_plate_inner.exec();
        electric_field_current_heat_frequency.exec();
        record_vector_verification(
            "ejq_closure", "E_real",
            [&](size_t i)
            { return plate_electric_field_real[i]; },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                Vecd a_imag(0.0, 0.0, ejq_a_imag_coeff * x_si[1]);
                Vecd e_real =
                    harmonic_angular_frequency_runtime * a_imag - ejq_grad_phi_real;
                return e_real;
            });
        record_vector_verification(
            "ejq_closure", "E_imag",
            [&](size_t i)
            { return plate_electric_field_imag[i]; },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                Vecd a_real(0.0, 0.0, ejq_a_real_coeff * x_si[0]);
                Vecd e_imag =
                    -harmonic_angular_frequency_runtime * a_real - ejq_grad_phi_imag;
                return e_imag;
            });
        record_vector_verification(
            "ejq_closure", "J_real",
            [&](size_t i)
            { return plate_current_density_real[i]; },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                Vecd a_imag(0.0, 0.0, ejq_a_imag_coeff * x_si[1]);
                Vecd e_real =
                    harmonic_angular_frequency_runtime * a_imag - ejq_grad_phi_real;
                Vecd j_real = plate_sigma[i] * e_real;
                return j_real;
            });
        record_vector_verification(
            "ejq_closure", "J_imag",
            [&](size_t i)
            { return plate_current_density_imag[i]; },
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                Vecd a_real(0.0, 0.0, ejq_a_real_coeff * x_si[0]);
                Vecd e_imag =
                    -harmonic_angular_frequency_runtime * a_real - ejq_grad_phi_imag;
                Vecd j_imag = plate_sigma[i] * e_imag;
                return j_imag;
            });
        record_scalar_verification(
            "ejq_closure", "joule_heat",
            plate_joule_heat,
            [&](size_t i)
            {
                Vecd x_si = plate_local_position_si(plate_positions[i]);
                Vecd a_real(0.0, 0.0, ejq_a_real_coeff * x_si[0]);
                Vecd a_imag(0.0, 0.0, ejq_a_imag_coeff * x_si[1]);
                Vecd e_real =
                    harmonic_angular_frequency_runtime * a_imag - ejq_grad_phi_real;
                Vecd e_imag =
                    -harmonic_angular_frequency_runtime * a_real - ejq_grad_phi_imag;
                return static_cast<Real>(0.5) * plate_sigma[i] *
                       (e_real.squaredNorm() + e_imag.squaredNorm());
            });

        std::cout << "[team7-verify] operator verification written to "
                  << io_environment.OutputFolder() + "/team7_operator_verification.csv"
                  << std::endl;
        return 0;
    };
    int operator_verification_status = run_frequency_operator_verification();
    if (operator_verification_status >= 0)
    {
        return operator_verification_status;
    }
    if (debug_coil_nan_detection && !use_frequency_aphi)
    {
        report_coil_nan("initialization");
        report_plate_air_nan("initialization");
    }

    std::ofstream plate_diagnostics_file(
        io_environment.OutputFolder() + "/team7_plate_diagnostics.csv",
        std::ios::out | std::ios::trunc);
    plate_diagnostics_file << std::setprecision(12);
    plate_diagnostics_file
        << "time,total_joule_power_si,avg_joule_heat,max_joule_heat,"
        << "avg_current_density_magnitude,max_current_density_magnitude,"
        << "avg_temperature,max_temperature,"
        << "max_a_rate,max_a_dot,max_grad_phi,max_electric_field,"
        << "em_iterations,em_max_rate,em_residual_ratio,dt_em_current,"
        << "total_volume_si,plate_avg_residual_norm,coil_avg_residual_norm,air_avg_residual_norm,"
        << "plate_interface_inner_avg_residual,plate_interface_air_avg_residual,"
        << "coil_interface_inner_avg_residual,coil_interface_air_avg_residual,"
        << "coil_equivalent_ampere_turns_real,coil_equivalent_ampere_turns_imag\n";

    std::ofstream team7_observables_file(
        io_environment.OutputFolder() + "/team7_observables.csv",
        std::ios::out | std::ios::trunc);
    team7_observables_file << std::setprecision(12);
    team7_observables_file
        << "time,total_joule_power_si,avg_temperature,max_temperature,"
        << "max_a_rate,max_a_dot,max_grad_phi,max_electric_field,"
        << "em_iterations,em_max_rate,em_residual_ratio,dt_em_current,"
        << "B_air_above_plate_x,B_air_above_plate_y,B_air_above_plate_z,B_air_above_plate_mag,"
        << "B_air_below_plate_x,B_air_below_plate_y,B_air_below_plate_z,B_air_below_plate_mag,"
        << "B_plate_center_x,B_plate_center_y,B_plate_center_z,B_plate_center_mag,"
        << "T_plate_center,T_plate_edge_xplus,"
        << "plate_avg_residual_norm,coil_avg_residual_norm,air_avg_residual_norm,"
        << "plate_interface_inner_avg_residual,plate_interface_air_avg_residual,"
        << "coil_interface_inner_avg_residual,coil_interface_air_avg_residual,"
        << "coil_equivalent_ampere_turns_real,coil_equivalent_ampere_turns_imag\n";

    std::ofstream em_iteration_diagnostics_file;
    if (em_only_diagnostics_mode)
    {
        em_iteration_diagnostics_file.open(
            io_environment.OutputFolder() + "/team7_em_iteration_diagnostics.csv",
            std::ios::out | std::ios::trunc);
        em_iteration_diagnostics_file << std::setprecision(12);
        em_iteration_diagnostics_file
            << "em_iteration_global,outer_iteration,outer_time,em_iteration_local,em_substep_time,"
            << "dt_em_step,em_max_rate,em_residual_ratio,total_joule_power_si,avg_joule_heat,max_joule_heat,"
            << "avg_current_density_magnitude,max_current_density_magnitude,"
            << "max_a_rate,max_a_dot,max_grad_phi,max_electric_field,"
            << "coil_configured_ampere_turns_real,coil_configured_ampere_turns_imag,"
            << "coil_equivalent_ampere_turns_real,coil_equivalent_ampere_turns_imag,"
            << "coil_avg_source_density_magnitude,coil_max_source_density_magnitude,"
            << "plate_avg_source,plate_avg_curl_nu_b,plate_avg_sigma_grad_phi,plate_avg_omega_sigma_a,plate_avg_residual,"
            << "coil_avg_source,coil_avg_curl_nu_b,coil_avg_sigma_grad_phi,coil_avg_omega_sigma_a,coil_avg_residual,"
            << "air_avg_curl_nu_b,air_avg_residual,"
            << "air_curl_inner_avg,air_curl_contact_avg,air_curl_inner_max,air_curl_contact_max,"
            << "air_use_inner_term,air_use_contact_term,"
            << "plate_max_residual,coil_max_residual,air_max_residual,"
            << "plate_interface_inner_avg_curl,plate_interface_inner_avg_residual,plate_interface_inner_volume,"
            << "plate_interface_air_avg_curl,plate_interface_air_avg_residual,plate_interface_air_volume,"
            << "coil_interface_inner_avg_curl,coil_interface_inner_avg_residual,coil_interface_inner_volume,"
            << "coil_interface_air_avg_curl,coil_interface_air_avg_residual,coil_interface_air_volume,"
            << "probe_b_air_above_mag,probe_b_air_below_mag,"
            << "probe_b_air_above_ref_mag,probe_b_air_below_ref_mag,"
            << "probe_b_air_above_sim_over_ref,probe_b_air_below_sim_over_ref,"
            << "runtime_source_scale,vacuum_reluctivity_scale\n";
    }

    Real next_output_time = 0.0;
    size_t iteration = 0;
    bool wrote_last_step_recording = false;
    size_t em_iteration_global = 0;
    Real dt_em_current = SMIN(dt_em_max, SMAX(dt_em_min, dt_em));
    size_t em_iterations_last = 0;
    Real em_max_rate_last = 0.0;
    Real em_residual_ratio_last = 0.0;
    bool frequency_em_converged = false;
    bool stop_due_to_frequency_hard_guard = false;
    Real em_max_relaxation_scaling =
        SMAX(coil_a_relaxation_scaling, SMAX(plate_a_relaxation_scaling, air_a_relaxation_scaling));
    Real em_rate_to_delta_scaling =
        em_rate_use_effective_update ? 1.0 : em_max_relaxation_scaling;
    Real coil_effective_rate_limit = coil_a_rate_limit;
    Real plate_effective_rate_limit = plate_a_rate_limit;
    Real air_effective_rate_limit = air_a_rate_limit;
    if (em_rate_use_effective_update)
    {
        coil_effective_rate_limit *= coil_a_relaxation_scaling;
        plate_effective_rate_limit *= plate_a_relaxation_scaling;
        air_effective_rate_limit *= air_a_relaxation_scaling;
    }
    Real effective_rate_limit_max =
        use_frequency_operator_coil_air_only
            ? coil_effective_rate_limit
            : SMAX(coil_effective_rate_limit, plate_effective_rate_limit);
    if (em_rate_include_air)
    {
        effective_rate_limit_max = SMAX(effective_rate_limit_max, air_effective_rate_limit);
    }
    Real em_frequency_rate_saturation_limit = effective_rate_limit_max;

    PlateDiagnostics initial_diagnostics = evaluate_plate_diagnostics();
    EmTermBreakdown initial_em_term_breakdown = evaluate_em_term_breakdown();
    EmInterfaceShellBreakdown initial_em_interface_breakdown = evaluate_em_interface_shell_breakdown();
    CoilSourceDiagnostics initial_coil_source_diagnostics = evaluate_coil_source_diagnostics();
    Vec3d b_air_above_initial = sample_magnetic_flux_density(air_vector_potential_curl_probe, probe_b_air_above_plate);
    Vec3d b_air_below_initial = sample_magnetic_flux_density(air_vector_potential_curl_probe, probe_b_air_below_plate);
    Vec3d b_plate_center_initial = sample_magnetic_flux_density(plate_vector_potential_curl_probe, probe_b_plate_center);
    Real b_air_above_initial_mag =
        sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_above_plate);
    Real b_air_below_initial_mag =
        sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_below_plate);
    Real b_plate_center_initial_mag =
        sample_magnetic_flux_density_magnitude(plate_vector_potential_curl_probe, plate_vector_potential_curl_imag, probe_b_plate_center);
    Real t_plate_center_initial = plate_temperature[probe_t_plate_center.particle_index];
    Real t_plate_edge_xplus_initial = plate_temperature[probe_t_plate_edge_xplus.particle_index];
    std::pair<Vec3d, Vec3d> b_air_above_initial_ref_field =
        compute_biot_savart_magnetic_flux_density(probe_b_air_above_plate.sampled_position);
    std::pair<Vec3d, Vec3d> b_air_below_initial_ref_field =
        compute_biot_savart_magnetic_flux_density(probe_b_air_below_plate.sampled_position);
    Real b_air_above_initial_ref_mag =
        complex_vector_magnitude(b_air_above_initial_ref_field.first,
                                 b_air_above_initial_ref_field.second);
    Real b_air_below_initial_ref_mag =
        complex_vector_magnitude(b_air_below_initial_ref_field.first,
                                 b_air_below_initial_ref_field.second);
    Real b_air_above_initial_sim_over_ref =
        b_air_above_initial_mag / (b_air_above_initial_ref_mag + TinyReal);
    Real b_air_below_initial_sim_over_ref =
        b_air_below_initial_mag / (b_air_below_initial_ref_mag + TinyReal);

    plate_diagnostics_file << 0.0 << ","
                           << initial_diagnostics.total_joule_power << ","
                           << initial_diagnostics.avg_joule_heat << ","
                           << initial_diagnostics.max_joule_heat << ","
                           << initial_diagnostics.avg_current_density_magnitude << ","
                           << initial_diagnostics.max_current_density_magnitude << ","
                           << initial_diagnostics.avg_temperature << ","
                           << initial_diagnostics.max_temperature << ","
                           << initial_diagnostics.max_a_rate << ","
                           << initial_diagnostics.max_a_dot << ","
                           << initial_diagnostics.max_grad_phi << ","
                           << initial_diagnostics.max_electric_field << ","
                           << 0 << ","
                           << 0.0 << ","
                           << 0.0 << ","
                           << dt_em_current << ","
                           << initial_diagnostics.total_volume << ","
                           << initial_em_term_breakdown.plate.avg_residual << ","
                           << initial_em_term_breakdown.coil.avg_residual << ","
                           << initial_em_term_breakdown.air.avg_residual << ","
                           << initial_em_interface_breakdown.plate_inner.avg_residual << ","
                           << initial_em_interface_breakdown.plate_air.avg_residual << ","
                           << initial_em_interface_breakdown.coil_inner.avg_residual << ","
                           << initial_em_interface_breakdown.coil_air.avg_residual << ","
                           << initial_coil_source_diagnostics.equivalent_ampere_turns_real << ","
                           << initial_coil_source_diagnostics.equivalent_ampere_turns_imag << "\n";
    team7_observables_file << 0.0 << ","
                           << initial_diagnostics.total_joule_power << ","
                           << initial_diagnostics.avg_temperature << ","
                           << initial_diagnostics.max_temperature << ","
                           << initial_diagnostics.max_a_rate << ","
                           << initial_diagnostics.max_a_dot << ","
                           << initial_diagnostics.max_grad_phi << ","
                           << initial_diagnostics.max_electric_field << ","
                           << 0 << ","
                           << 0.0 << ","
                           << 0.0 << ","
                           << dt_em_current << ","
                           << b_air_above_initial[0] << ","
                           << b_air_above_initial[1] << ","
                           << b_air_above_initial[2] << ","
                           << b_air_above_initial_mag << ","
                           << b_air_below_initial[0] << ","
                           << b_air_below_initial[1] << ","
                           << b_air_below_initial[2] << ","
                           << b_air_below_initial_mag << ","
                           << b_plate_center_initial[0] << ","
                           << b_plate_center_initial[1] << ","
                           << b_plate_center_initial[2] << ","
                           << b_plate_center_initial_mag << ","
                           << t_plate_center_initial << ","
                           << t_plate_edge_xplus_initial << ","
                           << initial_em_term_breakdown.plate.avg_residual << ","
                           << initial_em_term_breakdown.coil.avg_residual << ","
                           << initial_em_term_breakdown.air.avg_residual << ","
                           << initial_em_interface_breakdown.plate_inner.avg_residual << ","
                           << initial_em_interface_breakdown.plate_air.avg_residual << ","
                           << initial_em_interface_breakdown.coil_inner.avg_residual << ","
                           << initial_em_interface_breakdown.coil_air.avg_residual << ","
                           << initial_coil_source_diagnostics.equivalent_ampere_turns_real << ","
                           << initial_coil_source_diagnostics.equivalent_ampere_turns_imag << "\n";
    if (em_only_diagnostics_mode)
    {
        em_iteration_diagnostics_file << 0 << ","
                                      << 0 << ","
                                      << 0.0 << ","
                                      << 0 << ","
                                      << 0.0 << ","
                                      << dt_em_current << ","
                                      << 0.0 << ","
                                      << 0.0 << ","
                                      << initial_diagnostics.total_joule_power << ","
                                      << initial_diagnostics.avg_joule_heat << ","
                                      << initial_diagnostics.max_joule_heat << ","
                                      << initial_diagnostics.avg_current_density_magnitude << ","
                                      << initial_diagnostics.max_current_density_magnitude << ","
                                      << initial_diagnostics.max_a_rate << ","
                                      << initial_diagnostics.max_a_dot << ","
                                      << initial_diagnostics.max_grad_phi << ","
                                      << initial_diagnostics.max_electric_field << ","
                                      << initial_coil_source_diagnostics.configured_ampere_turns_real << ","
                                      << initial_coil_source_diagnostics.configured_ampere_turns_imag << ","
                                      << initial_coil_source_diagnostics.equivalent_ampere_turns_real << ","
                                      << initial_coil_source_diagnostics.equivalent_ampere_turns_imag << ","
                                      << initial_coil_source_diagnostics.avg_source_density_magnitude << ","
                                      << initial_coil_source_diagnostics.max_source_density_magnitude << ","
                                      << initial_em_term_breakdown.plate.avg_source << ","
                                      << initial_em_term_breakdown.plate.avg_curl_nu_b << ","
                                      << initial_em_term_breakdown.plate.avg_sigma_grad_phi << ","
                                      << initial_em_term_breakdown.plate.avg_omega_sigma_a << ","
                                      << initial_em_term_breakdown.plate.avg_residual << ","
                                      << initial_em_term_breakdown.coil.avg_source << ","
                                      << initial_em_term_breakdown.coil.avg_curl_nu_b << ","
                                      << initial_em_term_breakdown.coil.avg_sigma_grad_phi << ","
                                      << initial_em_term_breakdown.coil.avg_omega_sigma_a << ","
                                      << initial_em_term_breakdown.coil.avg_residual << ","
                                      << initial_em_term_breakdown.air.avg_curl_nu_b << ","
                                      << initial_em_term_breakdown.air.avg_residual << ","
                                      << air_curl_inner_avg_norm_last << ","
                                      << air_curl_contact_avg_norm_last << ","
                                      << air_curl_inner_max_norm_last << ","
                                      << air_curl_contact_max_norm_last << ","
                                      << static_cast<int>(use_frequency_air_inner_term) << ","
                                      << static_cast<int>(use_frequency_air_contact_term) << ","
                                      << initial_em_term_breakdown.plate.max_residual << ","
                                      << initial_em_term_breakdown.coil.max_residual << ","
                                      << initial_em_term_breakdown.air.max_residual << ","
                                      << initial_em_interface_breakdown.plate_inner.avg_curl_nu_b << ","
                                      << initial_em_interface_breakdown.plate_inner.avg_residual << ","
                                      << initial_em_interface_breakdown.plate_inner.total_volume << ","
                                      << initial_em_interface_breakdown.plate_air.avg_curl_nu_b << ","
                                      << initial_em_interface_breakdown.plate_air.avg_residual << ","
                                      << initial_em_interface_breakdown.plate_air.total_volume << ","
                                      << initial_em_interface_breakdown.coil_inner.avg_curl_nu_b << ","
                                      << initial_em_interface_breakdown.coil_inner.avg_residual << ","
                                      << initial_em_interface_breakdown.coil_inner.total_volume << ","
                                      << initial_em_interface_breakdown.coil_air.avg_curl_nu_b << ","
                                      << initial_em_interface_breakdown.coil_air.avg_residual << ","
                                      << initial_em_interface_breakdown.coil_air.total_volume << ","
                                      << b_air_above_initial_mag << ","
                                      << b_air_below_initial_mag << ","
                                      << b_air_above_initial_ref_mag << ","
                                      << b_air_below_initial_ref_mag << ","
                                      << b_air_above_initial_sim_over_ref << ","
                                      << b_air_below_initial_sim_over_ref << ","
                                      << runtime_source_scale << ","
                                      << vacuum_reluctivity_scale_runtime << "\n";
    }
    auto write_em_iteration_diagnostics_row =
        [&](size_t em_iteration_local,
            Real em_substep_time,
            Real dt_em_step,
            Real em_max_rate_value,
            Real em_residual_ratio_value)
    {
        if (!em_only_diagnostics_mode)
        {
            return;
        }
        PlateDiagnostics diagnostics = evaluate_plate_diagnostics();
        CoilSourceDiagnostics coil_source_diagnostics = evaluate_coil_source_diagnostics();
        EmTermBreakdown em_term_breakdown = evaluate_em_term_breakdown();
        EmInterfaceShellBreakdown em_interface_breakdown = evaluate_em_interface_shell_breakdown();
        Real b_air_above_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_above_plate);
        Real b_air_below_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe,
                                                   air_vector_potential_curl_imag,
                                                   probe_b_air_below_plate);
        std::pair<Vec3d, Vec3d> b_air_above_ref_field =
            compute_biot_savart_magnetic_flux_density(probe_b_air_above_plate.sampled_position);
        std::pair<Vec3d, Vec3d> b_air_below_ref_field =
            compute_biot_savart_magnetic_flux_density(probe_b_air_below_plate.sampled_position);
        Real b_air_above_ref_mag =
            complex_vector_magnitude(b_air_above_ref_field.first,
                                     b_air_above_ref_field.second);
        Real b_air_below_ref_mag =
            complex_vector_magnitude(b_air_below_ref_field.first,
                                     b_air_below_ref_field.second);
        Real b_air_above_sim_over_ref =
            b_air_above_mag / (b_air_above_ref_mag + TinyReal);
        Real b_air_below_sim_over_ref =
            b_air_below_mag / (b_air_below_ref_mag + TinyReal);
        em_iteration_diagnostics_file << em_iteration_global << ","
                                      << iteration << ","
                                      << physical_time << ","
                                      << em_iteration_local << ","
                                      << em_substep_time << ","
                                      << dt_em_step << ","
                                      << em_max_rate_value << ","
                                      << em_residual_ratio_value << ","
                                      << diagnostics.total_joule_power << ","
                                      << diagnostics.avg_joule_heat << ","
                                      << diagnostics.max_joule_heat << ","
                                      << diagnostics.avg_current_density_magnitude << ","
                                      << diagnostics.max_current_density_magnitude << ","
                                      << diagnostics.max_a_rate << ","
                                      << diagnostics.max_a_dot << ","
                                      << diagnostics.max_grad_phi << ","
                                      << diagnostics.max_electric_field << ","
                                      << coil_source_diagnostics.configured_ampere_turns_real << ","
                                      << coil_source_diagnostics.configured_ampere_turns_imag << ","
                                      << coil_source_diagnostics.equivalent_ampere_turns_real << ","
                                      << coil_source_diagnostics.equivalent_ampere_turns_imag << ","
                                      << coil_source_diagnostics.avg_source_density_magnitude << ","
                                      << coil_source_diagnostics.max_source_density_magnitude << ","
                                      << em_term_breakdown.plate.avg_source << ","
                                      << em_term_breakdown.plate.avg_curl_nu_b << ","
                                      << em_term_breakdown.plate.avg_sigma_grad_phi << ","
                                      << em_term_breakdown.plate.avg_omega_sigma_a << ","
                                      << em_term_breakdown.plate.avg_residual << ","
                                      << em_term_breakdown.coil.avg_source << ","
                                      << em_term_breakdown.coil.avg_curl_nu_b << ","
                                      << em_term_breakdown.coil.avg_sigma_grad_phi << ","
                                      << em_term_breakdown.coil.avg_omega_sigma_a << ","
                                      << em_term_breakdown.coil.avg_residual << ","
                                      << em_term_breakdown.air.avg_curl_nu_b << ","
                                      << em_term_breakdown.air.avg_residual << ","
                                      << air_curl_inner_avg_norm_last << ","
                                      << air_curl_contact_avg_norm_last << ","
                                      << air_curl_inner_max_norm_last << ","
                                      << air_curl_contact_max_norm_last << ","
                                      << static_cast<int>(use_frequency_air_inner_term) << ","
                                      << static_cast<int>(use_frequency_air_contact_term) << ","
                                      << em_term_breakdown.plate.max_residual << ","
                                      << em_term_breakdown.coil.max_residual << ","
                                      << em_term_breakdown.air.max_residual << ","
                                      << em_interface_breakdown.plate_inner.avg_curl_nu_b << ","
                                      << em_interface_breakdown.plate_inner.avg_residual << ","
                                      << em_interface_breakdown.plate_inner.total_volume << ","
                                      << em_interface_breakdown.plate_air.avg_curl_nu_b << ","
                                      << em_interface_breakdown.plate_air.avg_residual << ","
                                      << em_interface_breakdown.plate_air.total_volume << ","
                                      << em_interface_breakdown.coil_inner.avg_curl_nu_b << ","
                                      << em_interface_breakdown.coil_inner.avg_residual << ","
                                      << em_interface_breakdown.coil_inner.total_volume << ","
                                      << em_interface_breakdown.coil_air.avg_curl_nu_b << ","
                                      << em_interface_breakdown.coil_air.avg_residual << ","
                                      << em_interface_breakdown.coil_air.total_volume << ","
                                      << b_air_above_mag << ","
                                      << b_air_below_mag << ","
                                      << b_air_above_ref_mag << ","
                                      << b_air_below_ref_mag << ","
                                      << b_air_above_sim_over_ref << ","
                                      << b_air_below_sim_over_ref << ","
                                      << runtime_source_scale << ","
                                      << vacuum_reluctivity_scale_runtime << "\n";
    };

    auto execute_frequency_em_substep = [&](Real dt_em_step)
    {
        for (size_t block_sweep = 0; block_sweep < frequency_block_gs_sweeps; ++block_sweep)
        {
        Real dt_block_step =
            dt_em_step / static_cast<Real>(SMAX(static_cast<size_t>(1), frequency_block_gs_sweeps));
        if (use_circular_coil_source)
        {
            set_frequency_source_on_coil_circular.exec();
        }
        else
        {
            set_frequency_source_on_coil.exec();
        }
        if (apply_runtime_source_scale)
        {
            scale_frequency_source_on_coil.exec();
        }

        if (use_frequency_coil_scalar_potential)
        {
            electric_potential_source_real_coil.exec();
        }
        if (!use_frequency_operator_coil_air_only)
        {
            electric_potential_source_real_plate.exec();
        }
        if (use_frequency_scalar_contact_coupling)
        {
            if (use_frequency_coil_scalar_potential)
            {
                electric_potential_source_real_coil_contact.exec();
            }
            if (!use_frequency_operator_coil_air_only)
            {
                electric_potential_source_real_plate_contact.exec();
            }
        }
        if (use_frequency_air_scalar_potential)
        {
            electric_potential_source_real_air.exec();
            if (use_frequency_scalar_contact_coupling)
            {
                electric_potential_source_real_air_contact.exec();
            }
        }

        if (use_frequency_coil_scalar_potential)
        {
            electric_potential_source_imag_coil.exec();
        }
        if (!use_frequency_operator_coil_air_only)
        {
            electric_potential_source_imag_plate.exec();
        }
        if (use_frequency_scalar_contact_coupling)
        {
            if (use_frequency_coil_scalar_potential)
            {
                electric_potential_source_imag_coil_contact.exec();
            }
            if (!use_frequency_operator_coil_air_only)
            {
                electric_potential_source_imag_plate_contact.exec();
            }
        }
        if (use_frequency_air_scalar_potential)
        {
            electric_potential_source_imag_air.exec();
            if (use_frequency_scalar_contact_coupling)
            {
                electric_potential_source_imag_air_contact.exec();
            }
        }

        for (size_t k = 0; k < phi_relax_iterations; ++k)
        {
            if (use_frequency_air_scalar_potential)
            {
                electric_potential_real_relaxation_air_inner.exec();
                electric_potential_real_relaxation_air_update.exec(dt_block_step);
                constrain_phi_real_reference_air.exec();
                if (constrain_phi_boundary_air_enabled)
                {
                    constrain_phi_real_boundary_air.exec();
                }
                if (constrain_phi_all_air_enabled)
                {
                    constrain_phi_real_all_air.exec();
                }

                electric_potential_imag_relaxation_air_inner.exec();
                electric_potential_imag_relaxation_air_update.exec(dt_block_step);
                constrain_phi_imag_reference_air.exec();
                if (constrain_phi_boundary_air_enabled)
                {
                    constrain_phi_imag_boundary_air.exec();
                }
                if (constrain_phi_all_air_enabled)
                {
                    constrain_phi_imag_all_air.exec();
                }
            }

            if (use_frequency_coil_scalar_potential)
            {
                electric_potential_real_relaxation_coil_inner.exec();
                electric_potential_real_relaxation_coil_update.exec(dt_block_step);
                constrain_phi_real_reference_coil.exec();
                if (constrain_phi_boundary_coil_enabled)
                {
                    constrain_phi_real_boundary_coil.exec();
                }

                electric_potential_imag_relaxation_coil_inner.exec();
                electric_potential_imag_relaxation_coil_update.exec(dt_block_step);
                constrain_phi_imag_reference_coil.exec();
                if (constrain_phi_boundary_coil_enabled)
                {
                    constrain_phi_imag_boundary_coil.exec();
                }
            }

            if (!use_frequency_operator_coil_air_only)
            {
                electric_potential_real_relaxation_plate_inner.exec();
                electric_potential_real_relaxation_plate_update.exec(dt_block_step);
                constrain_phi_real_reference_plate.exec();
                if (constrain_phi_boundary_plate_enabled)
                {
                    constrain_phi_real_boundary_plate.exec();
                }

                electric_potential_imag_relaxation_plate_inner.exec();
                electric_potential_imag_relaxation_plate_update.exec(dt_block_step);
                constrain_phi_imag_reference_plate.exec();
                if (constrain_phi_boundary_plate_enabled)
                {
                    constrain_phi_imag_boundary_plate.exec();
                }
            }
        }

        for (size_t k = 0; k < a_relax_iterations; ++k)
        {
            Real air_pre_dt = use_frequency_split_air_sweep ? static_cast<Real>(0.5) * dt_block_step
                                                            : (use_frequency_air_post_sweep ? 0.0 : dt_block_step);
            Real air_post_dt = use_frequency_split_air_sweep ? static_cast<Real>(0.5) * dt_block_step
                                                             : (use_frequency_air_post_sweep ? dt_block_step : 0.0);

            refresh_frequency_coil_scalar_gradients();
            refresh_frequency_coil_magnetic_terms();
            if (use_frequency_coil_magnetic_only)
            {
                if (use_frequency_coil_magnetic_block_solver)
                {
                    if (use_frequency_operator_coil_air_only)
                    {
                        vector_potential_magnetic_only_block_equation_real_coil_air_only.exec(dt_block_step);
                        vector_potential_magnetic_only_block_equation_imag_coil_air_only.exec(dt_block_step);
                    }
                    else
                    {
                        vector_potential_magnetic_only_block_equation_real_coil.exec(dt_block_step);
                        vector_potential_magnetic_only_block_equation_imag_coil.exec(dt_block_step);
                    }
                }
                else
                {
                    if (use_frequency_operator_coil_air_only)
                    {
                        vector_potential_magnetic_only_equation_real_coil_air_only.exec(dt_block_step);
                        vector_potential_magnetic_only_equation_imag_coil_air_only.exec(dt_block_step);
                    }
                    else
                    {
                        vector_potential_magnetic_only_equation_real_coil.exec(dt_block_step);
                        vector_potential_magnetic_only_equation_imag_coil.exec(dt_block_step);
                    }
                }
            }
            else if (use_frequency_coupled_implicit_solver)
            {
                if (use_frequency_operator_coil_air_only)
                {
                    vector_potential_equation_coupled_coil_air_only.exec(dt_block_step);
                }
                else
                {
                    vector_potential_equation_coupled_coil.exec(dt_block_step);
                }
            }
            else
            {
                if (use_frequency_operator_coil_air_only)
                {
                    vector_potential_equation_real_coil_air_only.exec(dt_block_step);
                    vector_potential_equation_imag_coil_air_only.exec(dt_block_step);
                }
                else
                {
                    vector_potential_equation_real_coil.exec(dt_block_step);
                    vector_potential_equation_imag_coil.exec(dt_block_step);
                }
            }

            execute_frequency_air_magnetic_sweep(air_pre_dt);
            if (!use_frequency_operator_coil_air_only)
            {
                refresh_frequency_plate_scalar_gradients();
                refresh_frequency_plate_magnetic_terms();
                if (use_frequency_plate_backtracking)
                {
                    backup_frequency_plate_local_state();

                    if (use_frequency_air_post_sweep)
                    {
                        refresh_frequency_air_magnetic_terms();
                    }
                    FrequencyPlateAcceptanceDiagnostics plate_metric_before =
                        evaluate_plate_acceptance_diagnostics();
                    Real plate_acceptance_metric_before =
                        evaluate_plate_acceptance_metric(plate_metric_before);
                    bool accepted_plate_update = false;
                    Real plate_trial_scale = 1.0;

                    for (size_t trial = 0;
                         trial < frequency_plate_backtracking_max_trials;
                         ++trial)
                    {
                        if (trial > 0)
                        {
                            restore_frequency_plate_local_state();
                            refresh_frequency_plate_scalar_gradients();
                            refresh_frequency_plate_magnetic_terms();
                            if (use_frequency_air_post_sweep)
                            {
                                refresh_frequency_air_magnetic_terms();
                            }
                        }

                        Real trial_dt = dt_block_step * plate_trial_scale;
                        if (use_frequency_coupled_implicit_solver)
                        {
                            vector_potential_equation_coupled_plate.exec(trial_dt);
                        }
                        else
                        {
                            vector_potential_equation_real_plate.exec(trial_dt);
                            vector_potential_equation_imag_plate.exec(trial_dt);
                        }
                        apply_frequency_plate_a_constraints();
                        refresh_frequency_plate_scalar_gradients();
                        refresh_frequency_plate_magnetic_terms();
                        if (use_frequency_air_post_sweep)
                        {
                            refresh_frequency_air_magnetic_terms();
                        }

                        FrequencyPlateAcceptanceDiagnostics plate_metric_after =
                            evaluate_plate_acceptance_diagnostics();
                        Real plate_acceptance_metric_after =
                            evaluate_plate_acceptance_metric(plate_metric_after);
                        bool accept_candidate =
                            std::isfinite(plate_acceptance_metric_after) &&
                            (!std::isfinite(plate_acceptance_metric_before) ||
                             plate_acceptance_metric_after <=
                                 plate_acceptance_metric_before + TinyReal);
                        if (accept_candidate)
                        {
                            accepted_plate_update = true;
                            break;
                        }
                        plate_trial_scale *= static_cast<Real>(0.5);
                    }

                    if (!accepted_plate_update)
                    {
                        restore_frequency_plate_local_state();
                        zero_frequency_plate_change_rate();
                        refresh_frequency_plate_scalar_gradients();
                        refresh_frequency_plate_magnetic_terms();
                        if (use_frequency_air_post_sweep)
                        {
                            refresh_frequency_air_magnetic_terms();
                        }
                    }
                }
                else if (use_frequency_coupled_implicit_solver)
                {
                    vector_potential_equation_coupled_plate.exec(dt_block_step);
                }
                else
                {
                    vector_potential_equation_real_plate.exec(dt_block_step);
                    vector_potential_equation_imag_plate.exec(dt_block_step);
                }

                if (!use_frequency_plate_backtracking &&
                    frequency_plate_block_sweeps > 1)
                {
                    Real plate_acceptance_metric_before =
                        use_frequency_adaptive_plate_block_sweeps
                            ? evaluate_plate_acceptance_metric(
                                  evaluate_plate_acceptance_diagnostics())
                            : std::numeric_limits<Real>::quiet_NaN();
                    for (size_t plate_block_sweep = 1;
                         plate_block_sweep < frequency_plate_block_sweeps;
                         ++plate_block_sweep)
                    {
                        if (use_frequency_adaptive_plate_block_sweeps)
                        {
                            backup_frequency_plate_local_state();
                        }
                        electric_potential_source_real_plate.exec();
                        if (use_frequency_scalar_contact_coupling)
                        {
                            electric_potential_source_real_plate_contact.exec();
                        }
                        electric_potential_source_imag_plate.exec();
                        if (use_frequency_scalar_contact_coupling)
                        {
                            electric_potential_source_imag_plate_contact.exec();
                        }

                        for (size_t phi_sweep = 0; phi_sweep < phi_relax_iterations; ++phi_sweep)
                        {
                            electric_potential_real_relaxation_plate_inner.exec();
                            electric_potential_real_relaxation_plate_update.exec(dt_block_step);
                            constrain_phi_real_reference_plate.exec();
                            if (constrain_phi_boundary_plate_enabled)
                            {
                                constrain_phi_real_boundary_plate.exec();
                            }

                            electric_potential_imag_relaxation_plate_inner.exec();
                            electric_potential_imag_relaxation_plate_update.exec(dt_block_step);
                            constrain_phi_imag_reference_plate.exec();
                            if (constrain_phi_boundary_plate_enabled)
                            {
                                constrain_phi_imag_boundary_plate.exec();
                            }
                        }

                        refresh_frequency_plate_scalar_gradients();
                        refresh_frequency_plate_magnetic_terms();
                        if (use_frequency_coupled_implicit_solver)
                        {
                            vector_potential_equation_coupled_plate.exec(dt_block_step);
                        }
                        else
                        {
                            vector_potential_equation_real_plate.exec(dt_block_step);
                            vector_potential_equation_imag_plate.exec(dt_block_step);
                        }
                        apply_frequency_plate_a_constraints();

                        if (use_frequency_adaptive_plate_block_sweeps)
                        {
                            Real plate_acceptance_metric_after =
                                evaluate_plate_acceptance_metric(
                                    evaluate_plate_acceptance_diagnostics());
                            bool accept_candidate =
                                std::isfinite(plate_acceptance_metric_after) &&
                                (!std::isfinite(plate_acceptance_metric_before) ||
                                 plate_acceptance_metric_after <=
                                     plate_acceptance_metric_before + TinyReal);
                            if (!accept_candidate)
                            {
                                restore_frequency_plate_local_state();
                                refresh_frequency_plate_scalar_gradients();
                                refresh_frequency_plate_magnetic_terms();
                                break;
                            }
                            plate_acceptance_metric_before = plate_acceptance_metric_after;
                        }
                    }
                }
            }

            apply_frequency_a_constraints();
            execute_frequency_air_magnetic_sweep(air_post_dt);
        }
        }
    };

    //------------------------------------------------------------------
    //  Main loop (single-way scaffold)
    //------------------------------------------------------------------
    while (physical_time < end_time)
    {
        if (em_only_diagnostics_mode && iteration >= em_only_max_steps)
        {
            std::cout << "[team7-em-only] reached TEAM7_EM_ONLY_MAX_STEPS="
                      << em_only_max_steps << std::endl;
            break;
        }
        Real remaining_time = end_time - physical_time;
        if (remaining_time <= (em_only_diagnostics_mode ? TinyReal : dt_thermal_min))
        {
            break;
        }
        Real dt_thermal = 0.0;
        if (em_only_diagnostics_mode)
        {
            size_t em_steps_for_macro = SMAX(static_cast<size_t>(1), em_substeps_per_thermal);
            Real dt_macro = dt_em_current * static_cast<Real>(em_steps_for_macro);
            dt_thermal = SMIN(remaining_time, SMAX(dt_em_min, dt_macro));
        }
        else
        {
            dt_thermal = SMIN(dt_thermal_max, SMIN(get_thermal_time_step.exec(), remaining_time));
        }
        if (dt_thermal <= (em_only_diagnostics_mode ? TinyReal : dt_thermal_min))
        {
            break;
        }

        bool stop_em_only_after_frequency_converged = false;
        if (use_frequency_aphi)
        {
            bool perform_frequency_em_update =
                !frequency_em_solve_once || !frequency_em_converged || auto_normalize_source;
            if (perform_frequency_em_update)
            {
                size_t em_iteration_limit = em_substeps_per_thermal;
                size_t em_iteration_min = em_substeps_per_thermal;
                if (use_frequency_pseudo_steady)
                {
                    em_iteration_limit =
                        SMAX(em_substeps_per_thermal, frequency_pseudo_max_iterations);
                    em_iteration_min =
                        SMAX(em_substeps_per_thermal, frequency_pseudo_min_iterations);
                }

                em_iterations_last = 0;
                em_max_rate_last = 0.0;
                em_residual_ratio_last = 0.0;
                Real em_max_delta_a_last = std::numeric_limits<Real>::infinity();
                Real em_joule_rel_change_last = std::numeric_limits<Real>::infinity();
                Real previous_joule_power = std::numeric_limits<Real>::quiet_NaN();
                Real previous_em_max_rate = std::numeric_limits<Real>::quiet_NaN();
                Real previous_em_residual_ratio = std::numeric_limits<Real>::quiet_NaN();
                Real previous_em_plate_residual_metric =
                    std::numeric_limits<Real>::quiet_NaN();
                Real previous_em_conductor_residual_metric =
                    std::numeric_limits<Real>::quiet_NaN();
                Real previous_em_air_residual_metric =
                    std::numeric_limits<Real>::quiet_NaN();
                size_t dt_instability_streak = 0;
                size_t dt_reduction_count = 0;
                size_t saturation_streak = 0;
                size_t hard_guard_reject_streak = 0;
                bool converged_by_rate = false;
                bool converged_by_enhanced = false;
                bool stopped_by_saturation = false;
                bool stopped_by_hard_guard = false;
                for (size_t em_step = 0; em_step < em_iteration_limit; ++em_step)
                {
                    Real dt_em_step = dt_em_current;
                    bool global_backtracking_active =
                        use_frequency_global_backtracking &&
                        em_iterations_last >= em_adaptive_dt_warmup_iterations;
                    if (global_backtracking_active)
                    {
                        refresh_frequency_coil_scalar_gradients();
                        refresh_frequency_coil_magnetic_terms();
                        refresh_frequency_air_scalar_gradients();
                        refresh_frequency_air_magnetic_terms();
                        refresh_frequency_plate_scalar_gradients();
                        refresh_frequency_plate_magnetic_terms();
                        FrequencyEmAcceptanceDiagnostics em_acceptance_before =
                            evaluate_frequency_em_backtracking_diagnostics();
                        backup_frequency_em_state();
                        bool accepted_frequency_substep = false;
                        Real trial_dt_em_step = dt_em_step;
                        Real last_attempt_dt_em_step = dt_em_step;

                        for (size_t trial = 0;
                             trial < frequency_global_backtracking_max_trials;
                             ++trial)
                        {
                            if (trial > 0)
                            {
                                restore_frequency_em_state();
                            }

                            last_attempt_dt_em_step = trial_dt_em_step;
                            execute_frequency_em_substep(trial_dt_em_step);
                            FrequencyEmAcceptanceDiagnostics em_acceptance_after =
                                evaluate_frequency_em_backtracking_diagnostics();
                            bool accept_candidate =
                                std::isfinite(em_acceptance_after.weighted_metric) &&
                                (!std::isfinite(em_acceptance_before.weighted_metric) ||
                                 em_acceptance_after.weighted_metric <=
                                     em_acceptance_before.weighted_metric *
                                         frequency_global_backtracking_metric_growth_limit +
                                         TinyReal);
                            if (accept_candidate)
                            {
                                accepted_frequency_substep = true;
                                dt_em_step = trial_dt_em_step;
                                if (dt_em_step < dt_em_current)
                                {
                                    dt_em_current = SMAX(dt_em_min, dt_em_step);
                                }
                                break;
                            }
                            trial_dt_em_step =
                                SMAX(dt_em_min, trial_dt_em_step * static_cast<Real>(0.5));
                        }

                        if (!accepted_frequency_substep)
                        {
                            restore_frequency_em_state();
                            dt_em_step = last_attempt_dt_em_step;
                            dt_em_current = SMAX(dt_em_min, trial_dt_em_step);
                        }
                    }
                    else
                    {
                        bool hard_guard_active =
                            use_frequency_hard_guard &&
                            em_iterations_last >= frequency_hard_guard_warmup_iterations;
                        if (!hard_guard_active)
                        {
                            execute_frequency_em_substep(dt_em_step);
                        }
                        else
                        {
                            FrequencyEmHardGuardDiagnostics hard_guard_before =
                                evaluate_frequency_em_hard_guard_diagnostics();
                            FrequencyEmHardGuardDiagnostics hard_guard_after_last =
                                hard_guard_before;
                            backup_frequency_em_state();
                            bool accepted_frequency_substep = false;
                            Real trial_dt_em_step = dt_em_step;

                            for (size_t trial = 0;
                                 trial < frequency_hard_guard_max_trials;
                                 ++trial)
                            {
                                if (trial > 0)
                                {
                                    restore_frequency_em_state();
                                }

                                execute_frequency_em_substep(trial_dt_em_step);
                                FrequencyEmHardGuardDiagnostics hard_guard_after =
                                    evaluate_frequency_em_hard_guard_diagnostics();
                                hard_guard_after_last = hard_guard_after;

                                bool metric_ok =
                                    std::isfinite(hard_guard_after.weighted_metric) &&
                                    (!std::isfinite(hard_guard_before.weighted_metric) ||
                                     hard_guard_after.weighted_metric <=
                                         hard_guard_before.weighted_metric *
                                             frequency_hard_guard_metric_growth_limit +
                                             TinyReal);
                                bool probe_ok =
                                    std::isfinite(hard_guard_after.probe_b_metric) &&
                                    (!std::isfinite(hard_guard_before.probe_b_metric) ||
                                     hard_guard_after.probe_b_metric <=
                                         hard_guard_before.probe_b_metric *
                                             frequency_hard_guard_probe_growth_limit +
                                             frequency_hard_guard_probe_abs_tol +
                                             TinyReal);
                                bool probe_ref_ok =
                                    !frequency_hard_guard_use_probe_ref_cap ||
                                    (std::isfinite(hard_guard_after.probe_ref_ratio_metric) &&
                                     hard_guard_after.probe_ref_ratio_metric <=
                                         frequency_hard_guard_probe_ref_ratio_limit);
                                bool accept_candidate =
                                    metric_ok && probe_ok && probe_ref_ok;
                                if (accept_candidate)
                                {
                                    accepted_frequency_substep = true;
                                    dt_em_step = trial_dt_em_step;
                                    if (dt_em_step < dt_em_current)
                                    {
                                        dt_em_current = SMAX(dt_em_min, dt_em_step);
                                    }
                                    hard_guard_reject_streak = 0;
                                    break;
                                }
                                trial_dt_em_step =
                                    SMAX(dt_em_min, trial_dt_em_step *
                                                     frequency_hard_guard_dt_reduce_factor);
                            }

                            if (!accepted_frequency_substep)
                            {
                                restore_frequency_em_state();
                                zero_frequency_em_change_rate();
                                dt_em_current = SMAX(dt_em_min, trial_dt_em_step);
                                std::cout << std::scientific << std::setprecision(6)
                                          << "[team7-em] hard-guard reject substep, "
                                          << "weighted_before=" << hard_guard_before.weighted_metric
                                          << ", weighted_after=" << hard_guard_after_last.weighted_metric
                                          << ", probe_before=" << hard_guard_before.probe_b_metric
                                          << ", probe_after=" << hard_guard_after_last.probe_b_metric
                                          << ", probe_ref_ratio_after="
                                          << hard_guard_after_last.probe_ref_ratio_metric
                                          << ", probe_ref_cap_enabled="
                                          << frequency_hard_guard_use_probe_ref_cap
                                          << ", probe_ref_ratio_limit="
                                          << frequency_hard_guard_probe_ref_ratio_limit
                                          << ", dt_em_current=" << dt_em_current
                                          << std::fixed << std::setprecision(6)
                                          << std::endl;
                                hard_guard_reject_streak++;
                                if (hard_guard_reject_streak >=
                                    frequency_hard_guard_failfast_streak)
                                {
                                    stopped_by_hard_guard = true;
                                    break;
                                }
                                continue;
                            }
                        }
                    }

                    if (em_rate_use_effective_update && em_rate_use_pseudo_dt_scaling)
                    {
                        em_runtime_pseudo_dt_scale =
                            electromagnetics::ComputeNormalizedPseudoTimeStep(dt_em_step, dt_em);
                    }
                    else
                    {
                        em_runtime_pseudo_dt_scale = 1.0;
                    }
                    em_iterations_last++;
                    em_max_rate_last = evaluate_max_em_change_rate();
                    em_residual_ratio_last = evaluate_max_em_residual_ratio();
                    EmResidualBreakdown em_residual_breakdown_last =
                        evaluate_em_residual_breakdown();
                    Real em_plate_residual_metric_last =
                        use_frequency_operator_coil_air_only
                            ? static_cast<Real>(0.0)
                            : em_residual_breakdown_last.plate;
                    Real em_conductor_residual_metric_last =
                        use_frequency_operator_coil_air_only
                            ? em_residual_breakdown_last.coil
                            : SMAX(em_residual_breakdown_last.plate,
                                   em_residual_breakdown_last.coil);
                    Real em_air_residual_metric_last = em_residual_breakdown_last.air;
                    em_max_delta_a_last =
                        use_frequency_aphi
                            ? em_max_rate_last
                            : dt_em_step * em_rate_to_delta_scaling * em_max_rate_last;

                    if (use_frequency_aphi && frequency_failfast_divergence_guard)
                    {
                        bool failfast_residual =
                            std::isfinite(em_residual_ratio_last) &&
                            em_residual_ratio_last > frequency_failfast_residual_ratio_limit;
                        bool failfast_rate =
                            std::isfinite(em_max_rate_last) &&
                            em_max_rate_last > frequency_failfast_em_rate_limit;
                        bool failfast_air_residual =
                            std::isfinite(em_air_residual_metric_last) &&
                            em_air_residual_metric_last > frequency_failfast_air_residual_limit;
                        if (failfast_residual || failfast_rate || failfast_air_residual)
                        {
                            dt_em_current = SMAX(dt_em_min,
                                                 dt_em_current * frequency_failfast_dt_reduce_factor);
                            std::cout << std::scientific << std::setprecision(6)
                                      << "[team7-em] failfast divergence guard triggered, "
                                      << "em_residual_ratio=" << em_residual_ratio_last
                                      << ", em_max_rate=" << em_max_rate_last
                                      << ", em_air_residual=" << em_air_residual_metric_last
                                      << ", dt_em_current=" << dt_em_current
                                      << std::fixed << std::setprecision(6)
                                      << std::endl;
                            stopped_by_hard_guard = true;
                            break;
                        }
                    }

                    if (use_frequency_pseudo_steady && use_enhanced_em_convergence)
                    {
                        if (!use_frequency_operator_coil_air_only)
                        {
                            // Keep Joule diagnostics consistent with the current pseudo-step state.
                            electric_potential_gradient_real_plate_inner.exec();
                            if (use_frequency_scalar_contact_coupling)
                            {
                                electric_potential_gradient_real_plate_contact.exec();
                            }
                            electric_potential_gradient_imag_plate_inner.exec();
                            if (use_frequency_scalar_contact_coupling)
                            {
                                electric_potential_gradient_imag_plate_contact.exec();
                            }
                            electric_field_current_heat_frequency.exec();
                        }

                        PlateDiagnostics em_iteration_diagnostics = evaluate_plate_diagnostics();
                        Real current_joule_power = em_iteration_diagnostics.total_joule_power;
                        if (std::isfinite(previous_joule_power))
                        {
                            Real denominator =
                                SMAX(fabs(current_joule_power), fabs(previous_joule_power));
                            if (denominator > TinyReal)
                            {
                                em_joule_rel_change_last =
                                    fabs(current_joule_power - previous_joule_power) /
                                    (denominator + TinyReal);
                            }
                            else
                            {
                                em_joule_rel_change_last =
                                    fabs(current_joule_power - previous_joule_power);
                            }
                        }
                        previous_joule_power = current_joule_power;
                    }

                    Real effective_em_rate_saturation_limit =
                        em_frequency_rate_saturation_limit;
                    if (em_rate_use_effective_update && em_rate_use_pseudo_dt_scaling)
                    {
                        effective_em_rate_saturation_limit *= em_runtime_pseudo_dt_scale;
                    }
                    bool near_saturation =
                        effective_em_rate_saturation_limit > TinyReal &&
                        em_max_rate_last >=
                            em_saturation_ratio * effective_em_rate_saturation_limit;
                    saturation_streak = near_saturation ? saturation_streak + 1 : 0;

                    bool converged_by_update_rate =
                        em_max_rate_last <= em_rate_convergence_tolerance;
                    bool converged_by_residual =
                        use_em_relative_residual &&
                        em_residual_ratio_last <= em_relative_residual_convergence_tolerance;
                    converged_by_rate = converged_by_update_rate;
                    if (use_frequency_aphi && use_em_relative_residual)
                    {
                        converged_by_rate =
                            converged_by_update_rate && converged_by_residual;
                    }
                    else if (use_em_relative_residual)
                    {
                        converged_by_rate =
                            converged_by_update_rate || converged_by_residual;
                    }
                    if (use_em_adaptive_dt)
                    {
                        bool reduce_dt_for_instability =
                            em_max_rate_last > em_rate_divergence_threshold;
                        if (use_frequency_aphi)
                        {
                            bool allow_growth_based_dt_control =
                                em_iterations_last > em_adaptive_dt_warmup_iterations;
                            Real em_rate_stability_limit =
                                em_rate_convergence_tolerance * em_rate_stability_limit_factor;
                            bool rate_above_stability_limit =
                                allow_growth_based_dt_control &&
                                std::isfinite(em_rate_stability_limit) &&
                                em_max_rate_last > em_rate_stability_limit;
                            bool rate_grew_too_fast =
                                allow_growth_based_dt_control &&
                                std::isfinite(previous_em_max_rate) &&
                                previous_em_max_rate > em_rate_convergence_tolerance &&
                                em_max_rate_last >
                                    previous_em_max_rate * em_rate_growth_trigger;
                            bool residual_grew_too_fast =
                                allow_growth_based_dt_control &&
                                use_em_relative_residual &&
                                std::isfinite(previous_em_residual_ratio) &&
                                previous_em_residual_ratio >
                                    em_relative_residual_convergence_tolerance &&
                                em_residual_ratio_last >
                                    previous_em_residual_ratio * em_residual_growth_trigger;
                            bool growth_with_worsening_residual =
                                rate_grew_too_fast &&
                                use_em_relative_residual &&
                                std::isfinite(previous_em_residual_ratio) &&
                                em_residual_ratio_last > previous_em_residual_ratio;
                            bool air_residual_outpaces_conductor =
                                allow_growth_based_dt_control &&
                                use_em_relative_residual &&
                                std::isfinite(previous_em_air_residual_metric) &&
                                previous_em_air_residual_metric >
                                    em_relative_residual_convergence_tolerance &&
                                em_air_residual_metric_last >
                                    previous_em_air_residual_metric *
                                        em_residual_growth_trigger &&
                                (!std::isfinite(previous_em_conductor_residual_metric) ||
                                 em_conductor_residual_metric_last >=
                                     previous_em_conductor_residual_metric - TinyReal);
                            bool plate_residual_grew_too_fast =
                                !use_frequency_operator_coil_air_only &&
                                allow_growth_based_dt_control &&
                                use_em_relative_residual &&
                                std::isfinite(previous_em_plate_residual_metric) &&
                                previous_em_plate_residual_metric >
                                    em_relative_residual_convergence_tolerance &&
                                em_plate_residual_metric_last >
                                    previous_em_plate_residual_metric *
                                        em_residual_growth_trigger;
                            bool air_residual_outpaces_plate =
                                !use_frequency_operator_coil_air_only &&
                                allow_growth_based_dt_control &&
                                use_em_relative_residual &&
                                std::isfinite(previous_em_air_residual_metric) &&
                                previous_em_air_residual_metric >
                                    em_relative_residual_convergence_tolerance &&
                                em_air_residual_metric_last >
                                    previous_em_air_residual_metric *
                                        em_residual_growth_trigger &&
                                (!std::isfinite(previous_em_plate_residual_metric) ||
                                 em_plate_residual_metric_last >=
                                     previous_em_plate_residual_metric - TinyReal);
                            reduce_dt_for_instability =
                                reduce_dt_for_instability ||
                                rate_above_stability_limit ||
                                near_saturation ||
                                growth_with_worsening_residual ||
                                residual_grew_too_fast ||
                                air_residual_outpaces_conductor ||
                                plate_residual_grew_too_fast ||
                                air_residual_outpaces_plate;
                        }
                        if (reduce_dt_for_instability)
                        {
                            dt_instability_streak++;
                        }
                        else
                        {
                            dt_instability_streak = 0;
                        }
                        bool allow_dt_reduction =
                            reduce_dt_for_instability &&
                            dt_instability_streak >= em_dt_instability_patience &&
                            dt_reduction_count < em_dt_max_reductions_per_solve;
                        if (allow_dt_reduction)
                        {
                            dt_em_current = SMAX(dt_em_min, dt_em_current * em_dt_reduce_factor);
                            dt_reduction_count++;
                            dt_instability_streak = 0;
                        }
                        else if (em_max_rate_last < em_rate_convergence_tolerance)
                        {
                            bool allow_dt_increase = true;
                            Real dt_recovery_threshold = SMAX(
                                em_dt_recovery_min_factor * dt_em_min,
                                em_dt_recovery_relative_threshold * dt_em);
                            bool stuck_at_tiny_dt =
                                dt_em_current <= dt_recovery_threshold + TinyReal;
                            if (use_frequency_aphi && use_em_relative_residual)
                            {
                                allow_dt_increase = converged_by_rate;
                                if (!allow_dt_increase && stuck_at_tiny_dt)
                                {
                                    // Recovery path: avoid getting trapped at dt_em_min
                                    // when updates are already small but residual convergence
                                    // still needs larger corrective steps.
                                    allow_dt_increase = true;
                                }
                            }
                            if (allow_dt_increase)
                            {
                                Real dt_growth_factor = em_dt_increase_factor;
                                if (stuck_at_tiny_dt)
                                {
                                    dt_growth_factor = em_dt_recovery_factor;
                                }
                                dt_em_current = SMIN(dt_em_max, dt_em_current * dt_growth_factor);
                            }
                        }
                    }
                    previous_em_max_rate = em_max_rate_last;
                    previous_em_residual_ratio = em_residual_ratio_last;
                    previous_em_plate_residual_metric = use_frequency_operator_coil_air_only
                                                            ? std::numeric_limits<Real>::quiet_NaN()
                                                            : em_plate_residual_metric_last;
                    previous_em_conductor_residual_metric =
                        em_conductor_residual_metric_last;
                    previous_em_air_residual_metric = em_air_residual_metric_last;
                    if (em_only_diagnostics_mode)
                    {
                        em_iteration_global++;
                        Real em_substep_time =
                            physical_time + static_cast<Real>(em_iterations_last) * dt_em_step;
                        write_em_iteration_diagnostics_row(
                            em_iterations_last, em_substep_time, dt_em_step,
                            em_max_rate_last, em_residual_ratio_last);
                    }

                    if (use_frequency_pseudo_steady && em_iterations_last >= em_iteration_min)
                    {
                        if (converged_by_rate)
                        {
                            break;
                        }
                        if (use_enhanced_em_convergence)
                        {
                            bool delta_converged =
                                em_max_delta_a_last <= em_delta_a_convergence_tolerance;
                            bool joule_converged =
                                std::isfinite(em_joule_rel_change_last) &&
                                em_joule_rel_change_last <= em_joule_rel_convergence_tolerance;
                            converged_by_enhanced = delta_converged && joule_converged;
                            if (converged_by_enhanced)
                            {
                                break;
                            }
                        }
                        if (em_stop_on_saturation &&
                            em_saturation_patience > 0 &&
                            saturation_streak >= em_saturation_patience)
                        {
                            stopped_by_saturation = true;
                            break;
                        }
                    }
                }
                if (use_frequency_pseudo_steady &&
                    em_iterations_last >= em_iteration_limit &&
                    !converged_by_rate &&
                    !converged_by_enhanced &&
                    !stopped_by_saturation)
                {
                    EmRateBreakdown em_rate_breakdown = evaluate_em_rate_breakdown();
                    EmResidualBreakdown em_residual_breakdown = evaluate_em_residual_breakdown();
                    bool reduced_dt = false;
                    if (use_em_adaptive_dt)
                    {
                        if (em_max_rate_last > em_rate_divergence_threshold)
                        {
                            dt_em_current = SMAX(dt_em_min, dt_em_current * em_dt_reduce_factor);
                            reduced_dt = true;
                        }
                    }
                    std::cout << std::scientific << std::setprecision(6)
                              << "[team7-em] pseudo-steady not converged, max_rate="
                              << em_max_rate_last
                              << ", residual_ratio=" << em_residual_ratio_last
                              << ", max_rate_plate=" << em_rate_breakdown.plate
                              << ", max_rate_coil=" << em_rate_breakdown.coil
                              << ", max_rate_air=" << em_rate_breakdown.air
                              << ", residual_plate=" << em_residual_breakdown.plate
                              << ", residual_coil=" << em_residual_breakdown.coil
                              << ", residual_air=" << em_residual_breakdown.air
                              << ", source_scale=" << runtime_source_scale
                              << ", vacuum_reluctivity_scale="
                              << vacuum_reluctivity_scale_runtime
                              << (reduced_dt ? ", reducing dt_em_current to "
                                             : ", keeping dt_em_current at ")
                              << dt_em_current
                              << std::fixed << std::setprecision(6)
                              << std::endl;
                }
                if (stopped_by_saturation)
                {
                    std::cout << std::scientific << std::setprecision(6)
                              << "[team7-em] saturation plateau detected, max_rate="
                              << em_max_rate_last
                              << ", residual_ratio=" << em_residual_ratio_last
                              << ", max_delta_A=" << em_max_delta_a_last
                              << ", joule_rel_change=" << em_joule_rel_change_last
                              << ", saturation_streak=" << saturation_streak
                              << std::fixed << std::setprecision(6)
                              << std::endl;
                }
                if (stopped_by_hard_guard)
                {
                    std::cout << std::scientific << std::setprecision(6)
                              << "[team7-em] hard-guard fail-fast triggered, "
                              << "em_iters_accepted=" << em_iterations_last
                              << ", dt_em_current=" << dt_em_current
                              << ", reject_streak=" << hard_guard_reject_streak
                              << ", source_scale=" << runtime_source_scale
                              << ", vacuum_reluctivity_scale="
                              << vacuum_reluctivity_scale_runtime
                              << std::fixed << std::setprecision(6)
                              << std::endl;
                    stop_due_to_frequency_hard_guard = true;
                }
                bool frequency_converged_candidate =
                    converged_by_rate || converged_by_enhanced || stopped_by_saturation;
                if (stopped_by_hard_guard)
                {
                    frequency_converged_candidate = false;
                }
                if (em_convergence_require_rate)
                {
                    frequency_converged_candidate =
                        frequency_converged_candidate && converged_by_rate;
                }
                if (em_convergence_require_enhanced)
                {
                    frequency_converged_candidate =
                        frequency_converged_candidate && converged_by_enhanced;
                }
                frequency_em_converged = frequency_converged_candidate;
                if (em_only_diagnostics_mode &&
                    frequency_em_solve_once &&
                    frequency_em_converged)
                {
                    stop_em_only_after_frequency_converged = true;
                }
            }
            else
            {
                em_iterations_last = 0;
                em_max_rate_last = 0.0;
                em_residual_ratio_last = 0.0;
            }

            if (!use_frequency_operator_coil_air_only)
            {
                electric_potential_gradient_real_plate_inner.exec();
                if (use_frequency_scalar_contact_coupling)
                {
                    electric_potential_gradient_real_plate_contact.exec();
                }
                electric_potential_gradient_imag_plate_inner.exec();
                if (use_frequency_scalar_contact_coupling)
                {
                    electric_potential_gradient_imag_plate_contact.exec();
                }
                electric_field_current_heat_frequency.exec();
            }
        }
        else
        {
            em_iterations_last = em_substeps_per_thermal;
            for (size_t em_step = 0; em_step < em_substeps_per_thermal; ++em_step)
            {
                if (use_circular_coil_source)
                {
                    set_harmonic_source_on_coil_circular.exec();
                }
                else
                {
                    set_harmonic_source_on_coil.exec();
                }
                if (apply_runtime_source_scale)
                {
                    scale_harmonic_source_on_coil.exec();
                }
                if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 2 && em_step == 0)
                {
                    report_coil_rhs("after_set_source");
                }
                electric_potential_source_coil.exec();
                electric_potential_source_coil_contact.exec();
                electric_potential_source_plate.exec();
                electric_potential_source_plate_contact.exec();
                electric_potential_source_air.exec();
                electric_potential_source_air_contact.exec();

                for (size_t k = 0; k < phi_relax_iterations; ++k)
                {
                    electric_potential_relaxation_air_inner.exec();
                    electric_potential_relaxation_air_contact.exec();
                    electric_potential_relaxation_air_update.exec(dt_em);
                    constrain_phi_reference_air.exec();
                    if (constrain_phi_boundary_air_enabled)
                    {
                        constrain_phi_boundary_air.exec();
                    }
                    if (constrain_phi_all_air_enabled)
                    {
                        constrain_phi_all_air.exec();
                    }

                    electric_potential_relaxation_coil_inner.exec();
                    electric_potential_relaxation_coil_contact.exec();
                    electric_potential_relaxation_coil_update.exec(dt_em);
                    constrain_phi_reference_coil.exec();
                    if (constrain_phi_boundary_coil_enabled)
                    {
                        constrain_phi_boundary_coil.exec();
                    }

                    electric_potential_relaxation_plate_inner.exec();
                    electric_potential_relaxation_plate_contact.exec();
                    electric_potential_relaxation_plate_update.exec(dt_em);
                    constrain_phi_reference_plate.exec();
                    if (constrain_phi_boundary_plate_enabled)
                    {
                        constrain_phi_boundary_plate.exec();
                    }
                }
                if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 2 && em_step == 0)
                {
                    report_coil_nan("after_phi_relaxation");
                    report_plate_air_nan("after_phi_relaxation");
                }

                for (size_t k = 0; k < a_relax_iterations; ++k)
                {
                    electric_potential_gradient_coil_inner.exec();
                    electric_potential_gradient_coil_contact.exec();
                    vector_potential_curl_coil_inner.exec();
                    vector_potential_curl_coil_contact.exec();
                    curl_nu_b_coil_inner.exec();
                    curl_nu_b_coil_contact.exec();
                    vector_potential_equation_coil.exec(dt_em);
                    if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 2 && em_step == 0 && k == 0)
                    {
                        report_coil_rhs("after_first_coil_A_update");
                    }

                    vector_potential_curl_air_inner.exec();
                    vector_potential_curl_air_contact.exec();
                    curl_nu_b_air_inner.exec();
                    curl_nu_b_air_contact.exec();
                    vector_potential_equation_air.exec(dt_em);

                    electric_potential_gradient_plate_inner.exec();
                    electric_potential_gradient_plate_contact.exec();
                    vector_potential_curl_plate_inner.exec();
                    vector_potential_curl_plate_contact.exec();
                    curl_nu_b_plate_inner.exec();
                    curl_nu_b_plate_contact.exec();
                    vector_potential_equation_plate.exec(dt_em);

                    if (constrain_a_reference_coil_enabled)
                    {
                        constrain_a_reference_coil.exec();
                    }
                    if (constrain_a_reference_air_enabled)
                    {
                        constrain_a_reference_air.exec();
                    }
                    if (constrain_a_boundary_air_enabled)
                    {
                        apply_time_domain_air_a_boundary_constraint();
                    }
                    if (constrain_a_reference_plate_enabled)
                    {
                        constrain_a_reference_plate.exec();
                    }
                }
                // Keep A_dot synchronized with the latest A update in each EM substep.
                update_a_dot_coil.exec(dt_em);
                update_a_dot_plate.exec(dt_em);
                update_a_dot_air.exec(dt_em);
                if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 2 && em_step == 0)
                {
                    report_coil_nan("after_a_relaxation");
                    report_plate_air_nan("after_a_relaxation");
                }
                if (em_only_diagnostics_mode)
                {
                    if (em_rate_use_effective_update && em_rate_use_pseudo_dt_scaling)
                    {
                        em_runtime_pseudo_dt_scale =
                            electromagnetics::ComputeNormalizedPseudoTimeStep(dt_em, dt_em);
                    }
                    else
                    {
                        em_runtime_pseudo_dt_scale = 1.0;
                    }
                    em_max_rate_last = evaluate_max_em_change_rate();
                    em_residual_ratio_last = 0.0;
                    em_iteration_global++;
                    Real em_substep_time =
                        physical_time + static_cast<Real>(em_step + 1) * dt_em;
                    write_em_iteration_diagnostics_row(
                        em_step + 1, em_substep_time, dt_em,
                        em_max_rate_last, em_residual_ratio_last);
                }
            }

            if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 3)
            {
                report_plate_air_nan("after_em_substeps");
            }
            electric_potential_gradient_plate_inner.exec();
            electric_potential_gradient_plate_contact.exec();
            electric_field_current_heat.exec();
            if (em_rate_use_effective_update && em_rate_use_pseudo_dt_scaling)
            {
                em_runtime_pseudo_dt_scale =
                    electromagnetics::ComputeNormalizedPseudoTimeStep(dt_em, dt_em);
            }
            else
            {
                em_runtime_pseudo_dt_scale = 1.0;
            }
            em_max_rate_last = evaluate_max_em_change_rate();
            em_residual_ratio_last = 0.0;
            if (debug_coil_nan_detection && !use_frequency_aphi && iteration < 3)
            {
                report_plate_air_nan("after_plate_heat_update");
            }
        }

        if (stop_due_to_frequency_hard_guard)
        {
            std::cout << "[team7-em] stopping outer time loop due to hard-guard fail-fast."
                      << std::endl;
            break;
        }

        if (use_frequency_aphi)
        {
            bool allow_vacuum_reluctivity_feedback_update =
                std::isfinite(em_residual_ratio_last) &&
                em_residual_ratio_last <=
                    vacuum_reluctivity_feedback_max_residual_ratio;
            if (allow_vacuum_reluctivity_feedback_update)
            {
                try_update_vacuum_reluctivity_scale_by_coil_air_probe(iteration + 1);
            }
            try_update_runtime_source_scale_by_coil_air_probe(iteration + 1);
            CoilSourceDiagnostics coil_source_diagnostics_for_control =
                evaluate_coil_source_diagnostics();
            try_update_runtime_source_scale_by_circuit(
                iteration + 1, coil_source_diagnostics_for_control);
            try_auto_normalize_runtime_source_scale(iteration + 1);
        }

        if (!em_only_diagnostics_mode)
        {
            thermal_relaxation.exec(dt_thermal);
            add_joule_heat_to_temperature.exec(dt_thermal);
        }

        physical_time += dt_thermal;
        iteration += 1;

        if (physical_time >= next_output_time - TinyReal)
        {
            Real harmonic_factor =
                sin(2.0 * Pi * harmonic_frequency_hz_runtime * physical_time + harmonic_phase_runtime);
            PlateDiagnostics diagnostics = evaluate_plate_diagnostics();
            CoilSourceDiagnostics coil_source_diagnostics = evaluate_coil_source_diagnostics();
            EmTermBreakdown em_term_breakdown = evaluate_em_term_breakdown();
            EmInterfaceShellBreakdown em_interface_breakdown = evaluate_em_interface_shell_breakdown();
            EmMagneticDiagonalBreakdown em_magnetic_diagonal_breakdown =
                evaluate_em_magnetic_diagonal_breakdown();
            std::cout << std::fixed << std::setprecision(6)
                      << "N=" << iteration
                      << ", t=" << physical_time
                      << ", dt_thermal=" << dt_thermal
                      << ", dt_em_current=" << dt_em_current
                      << ", em_iters=" << em_iterations_last
                      << ", max_em_rate=" << std::scientific << em_max_rate_last
                      << ", residual_ratio=" << em_residual_ratio_last
                      << std::fixed
                      << ", sin(omega t)=" << harmonic_factor
                      << ", source_scale=" << std::scientific << runtime_source_scale
                      << ", vacuum_reluctivity_scale="
                      << vacuum_reluctivity_scale_runtime
                      << std::fixed
                      << ", P_joule_si=" << std::scientific << diagnostics.total_joule_power
                      << std::fixed
                      << std::endl;
            EmRateBreakdown em_rate_breakdown = evaluate_em_rate_breakdown();
            EmResidualBreakdown em_residual_breakdown = evaluate_em_residual_breakdown();
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-em-rate] plate=" << em_rate_breakdown.plate
                      << ", coil=" << em_rate_breakdown.coil
                      << ", air=" << em_rate_breakdown.air
                      << ", residual_plate=" << em_residual_breakdown.plate
                      << ", residual_coil=" << em_residual_breakdown.coil
                      << ", residual_air=" << em_residual_breakdown.air
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-em-source] coil_NI_cfg_real=" << coil_source_diagnostics.configured_ampere_turns_real
                      << ", coil_NI_cfg_imag=" << coil_source_diagnostics.configured_ampere_turns_imag
                      << ", coil_NI_eq_real=" << coil_source_diagnostics.equivalent_ampere_turns_real
                      << ", coil_NI_eq_imag=" << coil_source_diagnostics.equivalent_ampere_turns_imag
                      << ", coil_avg_|Js|=" << coil_source_diagnostics.avg_source_density_magnitude
                      << ", coil_max_|Js|=" << coil_source_diagnostics.max_source_density_magnitude
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-em-terms] plate(avg curl=" << em_term_breakdown.plate.avg_curl_nu_b
                      << ", avg sigmaGradPhi=" << em_term_breakdown.plate.avg_sigma_grad_phi
                      << ", avg omegaSigmaA=" << em_term_breakdown.plate.avg_omega_sigma_a
                      << ", avg residual=" << em_term_breakdown.plate.avg_residual
                      << "), coil(avg source=" << em_term_breakdown.coil.avg_source
                      << ", avg curl=" << em_term_breakdown.coil.avg_curl_nu_b
                      << ", avg sigmaGradPhi=" << em_term_breakdown.coil.avg_sigma_grad_phi
                      << ", avg omegaSigmaA=" << em_term_breakdown.coil.avg_omega_sigma_a
                      << ", avg residual=" << em_term_breakdown.coil.avg_residual
                      << "), air(avg curl=" << em_term_breakdown.air.avg_curl_nu_b
                      << ", avg residual=" << em_term_breakdown.air.avg_residual
                      << ")"
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-em-interface] plate(inner curl=" << em_interface_breakdown.plate_inner.avg_curl_nu_b
                      << ", inner residual=" << em_interface_breakdown.plate_inner.avg_residual
                      << ", air curl=" << em_interface_breakdown.plate_air.avg_curl_nu_b
                      << ", air residual=" << em_interface_breakdown.plate_air.avg_residual
                      << "), coil(inner curl=" << em_interface_breakdown.coil_inner.avg_curl_nu_b
                      << ", inner residual=" << em_interface_breakdown.coil_inner.avg_residual
                      << ", air curl=" << em_interface_breakdown.coil_air.avg_curl_nu_b
                      << ", air residual=" << em_interface_breakdown.coil_air.avg_residual
                      << ")"
                      << std::fixed << std::setprecision(6)
                      << std::endl;
            if (enable_coil_air_interface_b_diag)
            {
                EmCoilAirInterfaceShellFieldDiagnostics coil_air_ab = evaluate_em_coil_air_interface_ab_diagnostics();
                const Real b_shell_ratio =
                    (coil_air_ab.coil_inner_avg_b_mag > TinyReal)
                        ? (coil_air_ab.coil_air_avg_b_mag / coil_air_ab.coil_inner_avg_b_mag)
                        : 0.0;
                const Real a_shell_ratio =
                    (coil_air_ab.coil_inner_avg_a_mag > TinyReal)
                        ? (coil_air_ab.coil_air_avg_a_mag / coil_air_ab.coil_inner_avg_a_mag)
                        : 0.0;
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-em-interface-ab] coil_inner(avg_|A|=" << coil_air_ab.coil_inner_avg_a_mag
                          << ", avg_|B|=" << coil_air_ab.coil_inner_avg_b_mag
                          << ", vol=" << coil_air_ab.coil_inner_volume
                          << "), air_near_coil(avg_|A|=" << coil_air_ab.coil_air_avg_a_mag
                          << ", avg_|B|=" << coil_air_ab.coil_air_avg_b_mag
                          << ", vol=" << coil_air_ab.coil_air_volume
                          << ", B_air/B_coil_inner=" << b_shell_ratio
                          << ", A_air/A_coil_inner=" << a_shell_ratio << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
            if (enable_contact_magnetic_diagonal_diagnostics)
            {
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-em-magdiag] plate(local=" << em_magnetic_diagonal_breakdown.plate.avg_local_diagonal
                          << ", contact=" << em_magnetic_diagonal_breakdown.plate.avg_contact_diagonal
                          << ", contact_ratio=" << em_magnetic_diagonal_breakdown.plate.avg_contact_ratio
                          << ", conservative/jacobi=" << em_magnetic_diagonal_breakdown.plate.avg_conservative_to_jacobi
                          << ", balanced/jacobi=" << em_magnetic_diagonal_breakdown.plate.avg_balanced_to_jacobi
                          << "), coil(local=" << em_magnetic_diagonal_breakdown.coil.avg_local_diagonal
                          << ", contact=" << em_magnetic_diagonal_breakdown.coil.avg_contact_diagonal
                          << ", contact_ratio=" << em_magnetic_diagonal_breakdown.coil.avg_contact_ratio
                          << ", conservative/jacobi=" << em_magnetic_diagonal_breakdown.coil.avg_conservative_to_jacobi
                          << ", balanced/jacobi=" << em_magnetic_diagonal_breakdown.coil.avg_balanced_to_jacobi
                          << "), air(local=" << em_magnetic_diagonal_breakdown.air.avg_local_diagonal
                          << ", contact=" << em_magnetic_diagonal_breakdown.air.avg_contact_diagonal
                          << ", contact_ratio=" << em_magnetic_diagonal_breakdown.air.avg_contact_ratio
                          << ", conservative/jacobi=" << em_magnetic_diagonal_breakdown.air.avg_conservative_to_jacobi
                          << ", balanced/jacobi=" << em_magnetic_diagonal_breakdown.air.avg_balanced_to_jacobi
                          << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-em-magdiag-interface] plate(inner contact_ratio="
                          << em_magnetic_diagonal_breakdown.plate_inner.avg_contact_ratio
                          << ", air contact_ratio="
                          << em_magnetic_diagonal_breakdown.plate_air.avg_contact_ratio
                          << "), coil(inner contact_ratio="
                          << em_magnetic_diagonal_breakdown.coil_inner.avg_contact_ratio
                          << ", air contact_ratio="
                          << em_magnetic_diagonal_breakdown.coil_air.avg_contact_ratio
                          << ")"
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
            plate_diagnostics_file << physical_time << ","
                                   << diagnostics.total_joule_power << ","
                                   << diagnostics.avg_joule_heat << ","
                                   << diagnostics.max_joule_heat << ","
                                   << diagnostics.avg_current_density_magnitude << ","
                                   << diagnostics.max_current_density_magnitude << ","
                                   << diagnostics.avg_temperature << ","
                                   << diagnostics.max_temperature << ","
                                   << diagnostics.max_a_rate << ","
                                   << diagnostics.max_a_dot << ","
                                   << diagnostics.max_grad_phi << ","
                                   << diagnostics.max_electric_field << ","
                                   << em_iterations_last << ","
                                   << em_max_rate_last << ","
                                   << em_residual_ratio_last << ","
                                   << dt_em_current << ","
                                   << diagnostics.total_volume << ","
                                   << em_term_breakdown.plate.avg_residual << ","
                                   << em_term_breakdown.coil.avg_residual << ","
                                   << em_term_breakdown.air.avg_residual << ","
                                   << em_interface_breakdown.plate_inner.avg_residual << ","
                                   << em_interface_breakdown.plate_air.avg_residual << ","
                                   << em_interface_breakdown.coil_inner.avg_residual << ","
                                   << em_interface_breakdown.coil_air.avg_residual << ","
                                   << coil_source_diagnostics.equivalent_ampere_turns_real << ","
                                   << coil_source_diagnostics.equivalent_ampere_turns_imag << "\n";
            Vec3d b_air_above = sample_magnetic_flux_density(air_vector_potential_curl_probe, probe_b_air_above_plate);
            Vec3d b_air_below = sample_magnetic_flux_density(air_vector_potential_curl_probe, probe_b_air_below_plate);
            Vec3d b_plate_center = sample_magnetic_flux_density(plate_vector_potential_curl_probe, probe_b_plate_center);
            Real b_air_above_mag =
                sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_above_plate);
            Real b_air_below_mag =
                sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_below_plate);
            Real b_plate_center_mag =
                sample_magnetic_flux_density_magnitude(plate_vector_potential_curl_probe, plate_vector_potential_curl_imag, probe_b_plate_center);
            Real t_plate_center = plate_temperature[probe_t_plate_center.particle_index];
            Real t_plate_edge_xplus = plate_temperature[probe_t_plate_edge_xplus.particle_index];
            team7_observables_file << physical_time << ","
                                   << diagnostics.total_joule_power << ","
                                   << diagnostics.avg_temperature << ","
                                   << diagnostics.max_temperature << ","
                                   << diagnostics.max_a_rate << ","
                                   << diagnostics.max_a_dot << ","
                                   << diagnostics.max_grad_phi << ","
                                   << diagnostics.max_electric_field << ","
                                   << em_iterations_last << ","
                                   << em_max_rate_last << ","
                                   << em_residual_ratio_last << ","
                                   << dt_em_current << ","
                                   << b_air_above[0] << ","
                                   << b_air_above[1] << ","
                                   << b_air_above[2] << ","
                                   << b_air_above_mag << ","
                                   << b_air_below[0] << ","
                                   << b_air_below[1] << ","
                                   << b_air_below[2] << ","
                                   << b_air_below_mag << ","
                                   << b_plate_center[0] << ","
                                   << b_plate_center[1] << ","
                                   << b_plate_center[2] << ","
                                   << b_plate_center_mag << ","
                                   << t_plate_center << ","
                                   << t_plate_edge_xplus << ","
                                   << em_term_breakdown.plate.avg_residual << ","
                                   << em_term_breakdown.coil.avg_residual << ","
                                   << em_term_breakdown.air.avg_residual << ","
                                   << em_interface_breakdown.plate_inner.avg_residual << ","
                                   << em_interface_breakdown.plate_air.avg_residual << ","
                                   << em_interface_breakdown.coil_inner.avg_residual << ","
                                   << em_interface_breakdown.coil_air.avg_residual << ","
                                   << coil_source_diagnostics.equivalent_ampere_turns_real << ","
                                   << coil_source_diagnostics.equivalent_ampere_turns_imag << "\n";
            plate_diagnostics_file.flush();
            team7_observables_file.flush();
            if (em_only_diagnostics_mode)
            {
                em_iteration_diagnostics_file.flush();
            }
            const bool is_last_step =
                (physical_time + TinyReal >= end_time) || (iteration >= output_steps);
            if (!record_last_step_only || is_last_step)
            {
                write_states.writeToFile();
                wrote_last_step_recording = is_last_step;
            }
            next_output_time += output_interval;
        }
        if (stop_em_only_after_frequency_converged)
        {
            std::cout << "[team7-em-only] frequency EM solve marked converged/stagnated; "
                         "stopping EM-only diagnostics loop."
                      << std::endl;
            break;
        }
    }

    if (record_last_step_only && !wrote_last_step_recording)
    {
        write_states.writeToFile();
    }

    if (enable_validation_output)
    {
        PlateDiagnostics final_diagnostics = evaluate_plate_diagnostics();
        Real b_air_above_final_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_above_plate);
        Real b_air_below_final_mag =
            sample_magnetic_flux_density_magnitude(air_vector_potential_curl_probe, air_vector_potential_curl_imag, probe_b_air_below_plate);
        Real b_plate_center_final_mag =
            sample_magnetic_flux_density_magnitude(plate_vector_potential_curl_probe, plate_vector_potential_curl_imag, probe_b_plate_center);
        Real t_plate_center_final = plate_temperature[probe_t_plate_center.particle_index];
        Real t_plate_edge_xplus_final = plate_temperature[probe_t_plate_edge_xplus.particle_index];
        Real delta_t_plate_center_final = t_plate_center_final - t_plate_center_initial;
        Real delta_t_plate_edge_xplus_final = t_plate_edge_xplus_final - t_plate_edge_xplus_initial;

        std::ofstream validation_file(
            io_environment.OutputFolder() + "/team7_validation.csv",
            std::ios::out | std::ios::trunc);
        validation_file << std::setprecision(12);
        validation_file << "metric,simulated,reference,abs_error,rel_error,pass\n";

        size_t total_checked_metrics = 0;
        size_t passed_metrics = 0;
        auto write_validation_metric = [&](const std::string &metric_name, Real simulated_value, Real reference_value)
        {
            validation_file << metric_name << "," << simulated_value << ",";
            if (!std::isfinite(reference_value))
            {
                validation_file << "NA,NA,NA,NA\n";
                return;
            }

            Real absolute_error = fabs(simulated_value - reference_value);
            Real reference_norm = fabs(reference_value);
            Real relative_error = std::numeric_limits<Real>::infinity();
            bool pass_abs = absolute_error <= validation_abs_tolerance;
            bool pass_rel = false;
            if (reference_norm > TinyReal)
            {
                relative_error = absolute_error / reference_norm;
                pass_rel = relative_error <= validation_rel_tolerance;
            }
            else if (pass_abs)
            {
                relative_error = 0.0;
            }
            bool pass = pass_abs || pass_rel;

            total_checked_metrics++;
            if (pass)
            {
                passed_metrics++;
            }

            validation_file << reference_value << ","
                            << absolute_error << ",";
            if (std::isfinite(relative_error))
            {
                validation_file << relative_error << ",";
            }
            else
            {
                validation_file << "INF,";
            }
            validation_file << (pass ? 1 : 0) << "\n";
        };

        write_validation_metric("total_joule_power", final_diagnostics.total_joule_power, ref_total_joule_power);
        write_validation_metric("avg_temperature", final_diagnostics.avg_temperature, ref_avg_temperature);
        write_validation_metric("max_temperature", final_diagnostics.max_temperature, ref_max_temperature);
        write_validation_metric("B_air_above_plate_mag", b_air_above_final_mag, ref_b_air_above_mag);
        write_validation_metric("B_air_below_plate_mag", b_air_below_final_mag, ref_b_air_below_mag);
        write_validation_metric("B_plate_center_mag", b_plate_center_final_mag, ref_b_plate_center_mag);
        write_validation_metric("T_plate_center", t_plate_center_final, ref_t_plate_center);
        write_validation_metric("T_plate_edge_xplus", t_plate_edge_xplus_final, ref_t_plate_edge_xplus);
        write_validation_metric("delta_T_plate_center", delta_t_plate_center_final, ref_delta_t_plate_center);
        write_validation_metric("delta_T_plate_edge_xplus", delta_t_plate_edge_xplus_final, ref_delta_t_plate_edge_xplus);
        validation_file.flush();

        std::cout << std::scientific << std::setprecision(6)
                  << "[team7-validation] checked_metrics=" << total_checked_metrics
                  << ", passed_metrics=" << passed_metrics
                  << ", pass_rate=" << (total_checked_metrics > 0
                                            ? static_cast<Real>(passed_metrics) /
                                                  static_cast<Real>(total_checked_metrics)
                                            : 0.0)
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    }

    if (enable_curve_validation && (!b_curve_sample_points.empty() || !j_curve_sample_points.empty()))
    {
        std::ofstream curve_validation_file(
            io_environment.OutputFolder() + "/team7_curve_validation.csv",
            std::ios::out | std::ios::trunc);
        curve_validation_file << std::setprecision(12);
        curve_validation_file
            << "quantity,point_id,x_mm,target_x,target_y,target_z,sampled_x,sampled_y,sampled_z,"
            << "complex_mode,simulated,ref_50,ref_200,ref_selected,abs_error,rel_error,"
            << "support_point_count,support_weight_sum,support_nearest_distance,support_weighted_mean_distance,"
            << "boundary_clearance_x,boundary_clearance_y,boundary_clearance_z,high_quality,pass\n";

        auto evaluate_curve_samples =
            [&](std::vector<Team7CurveSamplePoint> &sample_points,
                bool is_bz_curve,
                ComplexCurveComponentMode complex_mode,
                bool high_quality_only) -> Team7CurveValidationSummary
        {
            Team7CurveValidationSummary summary;
            summary.quantity_name = is_bz_curve ? "Bz" : "Jy";
            summary.complex_mode_name =
                complex_curve_component_mode_name(complex_mode) +
                std::string(high_quality_only ? "_high_quality" : "_all");
            summary.total_points = sample_points.size();

            Real sum_abs_error = 0.0;
            Real sum_squared_error = 0.0;
            Real sum_relative_error = 0.0;
            size_t relative_error_count = 0;
            Real sum_abs_reference = 0.0;
            Real sum_abs_simulated = 0.0;

            for (size_t i = 0; i != sample_points.size(); ++i)
            {
                Team7CurveSamplePoint &sample = sample_points[i];
                if (!high_quality_only)
                {
                    assign_curve_sample_quality(
                        sample,
                        is_bz_curve,
                        is_bz_curve ? dp_air_finest : dp_plate,
                        is_bz_curve ? b_curve_interpolation_radius : j_curve_interpolation_radius);
                }
                Real simulated_value = 0.0;
                if (is_bz_curve)
                {
                    if (use_frequency_aphi && air_vector_potential_curl_real != nullptr)
                    {
                        simulated_value =
                            interpolate_complex_curve_component(air_positions,
                                                                air_vol,
                                                                air_particles.TotalRealParticles(),
                                                                air_vector_potential_curl_real,
                                                                air_vector_potential_curl_imag,
                                                                sample,
                                                                b_curve_component_index,
                                                                b_curve_interpolation_radius,
                                                                complex_mode);
                    }
                    else
                    {
                        simulated_value =
                            interpolate_complex_curve_component(air_positions,
                                                                air_vol,
                                                                air_particles.TotalRealParticles(),
                                                                air_vector_potential_curl_probe,
                                                                nullptr,
                                                                sample,
                                                                b_curve_component_index,
                                                                b_curve_interpolation_radius,
                                                                complex_mode);
                    }
                    simulated_value *= b_curve_unit_scale;
                }
                else
                {
                    if (use_frequency_aphi && plate_current_density_real != nullptr)
                    {
                        simulated_value =
                            interpolate_complex_curve_component(plate_positions,
                                                                plate_vol,
                                                                plate_particles.TotalRealParticles(),
                                                                plate_current_density_real,
                                                                plate_current_density_imag,
                                                                sample,
                                                                j_curve_component_index,
                                                                j_curve_interpolation_radius,
                                                                complex_mode);
                    }
                    else
                    {
                        simulated_value =
                            interpolate_complex_curve_component(plate_positions,
                                                                plate_vol,
                                                                plate_particles.TotalRealParticles(),
                                                                plate_current_density,
                                                                nullptr,
                                                                sample,
                                                                j_curve_component_index,
                                                                j_curve_interpolation_radius,
                                                                complex_mode);
                    }
                    simulated_value *= j_curve_unit_scale;
                }

                sample.simulated_value = simulated_value;
                sample.reference_selected =
                    std::isfinite(sample.reference_selected)
                        ? sample.reference_selected
                        : select_reference_value_by_frequency(
                              Team7CurveReferencePoint{sample.x_mm, sample.reference_50, sample.reference_200},
                              curve_reference_frequency_hz);

                if (std::isfinite(sample.reference_selected))
                {
                    if (high_quality_only && !sample.high_quality)
                    {
                        continue;
                    }
                    Real absolute_error = fabs(sample.simulated_value - sample.reference_selected);
                    Real reference_norm = fabs(sample.reference_selected);
                    Real relative_error = std::numeric_limits<Real>::infinity();
                    bool pass_abs = absolute_error <= curve_validation_abs_tolerance;
                    bool pass_rel = false;
                    if (reference_norm > TinyReal)
                    {
                        relative_error = absolute_error / reference_norm;
                        pass_rel = relative_error <= curve_validation_rel_tolerance;
                    }
                    else if (pass_abs)
                    {
                        relative_error = 0.0;
                    }
                    bool pass = pass_abs || pass_rel;

                    sample.absolute_error = absolute_error;
                    sample.relative_error = relative_error;
                    sample.pass = pass;

                    summary.checked_points++;
                    if (pass)
                    {
                        summary.passed_points++;
                    }
                    sum_abs_error += absolute_error;
                    sum_squared_error += absolute_error * absolute_error;
                    summary.max_abs_error = SMAX(summary.max_abs_error, absolute_error);
                    sum_abs_reference += reference_norm;
                    sum_abs_simulated += fabs(sample.simulated_value);
                    if (std::isfinite(relative_error))
                    {
                        sum_relative_error += relative_error;
                        relative_error_count++;
                    }
                }
                else
                {
                    sample.absolute_error = std::numeric_limits<Real>::quiet_NaN();
                    sample.relative_error = std::numeric_limits<Real>::quiet_NaN();
                    sample.pass = false;
                }

                curve_validation_file << sample.quantity_name << ","
                                      << i << ","
                                      << sample.x_mm << ","
                                      << sample.target_position[0] << ","
                                      << sample.target_position[1] << ","
                                      << sample.target_position[2] << ","
                                      << sample.sampled_position[0] << ","
                                      << sample.sampled_position[1] << ","
                                      << sample.sampled_position[2] << ","
                                      << summary.complex_mode_name << ","
                                      << sample.simulated_value << ","
                                      << sample.reference_50 << ","
                                      << sample.reference_200 << ","
                                      << sample.reference_selected << ","
                                      << sample.absolute_error << ","
                                      << sample.relative_error << ","
                                      << sample.support_point_count << ","
                                      << sample.support_weight_sum << ","
                                      << sample.support_nearest_distance << ","
                                      << sample.support_weighted_mean_distance << ","
                                      << sample.boundary_clearance_x << ","
                                      << sample.boundary_clearance_y << ","
                                      << sample.boundary_clearance_z << ","
                                      << (sample.high_quality ? 1 : 0) << ","
                                      << (sample.pass ? 1 : 0) << "\n";
            }

            if (summary.checked_points > 0)
            {
                summary.mean_abs_error = sum_abs_error / static_cast<Real>(summary.checked_points);
                summary.rmse = sqrt(sum_squared_error / static_cast<Real>(summary.checked_points));
                summary.mean_rel_error = relative_error_count > 0
                                             ? sum_relative_error / static_cast<Real>(relative_error_count)
                                             : 0.0;
                summary.mean_abs_reference =
                    sum_abs_reference / static_cast<Real>(summary.checked_points);
                summary.mean_abs_simulated =
                    sum_abs_simulated / static_cast<Real>(summary.checked_points);
                summary.mean_abs_ref_over_mean_abs_sim =
                    summary.mean_abs_reference / (summary.mean_abs_simulated + TinyReal);
            }
            return summary;
        };

        Team7CurveValidationSummary b_summary =
            evaluate_curve_samples(b_curve_sample_points, true, b_curve_complex_mode, false);
        Team7CurveValidationSummary j_summary =
            evaluate_curve_samples(j_curve_sample_points, false, j_curve_complex_mode, false);
        curve_validation_file.flush();

        std::ofstream curve_validation_summary_file(
            io_environment.OutputFolder() + "/team7_curve_validation_summary.csv",
            std::ios::out | std::ios::trunc);
        curve_validation_summary_file << std::setprecision(12);
        curve_validation_summary_file
            << "quantity,complex_mode,total_points,checked_points,passed_points,pass_rate,"
            << "mean_abs_error,rmse,max_abs_error,mean_rel_error,"
            << "mean_abs_reference,mean_abs_simulated,mean_abs_ref_over_mean_abs_sim\n";

        auto write_curve_summary_row = [&](const Team7CurveValidationSummary &summary)
        {
            Real pass_rate = summary.checked_points > 0
                                 ? static_cast<Real>(summary.passed_points) /
                                       static_cast<Real>(summary.checked_points)
                                 : 0.0;
            curve_validation_summary_file << summary.quantity_name << ","
                                          << summary.complex_mode_name << ","
                                          << summary.total_points << ","
                                          << summary.checked_points << ","
                                          << summary.passed_points << ","
                                          << pass_rate << ","
                                          << summary.mean_abs_error << ","
                                          << summary.rmse << ","
                                          << summary.max_abs_error << ","
                                          << summary.mean_rel_error << ","
                                          << summary.mean_abs_reference << ","
                                          << summary.mean_abs_simulated << ","
                                          << summary.mean_abs_ref_over_mean_abs_sim << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-curve-validation] quantity=" << summary.quantity_name
                      << ", complex_mode=" << summary.complex_mode_name
                      << ", checked=" << summary.checked_points
                      << ", passed=" << summary.passed_points
                      << ", mean_abs_error=" << summary.mean_abs_error
                      << ", rmse=" << summary.rmse
                      << ", max_abs_error=" << summary.max_abs_error
                      << ", mean_rel_error=" << summary.mean_rel_error
                      << ", mean_abs_ref=" << summary.mean_abs_reference
                      << ", mean_abs_sim=" << summary.mean_abs_simulated
                      << ", ref_over_sim=" << summary.mean_abs_ref_over_mean_abs_sim
                      << std::fixed << std::setprecision(6)
                      << std::endl;
        };
        write_curve_summary_row(b_summary);
        write_curve_summary_row(j_summary);
        curve_validation_summary_file.flush();

        std::ofstream curve_validation_quality_summary_file(
            io_environment.OutputFolder() + "/team7_curve_validation_summary_high_quality.csv",
            std::ios::out | std::ios::trunc);
        curve_validation_quality_summary_file << std::setprecision(12);
        curve_validation_quality_summary_file
            << "quantity,complex_mode,total_points,checked_points,passed_points,pass_rate,"
            << "mean_abs_error,rmse,max_abs_error,mean_rel_error,"
            << "mean_abs_reference,mean_abs_simulated,mean_abs_ref_over_mean_abs_sim\n";
        auto write_curve_quality_summary_row = [&](const Team7CurveValidationSummary &summary)
        {
            Real pass_rate = summary.checked_points > 0
                                 ? static_cast<Real>(summary.passed_points) /
                                       static_cast<Real>(summary.checked_points)
                                 : 0.0;
            curve_validation_quality_summary_file
                << summary.quantity_name << ","
                << summary.complex_mode_name << ","
                << summary.total_points << ","
                << summary.checked_points << ","
                << summary.passed_points << ","
                << pass_rate << ","
                << summary.mean_abs_error << ","
                << summary.rmse << ","
                << summary.max_abs_error << ","
                << summary.mean_rel_error << ","
                << summary.mean_abs_reference << ","
                << summary.mean_abs_simulated << ","
                << summary.mean_abs_ref_over_mean_abs_sim << "\n";
        };
        Team7CurveValidationSummary b_summary_high_quality =
            evaluate_curve_samples(b_curve_sample_points, true, b_curve_complex_mode, true);
        Team7CurveValidationSummary j_summary_high_quality =
            evaluate_curve_samples(j_curve_sample_points, false, j_curve_complex_mode, true);
        write_curve_quality_summary_row(b_summary_high_quality);
        write_curve_quality_summary_row(j_summary_high_quality);
        curve_validation_quality_summary_file.flush();
        if (enable_j_curve_decomposition_output &&
            use_frequency_aphi &&
            !j_curve_sample_points.empty() &&
            plate_vector_potential_real != nullptr &&
            plate_vector_potential_imag != nullptr &&
            plate_electric_potential_gradient_real != nullptr &&
            plate_electric_potential_gradient_imag != nullptr &&
            plate_current_density_real != nullptr &&
            plate_current_density_imag != nullptr &&
            plate_sigma != nullptr)
        {
            std::ofstream j_curve_decomposition_file(
                io_environment.OutputFolder() + "/team7_j_curve_decomposition.csv",
                std::ios::out | std::ios::trunc);
            j_curve_decomposition_file << std::setprecision(12);
            j_curve_decomposition_file
                << "point_id,x_mm,target_x,target_y,target_z,sampled_x,sampled_y,sampled_z,"
                << "support_point_count,support_weight_sum,support_nearest_distance,support_weighted_mean_distance,"
                << "boundary_clearance_x,boundary_clearance_y,boundary_clearance_z,high_quality,"
                << "component,complex_mode,unit_scale,sigma,"
                << "a_real_comp,a_imag_comp,grad_phi_real_comp,grad_phi_imag_comp,"
                << "e_real_comp,e_imag_comp,j_real_comp,j_imag_comp,"
                << "j_from_jomega_a_real,j_from_jomega_a_imag,"
                << "j_from_minus_grad_phi_real,j_from_minus_grad_phi_imag,"
                << "j_simulated_projected,j_from_jomega_a_projected,j_from_minus_grad_phi_projected,"
                << "j_reconstruction_error_real,j_reconstruction_error_imag\n";

            for (size_t i = 0; i != j_curve_sample_points.size(); ++i)
            {
                const Team7CurveSamplePoint &sample = j_curve_sample_points[i];
                Real sigma_interp = interpolate_scalar_curve_value(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_sigma,
                    sample,
                    j_curve_interpolation_radius);

                Real a_real_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_vector_potential_real,
                    plate_vector_potential_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::RealPart);
                Real a_imag_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_vector_potential_real,
                    plate_vector_potential_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::ImagPart);
                Real grad_phi_real_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_electric_potential_gradient_real,
                    plate_electric_potential_gradient_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::RealPart);
                Real grad_phi_imag_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_electric_potential_gradient_real,
                    plate_electric_potential_gradient_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::ImagPart);
                Real j_real_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_current_density_real,
                    plate_current_density_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::RealPart);
                Real j_imag_comp = interpolate_complex_curve_component(
                    plate_positions,
                    plate_vol,
                    plate_particles.TotalRealParticles(),
                    plate_current_density_real,
                    plate_current_density_imag,
                    sample,
                    j_curve_component_index,
                    j_curve_interpolation_radius,
                    ComplexCurveComponentMode::ImagPart);

                Real e_real_comp = harmonic_angular_frequency * a_imag_comp - grad_phi_real_comp;
                Real e_imag_comp = -harmonic_angular_frequency * a_real_comp - grad_phi_imag_comp;
                Real j_from_jomega_a_real = sigma_interp * harmonic_angular_frequency * a_imag_comp;
                Real j_from_jomega_a_imag = -sigma_interp * harmonic_angular_frequency * a_real_comp;
                Real j_from_minus_grad_phi_real = -sigma_interp * grad_phi_real_comp;
                Real j_from_minus_grad_phi_imag = -sigma_interp * grad_phi_imag_comp;
                Real j_simulated_projected = project_complex_component(
                    j_real_comp, j_imag_comp, j_curve_complex_mode) * j_curve_unit_scale;
                Real j_from_jomega_a_projected = project_complex_component(
                    j_from_jomega_a_real, j_from_jomega_a_imag, j_curve_complex_mode) *
                    j_curve_unit_scale;
                Real j_from_minus_grad_phi_projected = project_complex_component(
                    j_from_minus_grad_phi_real, j_from_minus_grad_phi_imag, j_curve_complex_mode) *
                    j_curve_unit_scale;
                Real j_reconstruction_error_real =
                    j_real_comp - (j_from_jomega_a_real + j_from_minus_grad_phi_real);
                Real j_reconstruction_error_imag =
                    j_imag_comp - (j_from_jomega_a_imag + j_from_minus_grad_phi_imag);

                j_curve_decomposition_file
                    << i << ","
                    << sample.x_mm << ","
                    << sample.target_position[0] << ","
                    << sample.target_position[1] << ","
                    << sample.target_position[2] << ","
                    << sample.sampled_position[0] << ","
                    << sample.sampled_position[1] << ","
                    << sample.sampled_position[2] << ","
                    << sample.support_point_count << ","
                    << sample.support_weight_sum << ","
                    << sample.support_nearest_distance << ","
                    << sample.support_weighted_mean_distance << ","
                    << sample.boundary_clearance_x << ","
                    << sample.boundary_clearance_y << ","
                    << sample.boundary_clearance_z << ","
                    << (sample.high_quality ? 1 : 0) << ","
                    << j_curve_component_index << ","
                    << complex_curve_component_mode_name(j_curve_complex_mode) << ","
                    << j_curve_unit_scale << ","
                    << sigma_interp << ","
                    << a_real_comp << ","
                    << a_imag_comp << ","
                    << grad_phi_real_comp << ","
                    << grad_phi_imag_comp << ","
                    << e_real_comp << ","
                    << e_imag_comp << ","
                    << j_real_comp << ","
                    << j_imag_comp << ","
                    << j_from_jomega_a_real << ","
                    << j_from_jomega_a_imag << ","
                    << j_from_minus_grad_phi_real << ","
                    << j_from_minus_grad_phi_imag << ","
                    << j_simulated_projected << ","
                    << j_from_jomega_a_projected << ","
                    << j_from_minus_grad_phi_projected << ","
                    << j_reconstruction_error_real << ","
                    << j_reconstruction_error_imag << "\n";
            }
            j_curve_decomposition_file.flush();
        }
    }
    if (enable_coil_air_biot_savart_validation)
    {
        std::ofstream coil_air_biot_savart_validation_file(
            io_environment.OutputFolder() + "/team7_coil_air_biot_savart_validation.csv",
            std::ios::out | std::ios::trunc);
        coil_air_biot_savart_validation_file << std::setprecision(12);
        coil_air_biot_savart_validation_file
            << "quantity,point_id,target_x,target_y,target_z,complex_mode,"
            << "simulated,reference,abs_error,rel_error,pass\n";

        auto evaluate_scalar_error =
            [&](Real simulated_value, Real reference_value) -> std::tuple<Real, Real, bool>
        {
            Real absolute_error = fabs(simulated_value - reference_value);
            Real reference_norm = fabs(reference_value);
            Real relative_error = std::numeric_limits<Real>::infinity();
            bool pass_abs = absolute_error <= coil_air_biot_savart_abs_tolerance;
            bool pass_rel = false;
            if (reference_norm > TinyReal)
            {
                relative_error = absolute_error / reference_norm;
                pass_rel = relative_error <= coil_air_biot_savart_rel_tolerance;
            }
            else if (pass_abs)
            {
                relative_error = 0.0;
            }
            return {absolute_error, relative_error, pass_abs || pass_rel};
        };

        ScalarValidationSummary probe_magnitude_summary;
        probe_magnitude_summary.quantity_name = "B_probe_magnitude";
        ScalarValidationSummary bz_curve_summary;
        bz_curve_summary.quantity_name = "Bz_curve_vs_biot_savart";
        size_t probe_passed_points = 0;
        size_t bz_curve_passed_points = 0;
        Real probe_sum_abs_error = 0.0;
        Real probe_sum_squared_error = 0.0;
        Real probe_sum_relative_error = 0.0;
        size_t probe_relative_error_count = 0;
        Real bz_curve_sum_abs_error = 0.0;
        Real bz_curve_sum_squared_error = 0.0;
        Real bz_curve_sum_relative_error = 0.0;
        size_t bz_curve_relative_error_count = 0;

        auto evaluate_probe_entry = [&](const ProbeDescriptor &probe, size_t point_id)
        {
            Vec3d simulated_real = sample_magnetic_flux_density(air_vector_potential_curl_probe, probe);
            Vec3d simulated_imag = Vec3d::Zero();
            if (use_frequency_aphi && air_vector_potential_curl_imag != nullptr)
            {
                simulated_imag = sample_magnetic_flux_density(air_vector_potential_curl_imag, probe);
            }
            Real simulated_magnitude = complex_vector_magnitude(simulated_real, simulated_imag);

            std::pair<Vec3d, Vec3d> biot_savart_field =
                compute_biot_savart_magnetic_flux_density(probe.sampled_position);
            Real biot_savart_magnitude =
                complex_vector_magnitude(biot_savart_field.first, biot_savart_field.second);

            auto [absolute_error, relative_error, pass] =
                evaluate_scalar_error(simulated_magnitude, biot_savart_magnitude);
            probe_magnitude_summary.checked_points++;
            if (pass)
            {
                probe_passed_points++;
            }
            probe_sum_abs_error += absolute_error;
            probe_sum_squared_error += absolute_error * absolute_error;
            probe_magnitude_summary.max_abs_error =
                SMAX(probe_magnitude_summary.max_abs_error, absolute_error);
            if (std::isfinite(relative_error))
            {
                probe_sum_relative_error += relative_error;
                probe_relative_error_count++;
            }

            coil_air_biot_savart_validation_file
                << probe.name << ","
                << point_id << ","
                << probe.sampled_position[0] << ","
                << probe.sampled_position[1] << ","
                << probe.sampled_position[2] << ","
                << "magnitude,"
                << simulated_magnitude << ","
                << biot_savart_magnitude << ","
                << absolute_error << ","
                << relative_error << ","
                << (pass ? 1 : 0) << "\n";
        };
        evaluate_probe_entry(probe_b_air_above_plate, 0);
        evaluate_probe_entry(probe_b_air_below_plate, 1);

        for (size_t i = 0; i != b_curve_sample_points.size(); ++i)
        {
            const Team7CurveSamplePoint &sample = b_curve_sample_points[i];
            Real simulated_value = 0.0;
            if (use_frequency_aphi && air_vector_potential_curl_real != nullptr)
            {
                simulated_value =
                    interpolate_complex_curve_component(air_positions,
                                                        air_vol,
                                                        air_particles.TotalRealParticles(),
                                                        air_vector_potential_curl_real,
                                                        air_vector_potential_curl_imag,
                                                        sample,
                                                        b_curve_component_index,
                                                        b_curve_interpolation_radius,
                                                        b_curve_complex_mode);
            }
            else
            {
                simulated_value =
                    interpolate_complex_curve_component(air_positions,
                                                        air_vol,
                                                        air_particles.TotalRealParticles(),
                                                        air_vector_potential_curl_probe,
                                                        nullptr,
                                                        sample,
                                                        b_curve_component_index,
                                                        b_curve_interpolation_radius,
                                                        b_curve_complex_mode);
            }
            simulated_value *= b_curve_unit_scale;

            std::pair<Vec3d, Vec3d> biot_savart_field =
                compute_biot_savart_magnetic_flux_density(sample.target_position);
            Real biot_savart_value = project_complex_component(
                                         biot_savart_field.first[b_curve_component_index],
                                         biot_savart_field.second[b_curve_component_index],
                                         b_curve_complex_mode) *
                                     b_curve_unit_scale;

            auto [absolute_error, relative_error, pass] =
                evaluate_scalar_error(simulated_value, biot_savart_value);
            bz_curve_summary.checked_points++;
            if (pass)
            {
                bz_curve_passed_points++;
            }
            bz_curve_sum_abs_error += absolute_error;
            bz_curve_sum_squared_error += absolute_error * absolute_error;
            bz_curve_summary.max_abs_error =
                SMAX(bz_curve_summary.max_abs_error, absolute_error);
            if (std::isfinite(relative_error))
            {
                bz_curve_sum_relative_error += relative_error;
                bz_curve_relative_error_count++;
            }

            coil_air_biot_savart_validation_file
                << "Bz_curve" << ","
                << i << ","
                << sample.target_position[0] << ","
                << sample.target_position[1] << ","
                << sample.target_position[2] << ","
                << complex_curve_component_mode_name(b_curve_complex_mode) << ","
                << simulated_value << ","
                << biot_savart_value << ","
                << absolute_error << ","
                << relative_error << ","
                << (pass ? 1 : 0) << "\n";
        }
        coil_air_biot_savart_validation_file.flush();

        auto finalize_scalar_summary =
            [&](ScalarValidationSummary &summary,
                Real sum_abs_error,
                Real sum_squared_error,
                Real sum_relative_error,
                size_t relative_error_count)
        {
            if (summary.checked_points == 0)
            {
                return;
            }
            Real checked_count = static_cast<Real>(summary.checked_points);
            summary.mean_abs_error = sum_abs_error / checked_count;
            summary.rmse = sqrt(sum_squared_error / checked_count);
            summary.mean_rel_error = relative_error_count > 0
                                         ? sum_relative_error / static_cast<Real>(relative_error_count)
                                         : 0.0;
        };
        finalize_scalar_summary(probe_magnitude_summary,
                                probe_sum_abs_error,
                                probe_sum_squared_error,
                                probe_sum_relative_error,
                                probe_relative_error_count);
        finalize_scalar_summary(bz_curve_summary,
                                bz_curve_sum_abs_error,
                                bz_curve_sum_squared_error,
                                bz_curve_sum_relative_error,
                                bz_curve_relative_error_count);

        std::ofstream coil_air_biot_savart_summary_file(
            io_environment.OutputFolder() + "/team7_coil_air_biot_savart_summary.csv",
            std::ios::out | std::ios::trunc);
        coil_air_biot_savart_summary_file << std::setprecision(12);
        coil_air_biot_savart_summary_file
            << "quantity,checked_points,passed_points,pass_rate,mean_abs_error,rmse,max_abs_error,mean_rel_error\n";

        auto write_scalar_summary_row =
            [&](const ScalarValidationSummary &summary, size_t passed_points)
        {
            Real pass_rate = summary.checked_points > 0
                                 ? static_cast<Real>(passed_points) /
                                       static_cast<Real>(summary.checked_points)
                                 : 0.0;
            coil_air_biot_savart_summary_file
                << summary.quantity_name << ","
                << summary.checked_points << ","
                << passed_points << ","
                << pass_rate << ","
                << summary.mean_abs_error << ","
                << summary.rmse << ","
                << summary.max_abs_error << ","
                << summary.mean_rel_error << "\n";
            std::cout << std::scientific << std::setprecision(6)
                      << "[team7-coil-air-biot-savart] quantity=" << summary.quantity_name
                      << ", checked=" << summary.checked_points
                      << ", passed=" << passed_points
                      << ", mean_abs_error=" << summary.mean_abs_error
                      << ", rmse=" << summary.rmse
                      << ", max_abs_error=" << summary.max_abs_error
                      << ", mean_rel_error=" << summary.mean_rel_error
                      << std::fixed << std::setprecision(6)
                      << std::endl;
        };
        write_scalar_summary_row(probe_magnitude_summary, probe_passed_points);
        write_scalar_summary_row(bz_curve_summary, bz_curve_passed_points);
        coil_air_biot_savart_summary_file.flush();

        if (enable_coil_air_probe_path_diagnostics)
        {
            struct ProbePathLineSpec
            {
                std::string line_name;
                Vec3d probe_point = Vec3d::Zero();
            };
            struct ProbePathLineSummary
            {
                std::string line_name;
                ScalarValidationSummary b_magnitude;
                ScalarValidationSummary b_z;
                size_t b_magnitude_passed = 0;
                size_t b_z_passed = 0;
            };
            auto signed_distance_to_coil_box = [&](const Vec3d &position) -> Real
            {
                Vec3d distance = (position - coil_center).cwiseAbs() - coil_halfsize;
                if (distance.maxCoeff() <= static_cast<Real>(0.0))
                {
                    Vec3d margin = coil_halfsize - (position - coil_center).cwiseAbs();
                    return -margin.minCoeff();
                }
                return point_box_distance(position, coil_center, coil_halfsize);
            };
            auto build_path_start_point = [&](const Vec3d &probe_point) -> Vec3d
            {
                Vec3d direction = probe_point - coil_center;
                Real direction_norm = direction.norm();
                if (direction_norm <= TinyReal)
                {
                    direction = Vec3d(0.0, 0.0, 1.0);
                    direction_norm = 1.0;
                }
                direction /= direction_norm;
                Real ray_to_surface = std::numeric_limits<Real>::infinity();
                for (size_t axis = 0; axis != 3; ++axis)
                {
                    Real direction_component = fabs(direction[axis]);
                    if (direction_component <= TinyReal)
                    {
                        continue;
                    }
                    Real candidate = coil_halfsize[axis] / direction_component;
                    ray_to_surface = SMIN(ray_to_surface, candidate);
                }
                if (!std::isfinite(ray_to_surface))
                {
                    ray_to_surface = coil_halfsize.maxCoeff();
                }
                return coil_center +
                       direction * (ray_to_surface + coil_air_probe_path_start_offset);
            };
            auto finalize_path_summary =
                [&](ScalarValidationSummary &summary,
                    Real sum_abs_error,
                    Real sum_squared_error,
                    Real sum_relative_error,
                    size_t relative_error_count)
            {
                if (summary.checked_points == 0)
                {
                    return;
                }
                Real checked_count = static_cast<Real>(summary.checked_points);
                summary.mean_abs_error = sum_abs_error / checked_count;
                summary.rmse = sqrt(sum_squared_error / checked_count);
                summary.mean_rel_error = relative_error_count > 0
                                             ? sum_relative_error /
                                                   static_cast<Real>(relative_error_count)
                                             : 0.0;
            };

            std::vector<ProbePathLineSpec> path_lines = {
                {"probe_path_air_above", probe_b_air_above_plate.sampled_position},
                {"probe_path_air_below", probe_b_air_below_plate.sampled_position}};

            std::ofstream probe_path_profile_file(
                io_environment.OutputFolder() + "/team7_coil_air_probe_path_profile.csv",
                std::ios::out | std::ios::trunc);
            probe_path_profile_file << std::setprecision(12);
            probe_path_profile_file
                << "line_name,point_id,path_s,target_x,target_y,target_z,"
                << "nearest_air_particle,sampled_x,sampled_y,sampled_z,"
                << "signed_distance_to_coil_box,"
                << "sim_b_mag,ref_b_mag,abs_err_b_mag,rel_err_b_mag,pass_b_mag,"
                << "sim_bz,ref_bz,abs_err_bz,rel_err_bz,pass_bz\n";

            std::vector<ProbePathLineSummary> path_line_summaries;
            path_line_summaries.reserve(path_lines.size());

            for (const ProbePathLineSpec &path_line : path_lines)
            {
                ProbePathLineSummary line_summary;
                line_summary.line_name = path_line.line_name;
                line_summary.b_magnitude.quantity_name =
                    path_line.line_name + "_Bmag";
                line_summary.b_z.quantity_name = path_line.line_name + "_Bz";
                Real b_magnitude_sum_abs_error = 0.0;
                Real b_magnitude_sum_squared_error = 0.0;
                Real b_magnitude_sum_relative_error = 0.0;
                size_t b_magnitude_relative_error_count = 0;
                Real b_z_sum_abs_error = 0.0;
                Real b_z_sum_squared_error = 0.0;
                Real b_z_sum_relative_error = 0.0;
                size_t b_z_relative_error_count = 0;

                Vec3d path_start_point = build_path_start_point(path_line.probe_point);
                Vec3d path_end_point = path_line.probe_point;

                for (size_t point_id = 0; point_id != coil_air_probe_path_points; ++point_id)
                {
                    Real path_s = coil_air_probe_path_points > 1
                                      ? static_cast<Real>(point_id) /
                                            static_cast<Real>(coil_air_probe_path_points - 1)
                                      : 0.0;
                    Vec3d target_position =
                        path_start_point + path_s * (path_end_point - path_start_point);
                    size_t nearest_air_particle =
                        find_nearest_particle_index(air_particles, air_positions, target_position);
                    Vec3d sampled_position = air_positions[nearest_air_particle];
                    Real signed_distance =
                        signed_distance_to_coil_box(target_position);

                    std::pair<Vec3d, Vec3d> simulated_field =
                        interpolate_complex_air_curl_at_position(
                            target_position, coil_air_probe_path_interp_radius);
                    std::pair<Vec3d, Vec3d> reference_field =
                        compute_biot_savart_magnetic_flux_density(target_position);

                    Real simulated_b_magnitude = complex_vector_magnitude(
                        simulated_field.first, simulated_field.second);
                    Real reference_b_magnitude = complex_vector_magnitude(
                        reference_field.first, reference_field.second);
                    auto [b_mag_abs_error, b_mag_rel_error, b_mag_pass] =
                        evaluate_scalar_error(simulated_b_magnitude, reference_b_magnitude);

                    Real simulated_b_z = project_complex_component(
                        simulated_field.first[2], simulated_field.second[2],
                        b_curve_complex_mode);
                    Real reference_b_z = project_complex_component(
                        reference_field.first[2], reference_field.second[2],
                        b_curve_complex_mode);
                    auto [b_z_abs_error, b_z_rel_error, b_z_pass] =
                        evaluate_scalar_error(simulated_b_z, reference_b_z);

                    line_summary.b_magnitude.checked_points++;
                    line_summary.b_z.checked_points++;
                    if (b_mag_pass)
                    {
                        line_summary.b_magnitude_passed++;
                    }
                    if (b_z_pass)
                    {
                        line_summary.b_z_passed++;
                    }
                    b_magnitude_sum_abs_error += b_mag_abs_error;
                    b_magnitude_sum_squared_error +=
                        b_mag_abs_error * b_mag_abs_error;
                    line_summary.b_magnitude.max_abs_error =
                        SMAX(line_summary.b_magnitude.max_abs_error, b_mag_abs_error);
                    if (std::isfinite(b_mag_rel_error))
                    {
                        b_magnitude_sum_relative_error += b_mag_rel_error;
                        b_magnitude_relative_error_count++;
                    }
                    b_z_sum_abs_error += b_z_abs_error;
                    b_z_sum_squared_error += b_z_abs_error * b_z_abs_error;
                    line_summary.b_z.max_abs_error =
                        SMAX(line_summary.b_z.max_abs_error, b_z_abs_error);
                    if (std::isfinite(b_z_rel_error))
                    {
                        b_z_sum_relative_error += b_z_rel_error;
                        b_z_relative_error_count++;
                    }

                    probe_path_profile_file
                        << path_line.line_name << ","
                        << point_id << ","
                        << path_s << ","
                        << target_position[0] << ","
                        << target_position[1] << ","
                        << target_position[2] << ","
                        << nearest_air_particle << ","
                        << sampled_position[0] << ","
                        << sampled_position[1] << ","
                        << sampled_position[2] << ","
                        << signed_distance << ","
                        << simulated_b_magnitude << ","
                        << reference_b_magnitude << ","
                        << b_mag_abs_error << ","
                        << b_mag_rel_error << ","
                        << (b_mag_pass ? 1 : 0) << ","
                        << simulated_b_z << ","
                        << reference_b_z << ","
                        << b_z_abs_error << ","
                        << b_z_rel_error << ","
                        << (b_z_pass ? 1 : 0) << "\n";
                }

                finalize_path_summary(line_summary.b_magnitude,
                                      b_magnitude_sum_abs_error,
                                      b_magnitude_sum_squared_error,
                                      b_magnitude_sum_relative_error,
                                      b_magnitude_relative_error_count);
                finalize_path_summary(line_summary.b_z,
                                      b_z_sum_abs_error,
                                      b_z_sum_squared_error,
                                      b_z_sum_relative_error,
                                      b_z_relative_error_count);
                path_line_summaries.push_back(line_summary);
            }
            probe_path_profile_file.flush();

            std::ofstream probe_path_summary_file(
                io_environment.OutputFolder() + "/team7_coil_air_probe_path_summary.csv",
                std::ios::out | std::ios::trunc);
            probe_path_summary_file << std::setprecision(12);
            probe_path_summary_file
                << "line_name,metric,checked_points,passed_points,pass_rate,"
                << "mean_abs_error,rmse,max_abs_error,mean_rel_error\n";
            for (const ProbePathLineSummary &line_summary : path_line_summaries)
            {
                Real b_magnitude_pass_rate =
                    line_summary.b_magnitude.checked_points > 0
                        ? static_cast<Real>(line_summary.b_magnitude_passed) /
                              static_cast<Real>(line_summary.b_magnitude.checked_points)
                        : 0.0;
                Real b_z_pass_rate =
                    line_summary.b_z.checked_points > 0
                        ? static_cast<Real>(line_summary.b_z_passed) /
                              static_cast<Real>(line_summary.b_z.checked_points)
                        : 0.0;
                probe_path_summary_file
                    << line_summary.line_name << ","
                    << "B_magnitude"
                    << ","
                    << line_summary.b_magnitude.checked_points << ","
                    << line_summary.b_magnitude_passed << ","
                    << b_magnitude_pass_rate << ","
                    << line_summary.b_magnitude.mean_abs_error << ","
                    << line_summary.b_magnitude.rmse << ","
                    << line_summary.b_magnitude.max_abs_error << ","
                    << line_summary.b_magnitude.mean_rel_error << "\n";
                probe_path_summary_file
                    << line_summary.line_name << ","
                    << "Bz"
                    << ","
                    << line_summary.b_z.checked_points << ","
                    << line_summary.b_z_passed << ","
                    << b_z_pass_rate << ","
                    << line_summary.b_z.mean_abs_error << ","
                    << line_summary.b_z.rmse << ","
                    << line_summary.b_z.max_abs_error << ","
                    << line_summary.b_z.mean_rel_error << "\n";
                std::cout << std::scientific << std::setprecision(6)
                          << "[team7-coil-air-probe-path] line="
                          << line_summary.line_name
                          << ", bmag_mean_rel_error="
                          << line_summary.b_magnitude.mean_rel_error
                          << ", bz_mean_rel_error="
                          << line_summary.b_z.mean_rel_error
                          << std::fixed << std::setprecision(6)
                          << std::endl;
            }
            probe_path_summary_file.flush();
        }
    }
    write_coil_air_transfer_diagnostics();
    if (em_only_diagnostics_mode)
    {
        em_iteration_diagnostics_file.flush();
    }

    std::cout << "TEAM7 single-way harmonic scaffold finished." << std::endl;
    return 0;
}
