/**
 * @file test_3d_em_aphi_laplace_eigen_team7_stl_em_only.cpp
 * @brief Formal TEAM7 geometry (coil.stl + plate.stl) with Laplace-structured
 * A-phi EM-only solve on a small air box. The first version intentionally keeps
 * the air domain small to limit particle count for SparseLU.
 */

#include "sphinxsys.h"
#include "electromagnetic_aphi_laplace_eigen.hpp"
#include "em_adaptive_cell_linked_list.h"
#include "electromagnetic_multiturn_coil_drive.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <string>
#include <chrono>

using namespace SPH;

namespace
{
using Complex = electromagnetics::Complex;
namespace fs = std::filesystem;
using extra_electromagnetics::NormalizeVectorOrDefault;
using extra_electromagnetics::TangentialDirectionAroundAxis;

std::string resolve_existing_path(const std::string &preferred,
                                  const StdVec<std::string> &fallbacks)
{
    if (fs::exists(preferred))
    {
        return preferred;
    }
    for (const std::string &candidate : fallbacks)
    {
        if (fs::exists(candidate))
        {
            return candidate;
        }
    }
    return preferred;
}

const std::string path_coil_stl = resolve_existing_path(
    "./input/coil.stl",
    {
        "./tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_laplace_eigen_team7_stl_em_only/bin/input/coil.stl",
        "../tests/extra_source_and_tests/3d_examples/particle_generation_em/data/coil.stl",
        "../build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_laplace_eigen_team7_stl_em_only/bin/input/coil.stl",
    });
const std::string path_plate_stl = resolve_existing_path(
    "./input/plate.stl",
    {
        "./tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_laplace_eigen_team7_stl_em_only/bin/input/plate.stl",
        "../tests/extra_source_and_tests/3d_examples/particle_generation_em/data/plate.stl",
        "../build/tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_laplace_eigen_team7_stl_em_only/bin/input/plate.stl",
    });

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

bool parse_bool_env_pre(const char *name, bool default_value)
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

Real parse_real_env_pre(const char *name, Real default_value)
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

struct Team7AirBoxParams
{
    Vec3d air_lower{-200.0, -200.0, -200.0};
    Vec3d air_upper{500.0, 500.0, 300.0};
    bool used_small_air_box = false;
    bool used_tight_stl_air_box = false;
    Real tight_stl_margin = 40.0;
};

Team7AirBoxParams parse_team7_air_box_pre_system()
{
    Team7AirBoxParams params;
    params.tight_stl_margin = parse_real_env_pre("TEAM7_TIGHT_STL_AIR_MARGIN", params.tight_stl_margin);
    if (parse_bool_env_pre("TEAM7_USE_TIGHT_STL_AIR_BOX", false))
    {
        params.used_tight_stl_air_box = true;
    }
    if (parse_bool_env_pre("TEAM7_USE_SMALL_AIR_BOX", true))
    {
        params.used_small_air_box = true;
        params.air_lower = Vec3d(-30.0, -80.0, -40.0);
        params.air_upper = Vec3d(340.0, 360.0, 180.0);
    }
    params.air_lower[0] = parse_real_env_pre("TEAM7_AIR_BOX_LOWER_X", params.air_lower[0]);
    params.air_lower[1] = parse_real_env_pre("TEAM7_AIR_BOX_LOWER_Y", params.air_lower[1]);
    params.air_lower[2] = parse_real_env_pre("TEAM7_AIR_BOX_LOWER_Z", params.air_lower[2]);
    params.air_upper[0] = parse_real_env_pre("TEAM7_AIR_BOX_UPPER_X", params.air_upper[0]);
    params.air_upper[1] = parse_real_env_pre("TEAM7_AIR_BOX_UPPER_Y", params.air_upper[1]);
    params.air_upper[2] = parse_real_env_pre("TEAM7_AIR_BOX_UPPER_Z", params.air_upper[2]);
    return params;
}

BoundingBoxd combine_bounds(const BoundingBoxd &a, const BoundingBoxd &b)
{
    return BoundingBoxd(a.lower_.cwiseMin(b.lower_), a.upper_.cwiseMax(b.upper_));
}

Team7AirBoxParams resolve_team7_air_box_from_geometry(const Team7AirBoxParams &requested,
                                                      const BoundingBoxd &coil_bounds,
                                                      const BoundingBoxd &plate_bounds)
{
    Team7AirBoxParams resolved = requested;
    if (requested.used_tight_stl_air_box)
    {
        BoundingBoxd union_bounds = combine_bounds(coil_bounds, plate_bounds);
        Vec3d margin = Vec3d::Ones() * requested.tight_stl_margin;
        resolved.air_lower = union_bounds.lower_ - margin;
        resolved.air_upper = union_bounds.upper_ + margin;
    }
    return resolved;
}

const Team7AirBoxParams requested_team7_air_box = parse_team7_air_box_pre_system();
const Real dp_0 = get_env_real_local("EM_APHI_LAPLACE_DP",
                                     get_env_real_local("TEAM7_SPH_REF_DP", 20.0));
const Real boundary_width = get_env_real_local("EM_APHI_LAPLACE_BOUNDARY_WIDTH", 2.0 * dp_0);
const Real boundary_shell_thickness =
    get_env_real_local("EM_APHI_LAPLACE_BOUNDARY_SHELL_THICKNESS", 1.25 * dp_0);

const Real sigma_air =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SIGMA_AIR", 1.0e-8);
const Real sigma_plate =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SIGMA_PLATE", 3.526e7);
const Real sigma_coil =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SIGMA_COIL", 1.0e-8);
const Real nu_air =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_NU_AIR", 1.0);
const Real nu_plate =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_NU_PLATE", 1.0);
const Real nu_coil =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_NU_COIL", 1.0);
const Real frequency_hz =
    get_env_real_local("EM_APHI_LAPLACE_FREQUENCY_HZ", 50.0);
const Real coil_current_rms =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_CURRENT_RMS", 0.1);
const Real coil_turns =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_TURNS", 2742.0);
const Real coil_current_peak_factor =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_CURRENT_PEAK_FACTOR", std::sqrt(2.0));
const Real coil_effective_area_input_geom =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_EFFECTIVE_AREA", -1.0);
const Real geom_length_to_m =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_GEOM_LENGTH_TO_M", 1.0e-3);
const bool use_circular_coil_source =
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_USE_CIRCULAR_COIL_SOURCE", true);
const Vecd coil_source_axis(
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_SOURCE_AXIS_X", 0.0),
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_SOURCE_AXIS_Y", 0.0),
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_COIL_SOURCE_AXIS_Z", 1.0));
const Vecd fallback_source_direction(
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SOURCE_DIR_X", 0.0),
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SOURCE_DIR_Y", 0.0),
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SOURCE_DIR_Z", 1.0));

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
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_WRITE_VTP", true);
const bool use_three_body_layout =
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_USE_THREE_BODY_LAYOUT", false);
const bool use_multi_resolution =
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_USE_MULTI_RESOLUTION", false);
const Real dp_coil_input =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_DP_COIL", dp_0);
const Real dp_plate_input =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_DP_PLATE", dp_0);
const Real dp_air_coarse_input =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_DP_AIR_COARSE", 2.0 * dp_0);
const int air_refinement_levels_input = static_cast<int>(
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_AIR_REFINEMENT_LEVELS", 1.0));
const bool use_magnetic_only_air =
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_USE_MAGNETIC_ONLY_AIR", true);
const bool use_magnetic_only_coil =
    get_env_bool_local("EM_APHI_LAPLACE_TEAM7_EMONLY_USE_MAGNETIC_ONLY_COIL", true);
const Real sigma_solve_regularization =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_SIGMA_SOLVE_REGULARIZATION", 1.0e-6);
const Real contact_gradient_scale =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_CONTACT_GRADIENT_SCALE", 1.0);
const Real contact_diffusion_scale =
    get_env_real_local("EM_APHI_LAPLACE_TEAM7_EMONLY_CONTACT_DIFFUSION_SCALE", 1.0);

class Team7AirBoxShape : public ComplexShape
{
  public:
    explicit Team7AirBoxShape(const std::string &shape_name,
                              const Vec3d &air_lower,
                              const Vec3d &air_upper)
        : ComplexShape(shape_name)
    {
        Vec3d halfsize = 0.5 * (air_upper - air_lower);
        Vec3d center = 0.5 * (air_upper + air_lower);
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class CoilShape : public ComplexShape
{
  public:
    explicit CoilShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
    }
};

class PlateShape : public ComplexShape
{
  public:
    explicit PlateShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
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
        add<GeometricShapeBox>(Transform(center), halfsize, "OuterBoundary");
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

class InnerCoilPlateUnionShape : public ComplexShape
{
  public:
    explicit InnerCoilPlateUnionShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
    }
};

class AdaptiveAroundInnerStructure
    : public extra_em_mesh::AdaptiveWithTolerantCellLinkedList<AdaptiveWithinShape>
{
    Shape *inner_shape_;

  public:
    AdaptiveAroundInnerStructure(Real global_resolution, Real h_spacing_ratio,
                                 Real refinement_to_global, int local_refinement_level,
                                 Shape *inner_shape)
        : extra_em_mesh::AdaptiveWithTolerantCellLinkedList<AdaptiveWithinShape>(
              global_resolution, h_spacing_ratio,
              refinement_to_global, local_refinement_level),
          inner_shape_(inner_shape) {}

    Real getLocalSpacing(Shape &shape, const Vecd &position) override
    {
        Real phi = inner_shape_->findSignedDistance(position);
        return phi < 0.0 ? finest_spacing_bound_ : smoothedSpacing(phi, 2.0 * spacing_ref_);
    }
};

bool is_boundary_particle(const Vec3d &position,
                          const Team7AirBoxParams &team7_air_box)
{
    return (position[0] - team7_air_box.air_lower[0] < boundary_shell_thickness) ||
           (team7_air_box.air_upper[0] - position[0] < boundary_shell_thickness) ||
           (position[1] - team7_air_box.air_lower[1] < boundary_shell_thickness) ||
           (team7_air_box.air_upper[1] - position[1] < boundary_shell_thickness) ||
           (position[2] - team7_air_box.air_lower[2] < boundary_shell_thickness) ||
           (team7_air_box.air_upper[2] - position[2] < boundary_shell_thickness);
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

struct GlobalBodyNeighborhoodSource
{
    size_t global_offset = 0;
    size_t particle_count = 0;
    Vecd *positions = nullptr;
    Real *volumetric_measure = nullptr;
    Real *smoothing_length_ratio = nullptr;
    const ParticleConfiguration *inner_configuration = nullptr;
    StdVec<std::pair<size_t, const ParticleConfiguration *>> contact_configurations;
};

struct GlobalDiscreteCloud
{
    StdVec<Vecd> positions_storage;
    StdVec<Real> volumetric_measure_storage;
    StdVec<Real> smoothing_length_ratio_storage;
    ParticleConfiguration particle_configuration;
    StdVec<StdVec<uint8_t>> neighbor_is_contact_storage;

    electromagnetics::LaplaceStructuredAPhiDiscreteView
    makeView(Real reference_smoothing_length)
    {
        return {positions_storage.size(),
                reference_smoothing_length,
                positions_storage.data(),
                volumetric_measure_storage.data(),
                smoothing_length_ratio_storage.data(),
                &particle_configuration,
                &neighbor_is_contact_storage};
    }
};

struct BodyOutputVariables
{
    size_t global_offset = 0;
    size_t particle_count = 0;
    Vecd *positions = nullptr;
    Real *volumetric_measure = nullptr;
    Real *smoothing_length_ratio = nullptr;

    Vecd *vtp_a_real = nullptr;
    Vecd *vtp_a_imag = nullptr;
    Real *vtp_phi_real = nullptr;
    Real *vtp_phi_imag = nullptr;
    Vecd *vtp_e_real = nullptr;
    Vecd *vtp_e_imag = nullptr;
    Vecd *vtp_j_real = nullptr;
    Vecd *vtp_j_imag = nullptr;
    Real *vtp_abs_a = nullptr;
    Real *vtp_abs_phi = nullptr;
    Real *vtp_abs_e = nullptr;
    Real *vtp_abs_j = nullptr;
    Real *vtp_joule = nullptr;
    Real *vtp_is_plate = nullptr;
    Real *vtp_is_coil = nullptr;
    Real *vtp_is_source = nullptr;
    Real *vtp_source_jz = nullptr;
};

BodyOutputVariables register_body_output_variables(RealBody &body,
                                                   size_t global_offset,
                                                   BodyStatesRecordingToVtp &write_states)
{
    BaseParticles &particles = body.getBaseParticles();
    BodyOutputVariables output_variables;
    output_variables.global_offset = global_offset;
    output_variables.particle_count = particles.TotalRealParticles();
    output_variables.positions = particles.getVariableDataByName<Vecd>("Position");
    output_variables.volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");
    output_variables.smoothing_length_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");

    output_variables.vtp_a_real =
        particles.registerStateVariableData<Vecd>("SolvedAReal", ZeroData<Vecd>::value);
    output_variables.vtp_a_imag =
        particles.registerStateVariableData<Vecd>("SolvedAImag", ZeroData<Vecd>::value);
    output_variables.vtp_phi_real =
        particles.registerStateVariableData<Real>("SolvedPhiReal", Real(0));
    output_variables.vtp_phi_imag =
        particles.registerStateVariableData<Real>("SolvedPhiImag", Real(0));
    output_variables.vtp_e_real =
        particles.registerStateVariableData<Vecd>("SolvedEReal", ZeroData<Vecd>::value);
    output_variables.vtp_e_imag =
        particles.registerStateVariableData<Vecd>("SolvedEImag", ZeroData<Vecd>::value);
    output_variables.vtp_j_real =
        particles.registerStateVariableData<Vecd>("SolvedJReal", ZeroData<Vecd>::value);
    output_variables.vtp_j_imag =
        particles.registerStateVariableData<Vecd>("SolvedJImag", ZeroData<Vecd>::value);
    output_variables.vtp_abs_a =
        particles.registerStateVariableData<Real>("SolvedAbsA", Real(0));
    output_variables.vtp_abs_phi =
        particles.registerStateVariableData<Real>("SolvedAbsPhi", Real(0));
    output_variables.vtp_abs_e =
        particles.registerStateVariableData<Real>("SolvedAbsE", Real(0));
    output_variables.vtp_abs_j =
        particles.registerStateVariableData<Real>("SolvedAbsJ", Real(0));
    output_variables.vtp_joule =
        particles.registerStateVariableData<Real>("SolvedJouleHeat", Real(0));
    output_variables.vtp_is_plate =
        particles.registerStateVariableData<Real>("RegionPlate", Real(0));
    output_variables.vtp_is_coil =
        particles.registerStateVariableData<Real>("RegionCoil", Real(0));
    output_variables.vtp_is_source =
        particles.registerStateVariableData<Real>("RegionSource", Real(0));
    output_variables.vtp_source_jz =
        particles.registerStateVariableData<Real>("SourceJz", Real(0));

    write_states.addToWrite<Vecd>(body, "SolvedAReal");
    write_states.addToWrite<Vecd>(body, "SolvedAImag");
    write_states.addToWrite<Real>(body, "SolvedPhiReal");
    write_states.addToWrite<Real>(body, "SolvedPhiImag");
    write_states.addToWrite<Vecd>(body, "SolvedEReal");
    write_states.addToWrite<Vecd>(body, "SolvedEImag");
    write_states.addToWrite<Vecd>(body, "SolvedJReal");
    write_states.addToWrite<Vecd>(body, "SolvedJImag");
    write_states.addToWrite<Real>(body, "SolvedAbsA");
    write_states.addToWrite<Real>(body, "SolvedAbsPhi");
    write_states.addToWrite<Real>(body, "SolvedAbsE");
    write_states.addToWrite<Real>(body, "SolvedAbsJ");
    write_states.addToWrite<Real>(body, "SolvedJouleHeat");
    write_states.addToWrite<Real>(body, "RegionPlate");
    write_states.addToWrite<Real>(body, "RegionCoil");
    write_states.addToWrite<Real>(body, "RegionSource");
    write_states.addToWrite<Real>(body, "SourceJz");
    write_states.addToWrite<Real>(body, "SmoothingLengthRatio");

    return output_variables;
}

void scatter_global_fields_to_body(const BodyOutputVariables &body_output,
                                   const electromagnetics::LaplaceStructuredAPhiFields &solved_fields,
                                   const StdVec<Real> &abs_a,
                                   const StdVec<Real> &abs_phi,
                                   const StdVec<Real> &abs_e,
                                   const StdVec<Real> &abs_j,
                                   const StdVec<Real> &joule,
                                   const StdVec<bool> &is_plate,
                                   const StdVec<bool> &is_coil,
                                   const StdVec<bool> &is_source,
                                   const StdVec<Vecd> &source_current_density)
{
    for (size_t local_i = 0; local_i != body_output.particle_count; ++local_i)
    {
        size_t global_i = body_output.global_offset + local_i;
        body_output.vtp_a_real[local_i] = Vecd(
            solved_fields.ax[global_i].real(),
            solved_fields.ay[global_i].real(),
            solved_fields.az[global_i].real());
        body_output.vtp_a_imag[local_i] = Vecd(
            solved_fields.ax[global_i].imag(),
            solved_fields.ay[global_i].imag(),
            solved_fields.az[global_i].imag());
        body_output.vtp_phi_real[local_i] = solved_fields.phi[global_i].real();
        body_output.vtp_phi_imag[local_i] = solved_fields.phi[global_i].imag();
        body_output.vtp_abs_a[local_i] = abs_a[global_i];
        body_output.vtp_abs_phi[local_i] = abs_phi[global_i];
        body_output.vtp_abs_e[local_i] = abs_e[global_i];
        body_output.vtp_abs_j[local_i] = abs_j[global_i];
        body_output.vtp_joule[local_i] = joule[global_i];
        body_output.vtp_is_plate[local_i] = is_plate[global_i] ? 1.0 : 0.0;
        body_output.vtp_is_coil[local_i] = is_coil[global_i] ? 1.0 : 0.0;
        body_output.vtp_is_source[local_i] = is_source[global_i] ? 1.0 : 0.0;
        body_output.vtp_source_jz[local_i] = source_current_density[global_i][2];

        if (!solved_fields.electric_field.empty())
        {
            const auto &electric_field = solved_fields.electric_field[global_i];
            body_output.vtp_e_real[local_i] = Vecd(
                electric_field[0].real(),
                electric_field[1].real(),
                electric_field[2].real());
            body_output.vtp_e_imag[local_i] = Vecd(
                electric_field[0].imag(),
                electric_field[1].imag(),
                electric_field[2].imag());
        }
        if (!solved_fields.current_density.empty())
        {
            const auto &current_density = solved_fields.current_density[global_i];
            body_output.vtp_j_real[local_i] = Vecd(
                current_density[0].real(),
                current_density[1].real(),
                current_density[2].real());
            body_output.vtp_j_imag[local_i] = Vecd(
                current_density[0].imag(),
                current_density[1].imag(),
                current_density[2].imag());
        }
    }
}

void append_neighborhood_entries(Neighborhood &target_neighborhood,
                                 const Neighborhood &source_neighborhood,
                                 StdVec<uint8_t> &target_neighbor_is_contact,
                                 size_t target_index_offset,
                                 bool is_contact_source)
{
    for (size_t n = 0; n != source_neighborhood.current_size_; ++n)
    {
        target_neighborhood.j_.push_back(source_neighborhood.j_[n] + target_index_offset);
        target_neighborhood.W_ij_.push_back(source_neighborhood.W_ij_[n]);
        target_neighborhood.dW_ij_.push_back(source_neighborhood.dW_ij_[n]);
        target_neighborhood.r_ij_.push_back(source_neighborhood.r_ij_[n]);
        target_neighborhood.e_ij_.push_back(source_neighborhood.e_ij_[n]);
        target_neighbor_is_contact.push_back(is_contact_source ? 1u : 0u);
        target_neighborhood.current_size_++;
        target_neighborhood.allocated_size_++;
    }
}

GlobalDiscreteCloud build_global_discrete_cloud(const StdVec<GlobalBodyNeighborhoodSource> &body_sources,
                                                Real reference_smoothing_length)
{
    GlobalDiscreteCloud global_cloud;
    size_t total_particles = 0;
    for (const auto &source : body_sources)
    {
        total_particles += source.particle_count;
    }

    global_cloud.positions_storage.resize(total_particles, ZeroData<Vecd>::value);
    global_cloud.volumetric_measure_storage.resize(total_particles, Real(0));
    global_cloud.smoothing_length_ratio_storage.resize(total_particles, Real(1));
    global_cloud.particle_configuration.resize(total_particles);
    global_cloud.neighbor_is_contact_storage.resize(total_particles);

    for (const auto &source : body_sources)
    {
        for (size_t local_i = 0; local_i != source.particle_count; ++local_i)
        {
            size_t global_i = source.global_offset + local_i;
            global_cloud.positions_storage[global_i] = source.positions[local_i];
            global_cloud.volumetric_measure_storage[global_i] = source.volumetric_measure[local_i];
            global_cloud.smoothing_length_ratio_storage[global_i] =
                source.smoothing_length_ratio != nullptr ? source.smoothing_length_ratio[local_i] : Real(1);

            if (source.inner_configuration)
            {
                append_neighborhood_entries(global_cloud.particle_configuration[global_i],
                                            (*source.inner_configuration)[local_i],
                                            global_cloud.neighbor_is_contact_storage[global_i],
                                            source.global_offset, false);
            }
            for (const auto &contact_source : source.contact_configurations)
            {
                append_neighborhood_entries(global_cloud.particle_configuration[global_i],
                                            (*contact_source.second)[local_i],
                                            global_cloud.neighbor_is_contact_storage[global_i],
                                            contact_source.first, true);
            }
        }
    }

    for (size_t global_i = 0; global_i != total_particles; ++global_i)
    {
        Neighborhood &neighborhood = global_cloud.particle_configuration[global_i];
        if (neighborhood.current_size_ == 0)
        {
            continue;
        }
        StdVec<size_t> permutation(neighborhood.current_size_);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(),
                  [&](size_t lhs, size_t rhs)
                  { return neighborhood.j_[lhs] < neighborhood.j_[rhs]; });

        Neighborhood sorted_neighborhood;
        StdVec<uint8_t> sorted_neighbor_is_contact;
        for (size_t perm_cursor = 0; perm_cursor != permutation.size();)
        {
            size_t base_idx = permutation[perm_cursor];
            size_t neighbor_index = neighborhood.j_[base_idx];
            Real merged_weight = neighborhood.W_ij_[base_idx];
            Real merged_dweight = neighborhood.dW_ij_[base_idx];
            Real merged_distance = neighborhood.r_ij_[base_idx];
            Vecd merged_direction = neighborhood.e_ij_[base_idx];
            uint8_t merged_is_contact = global_cloud.neighbor_is_contact_storage[global_i][base_idx];

            size_t next_cursor = perm_cursor + 1;
            while (next_cursor != permutation.size() &&
                   neighborhood.j_[permutation[next_cursor]] == neighbor_index)
            {
                size_t duplicate_idx = permutation[next_cursor];
                merged_weight += neighborhood.W_ij_[duplicate_idx];
                merged_dweight += neighborhood.dW_ij_[duplicate_idx];
                merged_distance = SMIN(merged_distance, neighborhood.r_ij_[duplicate_idx]);
                merged_direction += neighborhood.e_ij_[duplicate_idx];
                merged_is_contact = SMAX(merged_is_contact, global_cloud.neighbor_is_contact_storage[global_i][duplicate_idx]);
                next_cursor++;
            }

            if (merged_direction.squaredNorm() > TinyReal)
            {
                merged_direction /= merged_direction.norm();
            }
            else
            {
                merged_direction = neighborhood.e_ij_[base_idx];
            }

            sorted_neighborhood.j_.push_back(neighbor_index);
            sorted_neighborhood.W_ij_.push_back(merged_weight);
            sorted_neighborhood.dW_ij_.push_back(merged_dweight);
            sorted_neighborhood.r_ij_.push_back(merged_distance);
            sorted_neighborhood.e_ij_.push_back(merged_direction);
            sorted_neighbor_is_contact.push_back(merged_is_contact);
            sorted_neighborhood.current_size_++;
            sorted_neighborhood.allocated_size_++;
            perm_cursor = next_cursor;
        }
        neighborhood = std::move(sorted_neighborhood);
        global_cloud.neighbor_is_contact_storage[global_i] = std::move(sorted_neighbor_is_contact);
    }

    (void)reference_smoothing_length;
    return global_cloud;
}

Vecd source_profile_current_density(const Vec3d &position,
                                    const Vec3d &coil_center,
                                    const Vec3d &coil_axis,
                                    const Vec3d &fallback_direction,
                                    Real amplitude,
                                    bool use_circular_source)
{
    if (use_circular_source)
    {
        Vecd tangent = TangentialDirectionAroundAxis(
            position,
            coil_center,
            NormalizeVectorOrDefault(coil_axis, Vecd(0.0, 0.0, 1.0)));
        if (tangent.squaredNorm() > TinyReal)
        {
            return amplitude * tangent;
        }
    }
    return amplitude * NormalizeVectorOrDefault(fallback_direction, Vecd(0.0, 0.0, 1.0));
}
} // namespace

int main(int ac, char *av[])
{
    using Clock = std::chrono::steady_clock;
    auto seconds_since = [](const Clock::time_point &begin,
                            const Clock::time_point &end) -> double
    {
        return std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count();
    };

    const auto time_begin = Clock::now();
    TriangleMeshShapeSTL coil_mesh(path_coil_stl, Vec3d::Zero(), 1.0, "Team7CoilMesh");
    TriangleMeshShapeSTL plate_mesh(path_plate_stl, Vec3d::Zero(), 1.0, "Team7PlateMesh");
    BoundingBoxd coil_bounds = coil_mesh.findBounds();
    BoundingBoxd plate_bounds = plate_mesh.findBounds();
    Vec3d coil_center = 0.5 * (coil_bounds.lower_ + coil_bounds.upper_);
    Team7AirBoxParams team7_air_box =
        resolve_team7_air_box_from_geometry(requested_team7_air_box, coil_bounds, plate_bounds);
    const auto time_geometry_ready = Clock::now();

    BoundingBoxd system_domain_bounds(
        team7_air_box.air_lower - Vec3d::Ones() * boundary_width,
        team7_air_box.air_upper + Vec3d::Ones() * boundary_width);

    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    std::string output_tag = get_env_string_local("EM_APHI_LAPLACE_OUTPUT_TAG", "");
    if (!output_tag.empty())
    {
        io_environment.appendOutputFolder(output_tag);
    }

    std::cout << std::scientific << std::setprecision(6)
              << "[em-aphi-laplace-team7-stl-em-only-config] dp=" << dp_0
              << ", air_lower=(" << team7_air_box.air_lower[0] << "," << team7_air_box.air_lower[1] << "," << team7_air_box.air_lower[2] << ")"
              << ", air_upper=(" << team7_air_box.air_upper[0] << "," << team7_air_box.air_upper[1] << "," << team7_air_box.air_upper[2] << ")"
              << ", used_small_air_box=" << (team7_air_box.used_small_air_box ? 1 : 0)
              << ", used_tight_stl_air_box=" << (team7_air_box.used_tight_stl_air_box ? 1 : 0)
              << ", tight_stl_margin=" << team7_air_box.tight_stl_margin
              << ", sigma_air=" << sigma_air
              << ", sigma_plate=" << sigma_plate
              << ", sigma_coil=" << sigma_coil
              << ", nu_air=" << nu_air
              << ", nu_plate=" << nu_plate
              << ", nu_coil=" << nu_coil
              << ", frequency_hz=" << frequency_hz
              << ", coil_current_rms=" << coil_current_rms
              << ", coil_turns=" << coil_turns
              << ", coil_current_peak_factor=" << coil_current_peak_factor
              << ", coil_effective_area_input_geom=" << coil_effective_area_input_geom
              << ", geom_length_to_m=" << geom_length_to_m
              << ", use_circular_coil_source=" << (use_circular_coil_source ? 1 : 0)
              << ", use_three_body_layout=" << (use_three_body_layout ? 1 : 0)
              << ", use_multi_resolution=" << (use_multi_resolution ? 1 : 0)
              << ", dp_air_coarse=" << (use_multi_resolution ? dp_air_coarse_input : dp_0)
              << ", air_refinement_levels=" << (use_multi_resolution ? air_refinement_levels_input : 0)
              << ", use_magnetic_only_air=" << (use_magnetic_only_air ? 1 : 0)
              << ", use_magnetic_only_coil=" << (use_magnetic_only_coil ? 1 : 0)
              << ", sigma_solve_regularization=" << sigma_solve_regularization
              << ", contact_gradient_scale=" << contact_gradient_scale
              << ", contact_diffusion_scale=" << contact_diffusion_scale
              << ", coil_source_axis=(" << coil_source_axis[0] << "," << coil_source_axis[1] << "," << coil_source_axis[2] << ")"
              << ", enable_block_scaling=" << enable_block_scaling
              << ", phi_block_scale=" << phi_block_scale
              << ", enable_diagonal_equilibration=" << enable_diagonal_equilibration
              << ", diagonal_equilibration_iterations=" << diagonal_equilibration_iterations
              << ", solver_backend=" << solver_backend
              << ", iterative_max_iterations=" << iterative_max_iterations
              << ", iterative_tolerance=" << iterative_tolerance
              << std::endl;

    if (use_three_body_layout)
    {
        SolidBody coil_body(sph_system, makeShared<CoilShape>("Team7STLEMOnlyCoil"));
        coil_body.defineAdaptation<SPHAdaptation>(1.15, dp_0 / dp_coil_input);
        coil_body.defineMaterial<Solid>();
        coil_body.defineBodyLevelSetShape();
        coil_body.generateParticles<BaseParticles, Lattice>();
        coil_body.getBaseParticles().registerStateVariableData<Real>("SmoothingLengthRatio", Real(1.0));

        SolidBody plate_body(sph_system, makeShared<PlateShape>("Team7STLEMOnlyPlate"));
        plate_body.defineAdaptation<SPHAdaptation>(1.15, dp_0 / dp_plate_input);
        plate_body.defineMaterial<Solid>();
        plate_body.defineBodyLevelSetShape();
        plate_body.generateParticles<BaseParticles, Lattice>();
        plate_body.getBaseParticles().registerStateVariableData<Real>("SmoothingLengthRatio", Real(1.0));

        auto &inner_boundary_shape = sph_system.addShape<InnerBoundaryShape>("Team7InnerBoundary");
        Real dp_air_coarse = use_multi_resolution ? dp_air_coarse_input : dp_0;
        int air_refinement_levels = use_multi_resolution ? air_refinement_levels_input : 0;
        AdaptiveNearInnerSurface air_adaptation(
            dp_air_coarse, 1.15, 1.0, air_refinement_levels, &inner_boundary_shape);
        auto &air_body = sph_system.addAdaptiveBody<FluidBody, AdaptiveNearInnerSurface>(
            air_adaptation, makeShared<AirShape>("Team7STLEMOnlyAir", team7_air_box.air_lower, team7_air_box.air_upper));
        air_body.defineMaterial<Solid>();
        air_body.generateParticles<BaseParticles, Lattice>();

        BodyStatesRecordingToVtp write_states(sph_system);

        InnerRelation coil_inner(coil_body);
        InnerRelation plate_inner(plate_body);
        AdaptiveInnerRelation air_inner(air_body);
        AdaptiveContactRelation coil_contact(coil_body, {&plate_body, &air_body});
        AdaptiveContactRelation plate_contact(plate_body, {&coil_body, &air_body});
        AdaptiveContactRelation air_contact(air_body, {&coil_body, &plate_body});

        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        const auto time_particles_ready = Clock::now();

        const size_t coil_particles = coil_body.getBaseParticles().TotalRealParticles();
        const size_t plate_particles = plate_body.getBaseParticles().TotalRealParticles();
        const size_t air_particles = air_body.getBaseParticles().TotalRealParticles();
        const size_t total_particles = coil_particles + plate_particles + air_particles;
        const size_t coil_offset = 0;
        const size_t plate_offset = coil_offset + coil_particles;
        const size_t air_offset = plate_offset + plate_particles;

        BodyOutputVariables coil_output =
            register_body_output_variables(coil_body, coil_offset, write_states);
        BodyOutputVariables plate_output =
            register_body_output_variables(plate_body, plate_offset, write_states);
        BodyOutputVariables air_output =
            register_body_output_variables(air_body, air_offset, write_states);

        StdVec<GlobalBodyNeighborhoodSource> body_sources(3);
        body_sources[0] = {coil_offset,
                           coil_particles,
                           coil_output.positions,
                           coil_output.volumetric_measure,
                           coil_output.smoothing_length_ratio,
                           &coil_inner.inner_configuration_,
                           {{plate_offset, &coil_contact.contact_configuration_[0]},
                            {air_offset, &coil_contact.contact_configuration_[1]}}};
        body_sources[1] = {plate_offset,
                           plate_particles,
                           plate_output.positions,
                           plate_output.volumetric_measure,
                           plate_output.smoothing_length_ratio,
                           &plate_inner.inner_configuration_,
                           {{coil_offset, &plate_contact.contact_configuration_[0]},
                            {air_offset, &plate_contact.contact_configuration_[1]}}};
        body_sources[2] = {air_offset,
                           air_particles,
                           air_output.positions,
                           air_output.volumetric_measure,
                           air_output.smoothing_length_ratio,
                           &air_inner.inner_configuration_,
                           {{coil_offset, &air_contact.contact_configuration_[0]},
                            {plate_offset, &air_contact.contact_configuration_[1]}}};

        GlobalDiscreteCloud global_cloud =
            build_global_discrete_cloud(body_sources, dp_0);
        electromagnetics::LaplaceStructuredAPhiDiscreteView discrete_view =
            global_cloud.makeView(dp_0);

        Vecd *positions = global_cloud.positions_storage.data();
        Real *volumetric_measure = global_cloud.volumetric_measure_storage.data();
        StdVec<Real> global_h_ratio(total_particles, 1.0);
        for (size_t i = 0; i != coil_particles; ++i)
        {
            global_h_ratio[coil_offset + i] = coil_output.smoothing_length_ratio[i];
        }
        for (size_t i = 0; i != plate_particles; ++i)
        {
            global_h_ratio[plate_offset + i] = plate_output.smoothing_length_ratio[i];
        }
        for (size_t i = 0; i != air_particles; ++i)
        {
            global_h_ratio[air_offset + i] = air_output.smoothing_length_ratio[i];
        }

        StdVec<Real> sigma_physical(total_particles, sigma_air);
        StdVec<Real> sigma_em(total_particles, use_magnetic_only_air ? 0.0 : sigma_air);
        StdVec<Real> sigma_solve(total_particles, use_magnetic_only_air ? sigma_solve_regularization : sigma_air);
        StdVec<Real> nu(total_particles, nu_air);
        StdVec<bool> is_plate(total_particles, false);
        StdVec<bool> is_coil(total_particles, false);
        StdVec<bool> is_source(total_particles, false);
        StdVec<Vecd> source_current_density(total_particles, ZeroData<Vecd>::value);

        Real coil_volume_total_geom = 0.0;
        Real coil_radius_volume_moment = 0.0;
        for (size_t i = 0; i != coil_particles; ++i)
        {
            size_t global_i = coil_offset + i;
            is_coil[global_i] = true;
            sigma_physical[global_i] = sigma_coil;
            sigma_em[global_i] = use_magnetic_only_coil ? 0.0 : sigma_coil;
            sigma_solve[global_i] = use_magnetic_only_coil ? sigma_solve_regularization : sigma_coil;
            nu[global_i] = nu_coil;
            coil_volume_total_geom += volumetric_measure[global_i];
            if (use_circular_coil_source)
            {
                Vecd axis_used = NormalizeVectorOrDefault(coil_source_axis, Vecd(0.0, 0.0, 1.0));
                Vecd relative = positions[global_i] - coil_center;
                Vecd radial = relative - relative.dot(axis_used) * axis_used;
                coil_radius_volume_moment += volumetric_measure[global_i] * radial.norm();
            }
        }
        for (size_t i = 0; i != plate_particles; ++i)
        {
            size_t global_i = plate_offset + i;
            is_plate[global_i] = true;
            sigma_physical[global_i] = sigma_plate;
            sigma_em[global_i] = sigma_plate;
            sigma_solve[global_i] = sigma_plate;
            nu[global_i] = nu_plate;
        }

        Real safe_geom_length_to_m = SMAX(geom_length_to_m, TinyReal);
        Real effective_area_geom = coil_effective_area_input_geom;
        if (!(effective_area_geom > TinyReal))
        {
            if (use_circular_coil_source)
            {
                Real mean_radius_geom =
                    coil_radius_volume_moment / (coil_volume_total_geom + TinyReal);
                Real projection_length_geom = 2.0 * Pi * mean_radius_geom;
                effective_area_geom =
                    coil_volume_total_geom / (projection_length_geom + TinyReal);
            }
            else
            {
                Vecd dir_used = NormalizeVectorOrDefault(fallback_source_direction, Vecd(0.0, 0.0, 1.0));
                Vecd coil_extent = coil_bounds.upper_ - coil_bounds.lower_;
                Real projection_length_geom =
                    fabs(dir_used[0]) * coil_extent[0] +
                    fabs(dir_used[1]) * coil_extent[1] +
                    fabs(dir_used[2]) * coil_extent[2];
                effective_area_geom =
                    coil_volume_total_geom / (projection_length_geom + TinyReal);
            }
        }
        Real effective_area_si = effective_area_geom * safe_geom_length_to_m * safe_geom_length_to_m;
        Real current_peak = SMAX(coil_current_peak_factor, static_cast<Real>(1.0)) *
                            SMAX(coil_current_rms, static_cast<Real>(0.0));
        Real source_current_density_magnitude =
            coil_turns * current_peak / (effective_area_si + TinyReal);

        std::cout << std::scientific << std::setprecision(6)
                  << "[em-aphi-laplace-team7-stl-em-only-source] coil_volume_total_geom=" << coil_volume_total_geom
                  << ", effective_area_geom=" << effective_area_geom
                  << ", effective_area_si=" << effective_area_si
                  << ", source_current_density_magnitude=" << source_current_density_magnitude
                  << std::endl;

        size_t source_particles = 0;
        for (size_t i = 0; i != coil_particles; ++i)
        {
            size_t global_i = coil_offset + i;
            Vecd js = source_profile_current_density(
                positions[global_i], coil_center, coil_source_axis,
                fallback_source_direction, source_current_density_magnitude,
                use_circular_coil_source);
            source_current_density[global_i] = js;
            is_source[global_i] = js.squaredNorm() > TinyReal;
            if (is_source[global_i])
            {
                source_particles++;
            }
        }

        electromagnetics::LaplaceStructuredAPhiBoundaryCondition boundary_condition;
        boundary_condition.is_dirichlet.resize(total_particles, false);
        boundary_condition.ax.resize(total_particles, Complex(0.0, 0.0));
        boundary_condition.ay.resize(total_particles, Complex(0.0, 0.0));
        boundary_condition.az.resize(total_particles, Complex(0.0, 0.0));
        boundary_condition.phi.resize(total_particles, Complex(0.0, 0.0));

        size_t reference_phi_index = plate_offset;
        for (size_t i = 0; i != total_particles; ++i)
        {
            boundary_condition.is_dirichlet[i] =
                is_boundary_particle(positions[i], team7_air_box);
            if (is_plate[i] && !boundary_condition.is_dirichlet[i])
            {
                reference_phi_index = i;
                break;
            }
        }

        electromagnetics::LaplaceStructuredAPhiParameters parameters;
        parameters.angular_frequency = 2.0 * Pi * frequency_hz;
        parameters.geom_length_to_m = geom_length_to_m;
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
        parameters.contact_gradient_scale = contact_gradient_scale;
        parameters.contact_diffusion_scale = contact_diffusion_scale;

        electromagnetics::LaplaceStructuredAPhiEigenSolver solver(discrete_view, parameters);
        electromagnetics::LaplaceStructuredAPhiEigenSolver *active_solver = &solver;
        std::unique_ptr<electromagnetics::LaplaceStructuredAPhiEigenSolver> fallback_solver;
        electromagnetics::VecC rhs(static_cast<int>(4 * total_particles));
        rhs.setZero();
        for (size_t i = 0; i != total_particles; ++i)
        {
            rhs[solver.idxAx(i)] = Complex(source_current_density[i][0], 0.0);
            rhs[solver.idxAy(i)] = Complex(source_current_density[i][1], 0.0);
            rhs[solver.idxAz(i)] = Complex(source_current_density[i][2], 0.0);
        }
        const auto time_setup_ready = Clock::now();

        solver.assembleSystem(sigma_solve, nu, rhs, boundary_condition);
        const auto time_assembled = Clock::now();

        electromagnetics::LaplaceStructuredAPhiFields solved_fields;
        std::string solver_message;
        bool solve_success = solver.solve(solved_fields, solver_message);
        if (!solve_success)
        {
            std::cerr << "[em-aphi-laplace-team7-stl-em-only] primary solve failed"
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
                    discrete_view, fallback_parameters);
                fallback_solver->assembleSystem(sigma_solve, nu, rhs, boundary_condition);

                solve_success = fallback_solver->solve(solved_fields, solver_message);
                if (!solve_success)
                {
                    std::cerr << "[em-aphi-laplace-team7-stl-em-only] fallback solve failed"
                              << ", backend=" << fallback_solver->LastSolverBackend()
                              << ", iterations=" << fallback_solver->LastSolverIterations()
                              << ", estimated_error=" << fallback_solver->LastSolverEstimatedError()
                              << ", message=" << solver_message
                              << std::endl;
                    return 1;
                }

                std::cout << "[em-aphi-laplace-team7-stl-em-only] fallback solve succeeded"
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
        const auto time_solved = Clock::now();

        electromagnetics::LaplaceStructuredAPhiDiagnostics solved_diagnostics =
            active_solver->computeDiagnostics(sigma_em, nu, solved_fields);

        PhysicalRegionSummary plate_summary;
        PhysicalRegionSummary air_summary;
        PhysicalRegionSummary coil_summary;
        PhysicalRegionSummary source_summary;
        StdVec<Real> abs_a(total_particles, 0.0);
        StdVec<Real> abs_phi(total_particles, 0.0);
        StdVec<Real> abs_e(total_particles, 0.0);
        StdVec<Real> abs_j(total_particles, 0.0);
        StdVec<Real> joule(total_particles, 0.0);

        for (size_t i = 0; i != total_particles; ++i)
        {
            abs_a[i] = std::sqrt(std::norm(solved_fields.ax[i]) +
                                 std::norm(solved_fields.ay[i]) +
                                 std::norm(solved_fields.az[i]));
            abs_phi[i] = std::abs(solved_fields.phi[i]);
            abs_e[i] = solved_fields.electric_field.empty() ? 0.0 :
                static_cast<Real>(solved_fields.electric_field[i].norm());
            abs_j[i] = solved_fields.current_density.empty() ? 0.0 :
                static_cast<Real>(solved_fields.current_density[i].norm());
            joule[i] = solved_fields.joule_heat.empty() ? 0.0 : solved_fields.joule_heat[i];

            auto accumulate_summary = [&](PhysicalRegionSummary &summary)
            {
                summary.particles++;
                summary.mean_abs_a += abs_a[i];
                summary.max_abs_a = SMAX(summary.max_abs_a, abs_a[i]);
                summary.mean_abs_phi += abs_phi[i];
                summary.max_abs_phi = SMAX(summary.max_abs_phi, abs_phi[i]);
                summary.mean_abs_e += abs_e[i];
                summary.max_abs_e = SMAX(summary.max_abs_e, abs_e[i]);
                summary.mean_abs_j += abs_j[i];
                summary.max_abs_j = SMAX(summary.max_abs_j, abs_j[i]);
                summary.mean_joule += joule[i];
                summary.max_joule = SMAX(summary.max_joule, joule[i]);
            };

            if (is_plate[i])
            {
                accumulate_summary(plate_summary);
            }
            else if (is_coil[i])
            {
                accumulate_summary(coil_summary);
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

        finalize_region_summary(plate_summary);
        finalize_region_summary(air_summary);
        finalize_region_summary(coil_summary);
        finalize_region_summary(source_summary);

        scatter_global_fields_to_body(coil_output, solved_fields, abs_a, abs_phi, abs_e, abs_j, joule,
                                      is_plate, is_coil, is_source, source_current_density);
        scatter_global_fields_to_body(plate_output, solved_fields, abs_a, abs_phi, abs_e, abs_j, joule,
                                      is_plate, is_coil, is_source, source_current_density);
        scatter_global_fields_to_body(air_output, solved_fields, abs_a, abs_phi, abs_e, abs_j, joule,
                                      is_plate, is_coil, is_source, source_current_density);
        const auto time_postprocess_ready = Clock::now();

        if (write_vtp)
        {
            coil_body.setNewlyUpdated();
            plate_body.setNewlyUpdated();
            air_body.setNewlyUpdated();
            write_states.writeToFile(0.0);
        }
        const auto time_output_ready = Clock::now();

        const double time_geometry_s = seconds_since(time_begin, time_geometry_ready);
        const double time_particle_generation_s = seconds_since(time_geometry_ready, time_particles_ready);
        const double time_region_setup_s = seconds_since(time_particles_ready, time_setup_ready);
        const double time_assemble_s = seconds_since(time_setup_ready, time_assembled);
        const double time_solve_s = seconds_since(time_assembled, time_solved);
        const double time_postprocess_s = seconds_since(time_solved, time_postprocess_ready);
        const double time_output_s = seconds_since(time_postprocess_ready, time_output_ready);
        const double time_total_s = seconds_since(time_begin, time_output_ready);

        std::ofstream summary_file(
            io_environment.OutputFolder() + "/em_aphi_laplace_team7_stl_em_only_summary.csv",
            std::ios::out | std::ios::trunc);
        summary_file << std::setprecision(12);
        summary_file << "solve_success,particles,plate_particles,coil_particles,air_particles,source_particles,"
                     << "dp,air_box_lower_x,air_box_lower_y,air_box_lower_z,air_box_upper_x,air_box_upper_y,air_box_upper_z,"
                     << "used_small_air_box,used_tight_stl_air_box,tight_stl_margin,sigma_air,sigma_plate,sigma_coil,sigma_solve_regularization,nu_air,nu_plate,nu_coil,frequency_hz,coil_current_rms,coil_turns,coil_current_peak_factor,coil_effective_area_input_geom,coil_effective_area_used_geom,coil_effective_area_used_m2,geom_length_to_m,use_circular_coil_source,use_three_body_layout,use_multi_resolution,dp_air_coarse,air_refinement_levels,use_magnetic_only_air,use_magnetic_only_coil,coil_source_axis_x,coil_source_axis_y,coil_source_axis_z,"
                     << "solver_backend,solver_iterations,solver_estimated_error,"
                     << "used_phi_block_scale,used_diagonal_equilibration_iterations,"
                     << "solved_div_a_l2,solved_current_continuity_l2,solved_phi_equation_residual_l2,"
                     << "solved_physical_current_continuity_l2,solved_relative_current_continuity_l2,"
                     << "solved_linear_residual_norm,solved_total_joule_power,solved_gauge_penalty_energy,solved_magnetic_energy,"
                     << "time_geometry_s,time_particle_generation_s,time_region_setup_s,time_assemble_s,time_solve_s,time_postprocess_s,time_output_s,time_total_s,"
                     << "solved_max_abs_a,solved_max_abs_e,solved_max_abs_j,"
                     << "plate_mean_abs_a,plate_max_abs_a,plate_mean_abs_phi,plate_max_abs_phi,"
                     << "plate_mean_abs_e,plate_max_abs_e,plate_mean_abs_j,plate_max_abs_j,plate_mean_joule,plate_max_joule,"
                     << "coil_mean_abs_a,coil_max_abs_a,coil_mean_abs_phi,coil_max_abs_phi,"
                     << "coil_mean_abs_e,coil_max_abs_e,coil_mean_abs_j,coil_max_abs_j,coil_mean_joule,coil_max_joule,"
                     << "air_mean_abs_a,air_max_abs_a,air_mean_abs_phi,air_max_abs_phi,"
                     << "air_mean_abs_e,air_max_abs_e,air_mean_abs_j,air_max_abs_j,air_mean_joule,air_max_joule,"
                     << "source_mean_abs_a,source_max_abs_a,source_mean_abs_phi,source_max_abs_phi,"
                     << "source_mean_abs_e,source_max_abs_e,source_mean_abs_j,source_max_abs_j,source_mean_joule,source_max_joule\n";

        summary_file << (solve_success ? 1 : 0) << ","
                     << total_particles << ","
                     << plate_particles << ","
                     << coil_particles << ","
                     << air_particles << ","
                     << source_particles << ","
                     << dp_0 << ","
                     << team7_air_box.air_lower[0] << ","
                     << team7_air_box.air_lower[1] << ","
                     << team7_air_box.air_lower[2] << ","
                     << team7_air_box.air_upper[0] << ","
                     << team7_air_box.air_upper[1] << ","
                     << team7_air_box.air_upper[2] << ","
                     << (team7_air_box.used_small_air_box ? 1 : 0) << ","
                     << (team7_air_box.used_tight_stl_air_box ? 1 : 0) << ","
                     << team7_air_box.tight_stl_margin << ","
                     << sigma_air << ","
                     << sigma_plate << ","
                     << sigma_coil << ","
                     << sigma_solve_regularization << ","
                     << nu_air << ","
                     << nu_plate << ","
                     << nu_coil << ","
                     << frequency_hz << ","
                     << coil_current_rms << ","
                     << coil_turns << ","
                     << coil_current_peak_factor << ","
                     << coil_effective_area_input_geom << ","
                     << effective_area_geom << ","
                     << effective_area_si << ","
                     << geom_length_to_m << ","
                     << (use_circular_coil_source ? 1 : 0) << ","
                     << (use_three_body_layout ? 1 : 0) << ","
                     << (use_multi_resolution ? 1 : 0) << ","
                     << dp_air_coarse << ","
                     << air_refinement_levels << ","
                     << (use_magnetic_only_air ? 1 : 0) << ","
                     << (use_magnetic_only_coil ? 1 : 0) << ","
                     << coil_source_axis[0] << ","
                     << coil_source_axis[1] << ","
                     << coil_source_axis[2] << ","
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
                     << time_geometry_s << ","
                     << time_particle_generation_s << ","
                     << time_region_setup_s << ","
                     << time_assemble_s << ","
                     << time_solve_s << ","
                     << time_postprocess_s << ","
                     << time_output_s << ","
                     << time_total_s << ","
                     << solved_diagnostics.max_abs_a << ","
                     << solved_diagnostics.max_abs_e << ","
                     << solved_diagnostics.max_abs_j << ","
                     << plate_summary.mean_abs_a << ","
                     << plate_summary.max_abs_a << ","
                     << plate_summary.mean_abs_phi << ","
                     << plate_summary.max_abs_phi << ","
                     << plate_summary.mean_abs_e << ","
                     << plate_summary.max_abs_e << ","
                     << plate_summary.mean_abs_j << ","
                     << plate_summary.max_abs_j << ","
                     << plate_summary.mean_joule << ","
                     << plate_summary.max_joule << ","
                     << coil_summary.mean_abs_a << ","
                     << coil_summary.max_abs_a << ","
                     << coil_summary.mean_abs_phi << ","
                     << coil_summary.max_abs_phi << ","
                     << coil_summary.mean_abs_e << ","
                     << coil_summary.max_abs_e << ","
                     << coil_summary.mean_abs_j << ","
                     << coil_summary.max_abs_j << ","
                     << coil_summary.mean_joule << ","
                     << coil_summary.max_joule << ","
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
                io_environment.OutputFolder() + "/em_aphi_laplace_team7_stl_em_only_particles.csv",
                std::ios::out | std::ios::trunc);
            particle_file << std::setprecision(12);
            particle_file << "particle_id,x,y,z,is_boundary,is_plate,is_coil,is_source,smoothing_length_ratio,sigma_physical,sigma_em,sigma_solve,nu,source_jx,source_jy,source_jz,"
                          << "ax_real,ax_imag,ay_real,ay_imag,az_real,az_imag,phi_real,phi_imag,"
                          << "abs_a,abs_phi,abs_e,abs_j,joule_heat\n";
            for (size_t i = 0; i != total_particles; ++i)
            {
                particle_file << i << ","
                              << positions[i][0] << ","
                              << positions[i][1] << ","
                              << positions[i][2] << ","
                              << (boundary_condition.is_dirichlet[i] ? 1 : 0) << ","
                              << (is_plate[i] ? 1 : 0) << ","
                              << (is_coil[i] ? 1 : 0) << ","
                              << (is_source[i] ? 1 : 0) << ","
                              << global_h_ratio[i] << ","
                              << sigma_physical[i] << ","
                              << sigma_em[i] << ","
                              << sigma_solve[i] << ","
                              << nu[i] << ","
                              << source_current_density[i][0] << ","
                              << source_current_density[i][1] << ","
                              << source_current_density[i][2] << ","
                              << solved_fields.ax[i].real() << ","
                              << solved_fields.ax[i].imag() << ","
                              << solved_fields.ay[i].real() << ","
                              << solved_fields.ay[i].imag() << ","
                              << solved_fields.az[i].real() << ","
                              << solved_fields.az[i].imag() << ","
                              << solved_fields.phi[i].real() << ","
                              << solved_fields.phi[i].imag() << ","
                              << abs_a[i] << ","
                              << abs_phi[i] << ","
                              << abs_e[i] << ","
                              << abs_j[i] << ","
                              << joule[i] << "\n";
            }
        }

        std::cout << "[em-aphi-laplace-team7-stl-em-only] " << solver_message
                  << ", solver_backend=" << active_solver->LastSolverBackend()
                  << ", solver_iterations=" << active_solver->LastSolverIterations()
                  << ", solver_estimated_error=" << active_solver->LastSolverEstimatedError()
                  << ", particles=" << total_particles
                  << ", plate_particles=" << plate_particles
                  << ", coil_particles=" << coil_particles
                  << ", source_particles=" << source_particles
                  << ", solved_total_joule_power=" << solved_diagnostics.total_joule_power
                  << ", solved_max_abs_j=" << solved_diagnostics.max_abs_j
                  << ", plate_mean_abs_j=" << plate_summary.mean_abs_j
                  << ", plate_mean_joule=" << plate_summary.mean_joule
                  << ", coil_mean_joule=" << coil_summary.mean_joule
                  << ", solved_linear_residual_norm=" << solved_diagnostics.linear_residual_norm
                  << ", time_geometry_s=" << time_geometry_s
                  << ", time_particle_generation_s=" << time_particle_generation_s
                  << ", time_region_setup_s=" << time_region_setup_s
                  << ", time_assemble_s=" << time_assemble_s
                  << ", time_solve_s=" << time_solve_s
                  << ", time_postprocess_s=" << time_postprocess_s
                  << ", time_output_s=" << time_output_s
                  << ", time_total_s=" << time_total_s
                  << std::endl;

        return 0;
    }

    auto &inner_structure_shape = sph_system.addShape<InnerCoilPlateUnionShape>("Team7InnerStructure");
    Real dp_air_coarse = use_multi_resolution ? dp_air_coarse_input : dp_0;
    int air_refinement_levels = use_multi_resolution ? air_refinement_levels_input : 0;
    AdaptiveAroundInnerStructure body_adaptation(
        dp_air_coarse, 1.15, 1.0, air_refinement_levels, &inner_structure_shape);
    auto &physical_body = sph_system.addAdaptiveBody<SolidBody, AdaptiveAroundInnerStructure>(
        body_adaptation,
        makeShared<Team7AirBoxShape>("Team7STLEMOnlyBox", team7_air_box.air_lower, team7_air_box.air_upper));
    physical_body.defineMaterial<Solid>();
    physical_body.defineBodyLevelSetShape();
    physical_body.generateParticles<BaseParticles, Lattice>();

    UniquePtr<BaseInnerRelation> inner_relation_ptr;
    if (use_multi_resolution)
    {
        inner_relation_ptr = makeUnique<AdaptiveInnerRelation>(physical_body);
    }
    else
    {
        inner_relation_ptr = makeUnique<InnerRelation>(physical_body);
    }
    BaseInnerRelation &inner_relation = *inner_relation_ptr;
    BodyStatesRecordingToVtp write_states(sph_system);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    const auto time_particles_ready = Clock::now();

    BaseParticles &particles = physical_body.getBaseParticles();
    size_t total_particles = particles.TotalRealParticles();
    Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *vtp_a_real = particles.registerStateVariableData<Vecd>("SolvedAReal", ZeroData<Vecd>::value);
    Vecd *vtp_a_imag = particles.registerStateVariableData<Vecd>("SolvedAImag", ZeroData<Vecd>::value);
    Real *vtp_phi_real = particles.registerStateVariableData<Real>("SolvedPhiReal", Real(0));
    Real *vtp_phi_imag = particles.registerStateVariableData<Real>("SolvedPhiImag", Real(0));
    Vecd *vtp_e_real = particles.registerStateVariableData<Vecd>("SolvedEReal", ZeroData<Vecd>::value);
    Vecd *vtp_e_imag = particles.registerStateVariableData<Vecd>("SolvedEImag", ZeroData<Vecd>::value);
    Vecd *vtp_j_real = particles.registerStateVariableData<Vecd>("SolvedJReal", ZeroData<Vecd>::value);
    Vecd *vtp_j_imag = particles.registerStateVariableData<Vecd>("SolvedJImag", ZeroData<Vecd>::value);
    Real *vtp_abs_a = particles.registerStateVariableData<Real>("SolvedAbsA", Real(0));
    Real *vtp_abs_phi = particles.registerStateVariableData<Real>("SolvedAbsPhi", Real(0));
    Real *vtp_abs_e = particles.registerStateVariableData<Real>("SolvedAbsE", Real(0));
    Real *vtp_abs_j = particles.registerStateVariableData<Real>("SolvedAbsJ", Real(0));
    Real *vtp_joule = particles.registerStateVariableData<Real>("SolvedJouleHeat", Real(0));
    Real *vtp_is_plate = particles.registerStateVariableData<Real>("RegionPlate", Real(0));
    Real *vtp_is_coil = particles.registerStateVariableData<Real>("RegionCoil", Real(0));
    Real *vtp_is_source = particles.registerStateVariableData<Real>("RegionSource", Real(0));
    Real *vtp_source_jz = particles.registerStateVariableData<Real>("SourceJz", Real(0));
    Real *vtp_h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");

    write_states.addToWrite<Vecd>(physical_body, "SolvedAReal");
    write_states.addToWrite<Vecd>(physical_body, "SolvedAImag");
    write_states.addToWrite<Real>(physical_body, "SolvedPhiReal");
    write_states.addToWrite<Real>(physical_body, "SolvedPhiImag");
    write_states.addToWrite<Vecd>(physical_body, "SolvedEReal");
    write_states.addToWrite<Vecd>(physical_body, "SolvedEImag");
    write_states.addToWrite<Vecd>(physical_body, "SolvedJReal");
    write_states.addToWrite<Vecd>(physical_body, "SolvedJImag");
    write_states.addToWrite<Real>(physical_body, "SolvedAbsA");
    write_states.addToWrite<Real>(physical_body, "SolvedAbsPhi");
    write_states.addToWrite<Real>(physical_body, "SolvedAbsE");
    write_states.addToWrite<Real>(physical_body, "SolvedAbsJ");
    write_states.addToWrite<Real>(physical_body, "SolvedJouleHeat");
    write_states.addToWrite<Real>(physical_body, "RegionPlate");
    write_states.addToWrite<Real>(physical_body, "RegionCoil");
    write_states.addToWrite<Real>(physical_body, "RegionSource");
    write_states.addToWrite<Real>(physical_body, "SourceJz");
    write_states.addToWrite<Real>(physical_body, "SmoothingLengthRatio");

    StdVec<Real> sigma_physical(total_particles, sigma_air);
    StdVec<Real> sigma_em(total_particles, use_magnetic_only_air ? 0.0 : sigma_air);
    StdVec<Real> sigma_solve(total_particles, use_magnetic_only_air ? sigma_solve_regularization : sigma_air);
    StdVec<Real> nu(total_particles, nu_air);
    StdVec<bool> is_plate(total_particles, false);
    StdVec<bool> is_coil(total_particles, false);
    StdVec<bool> is_source(total_particles, false);
    StdVec<Vecd> source_current_density(total_particles, ZeroData<Vecd>::value);

    size_t plate_particles = 0;
    size_t coil_particles = 0;
    size_t air_particles = 0;
    size_t source_particles = 0;
    Real coil_volume_total_geom = 0.0;
    Real coil_radius_volume_moment = 0.0;

    for (size_t i = 0; i != total_particles; ++i)
    {
        const Vec3d position = positions[i];
        bool plate_region = plate_mesh.findSignedDistance(position) <= 0.0;
        bool coil_region = coil_mesh.findSignedDistance(position) <= 0.0;

        is_plate[i] = plate_region;
        is_coil[i] = coil_region;

        if (plate_region)
        {
            sigma_physical[i] = sigma_plate;
            sigma_em[i] = sigma_plate;
            sigma_solve[i] = sigma_plate;
            nu[i] = nu_plate;
            plate_particles++;
        }
        else if (coil_region)
        {
            sigma_physical[i] = sigma_coil;
            sigma_em[i] = use_magnetic_only_coil ? 0.0 : sigma_coil;
            sigma_solve[i] = use_magnetic_only_coil ? sigma_solve_regularization : sigma_coil;
            nu[i] = nu_coil;
            coil_particles++;
            coil_volume_total_geom += volumetric_measure[i];
            if (use_circular_coil_source)
            {
                Vecd relative = position - coil_center;
                Vecd axis_used = NormalizeVectorOrDefault(coil_source_axis, Vecd(0.0, 0.0, 1.0));
                Vecd radial = relative - relative.dot(axis_used) * axis_used;
                coil_radius_volume_moment += volumetric_measure[i] * radial.norm();
            }
        }
        else
        {
            air_particles++;
        }
        vtp_is_plate[i] = is_plate[i] ? 1.0 : 0.0;
        vtp_is_coil[i] = is_coil[i] ? 1.0 : 0.0;
    }

    Real safe_geom_length_to_m = SMAX(geom_length_to_m, TinyReal);
    Real effective_area_geom = coil_effective_area_input_geom;
    if (!(effective_area_geom > TinyReal))
    {
        if (use_circular_coil_source)
        {
            Real mean_radius_geom =
                coil_radius_volume_moment / (coil_volume_total_geom + TinyReal);
            Real projection_length_geom = 2.0 * Pi * mean_radius_geom;
            effective_area_geom =
                coil_volume_total_geom / (projection_length_geom + TinyReal);
        }
        else
        {
            Vecd dir_used = NormalizeVectorOrDefault(fallback_source_direction, Vecd(0.0, 0.0, 1.0));
            Vecd coil_extent = coil_bounds.upper_ - coil_bounds.lower_;
            Real projection_length_geom =
                fabs(dir_used[0]) * coil_extent[0] +
                fabs(dir_used[1]) * coil_extent[1] +
                fabs(dir_used[2]) * coil_extent[2];
            effective_area_geom =
                coil_volume_total_geom / (projection_length_geom + TinyReal);
        }
    }
    Real effective_area_si = effective_area_geom * safe_geom_length_to_m * safe_geom_length_to_m;
    Real current_peak = SMAX(coil_current_peak_factor, static_cast<Real>(1.0)) *
                        SMAX(coil_current_rms, static_cast<Real>(0.0));
    Real source_current_density_magnitude =
        coil_turns * current_peak / (effective_area_si + TinyReal);

    std::cout << std::scientific << std::setprecision(6)
              << "[em-aphi-laplace-team7-stl-em-only-source] coil_volume_total_geom=" << coil_volume_total_geom
              << ", effective_area_geom=" << effective_area_geom
              << ", effective_area_si=" << effective_area_si
              << ", source_current_density_magnitude=" << source_current_density_magnitude
              << std::endl;

    for (size_t i = 0; i != total_particles; ++i)
    {
        Vecd js = is_coil[i]
                      ? source_profile_current_density(positions[i], coil_center, coil_source_axis,
                                                       fallback_source_direction, source_current_density_magnitude,
                                                       use_circular_coil_source)
                      : ZeroData<Vecd>::value;
        source_current_density[i] = js;
        is_source[i] = js.squaredNorm() > TinyReal;
        vtp_is_source[i] = is_source[i] ? 1.0 : 0.0;
        vtp_source_jz[i] = js[2];
        if (is_source[i])
        {
            source_particles++;
        }
    }

    electromagnetics::LaplaceStructuredAPhiBoundaryCondition boundary_condition;
    boundary_condition.is_dirichlet.resize(total_particles, false);
    boundary_condition.ax.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.ay.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.az.resize(total_particles, Complex(0.0, 0.0));
    boundary_condition.phi.resize(total_particles, Complex(0.0, 0.0));

    size_t reference_phi_index = 0;
    for (size_t i = 0; i != total_particles; ++i)
    {
        boundary_condition.is_dirichlet[i] = is_boundary_particle(positions[i], team7_air_box);
        if (is_plate[i] && !boundary_condition.is_dirichlet[i])
        {
            reference_phi_index = i;
        }
    }

    electromagnetics::LaplaceStructuredAPhiParameters parameters;
    parameters.angular_frequency = 2.0 * Pi * frequency_hz;
    parameters.geom_length_to_m = geom_length_to_m;
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
    parameters.contact_gradient_scale = contact_gradient_scale;
    parameters.contact_diffusion_scale = contact_diffusion_scale;

    electromagnetics::LaplaceStructuredAPhiDiscreteView single_body_view{
        total_particles,
        physical_body.getSPHAdaptation().ReferenceSmoothingLength(),
        positions,
        volumetric_measure,
        vtp_h_ratio,
        &inner_relation.inner_configuration_,
        nullptr};
    electromagnetics::LaplaceStructuredAPhiEigenSolver solver(single_body_view, parameters);
    electromagnetics::LaplaceStructuredAPhiEigenSolver *active_solver = &solver;
    std::unique_ptr<electromagnetics::LaplaceStructuredAPhiEigenSolver> fallback_solver;
    electromagnetics::VecC rhs(static_cast<int>(4 * total_particles));
    rhs.setZero();
    for (size_t i = 0; i != total_particles; ++i)
    {
        rhs[solver.idxAx(i)] = Complex(source_current_density[i][0], 0.0);
        rhs[solver.idxAy(i)] = Complex(source_current_density[i][1], 0.0);
        rhs[solver.idxAz(i)] = Complex(source_current_density[i][2], 0.0);
    }
    const auto time_setup_ready = Clock::now();

    solver.assembleSystem(sigma_solve, nu, rhs, boundary_condition);
    const auto time_assembled = Clock::now();

    electromagnetics::LaplaceStructuredAPhiFields solved_fields;
    std::string solver_message;
    bool solve_success = solver.solve(solved_fields, solver_message);
    if (!solve_success)
    {
        std::cerr << "[em-aphi-laplace-team7-stl-em-only] primary solve failed"
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
                single_body_view, fallback_parameters);
            fallback_solver->assembleSystem(sigma_solve, nu, rhs, boundary_condition);

            solve_success = fallback_solver->solve(solved_fields, solver_message);
            if (!solve_success)
            {
                std::cerr << "[em-aphi-laplace-team7-stl-em-only] fallback solve failed"
                          << ", backend=" << fallback_solver->LastSolverBackend()
                          << ", iterations=" << fallback_solver->LastSolverIterations()
                          << ", estimated_error=" << fallback_solver->LastSolverEstimatedError()
                          << ", message=" << solver_message
                          << std::endl;
                return 1;
            }

            std::cout << "[em-aphi-laplace-team7-stl-em-only] fallback solve succeeded"
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
    const auto time_solved = Clock::now();

    electromagnetics::LaplaceStructuredAPhiDiagnostics solved_diagnostics =
        active_solver->computeDiagnostics(sigma_em, nu, solved_fields);

    PhysicalRegionSummary plate_summary;
    PhysicalRegionSummary air_summary;
    PhysicalRegionSummary coil_summary;
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

        vtp_a_real[i] = Vecd(solved_fields.ax[i].real(), solved_fields.ay[i].real(), solved_fields.az[i].real());
        vtp_a_imag[i] = Vecd(solved_fields.ax[i].imag(), solved_fields.ay[i].imag(), solved_fields.az[i].imag());
        vtp_phi_real[i] = solved_fields.phi[i].real();
        vtp_phi_imag[i] = solved_fields.phi[i].imag();
        vtp_abs_a[i] = abs_a;
        vtp_abs_phi[i] = abs_phi;
        vtp_abs_e[i] = abs_e;
        vtp_abs_j[i] = abs_j;
        vtp_joule[i] = joule;

        if (!solved_fields.electric_field.empty())
        {
            const auto &electric_field = solved_fields.electric_field[i];
            vtp_e_real[i] = Vecd(electric_field[0].real(), electric_field[1].real(), electric_field[2].real());
            vtp_e_imag[i] = Vecd(electric_field[0].imag(), electric_field[1].imag(), electric_field[2].imag());
        }
        if (!solved_fields.current_density.empty())
        {
            const auto &current_density = solved_fields.current_density[i];
            vtp_j_real[i] = Vecd(current_density[0].real(), current_density[1].real(), current_density[2].real());
            vtp_j_imag[i] = Vecd(current_density[0].imag(), current_density[1].imag(), current_density[2].imag());
        }

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

        if (is_plate[i])
        {
            accumulate_summary(plate_summary);
        }
        else if (is_coil[i])
        {
            accumulate_summary(coil_summary);
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

    finalize_region_summary(plate_summary);
    finalize_region_summary(air_summary);
    finalize_region_summary(coil_summary);
    finalize_region_summary(source_summary);
    const auto time_postprocess_ready = Clock::now();

    if (write_vtp)
    {
        physical_body.setNewlyUpdated();
        write_states.writeToFile(0.0);
    }
    const auto time_output_ready = Clock::now();

    const double time_geometry_s = seconds_since(time_begin, time_geometry_ready);
    const double time_particle_generation_s = seconds_since(time_geometry_ready, time_particles_ready);
    const double time_region_setup_s = seconds_since(time_particles_ready, time_setup_ready);
    const double time_assemble_s = seconds_since(time_setup_ready, time_assembled);
    const double time_solve_s = seconds_since(time_assembled, time_solved);
    const double time_postprocess_s = seconds_since(time_solved, time_postprocess_ready);
    const double time_output_s = seconds_since(time_postprocess_ready, time_output_ready);
    const double time_total_s = seconds_since(time_begin, time_output_ready);

    std::ofstream summary_file(
        io_environment.OutputFolder() + "/em_aphi_laplace_team7_stl_em_only_summary.csv",
        std::ios::out | std::ios::trunc);
    summary_file << std::setprecision(12);
    summary_file << "solve_success,particles,plate_particles,coil_particles,air_particles,source_particles,"
                 << "dp,air_box_lower_x,air_box_lower_y,air_box_lower_z,air_box_upper_x,air_box_upper_y,air_box_upper_z,"
                 << "used_small_air_box,used_tight_stl_air_box,tight_stl_margin,sigma_air,sigma_plate,sigma_coil,sigma_solve_regularization,nu_air,nu_plate,nu_coil,frequency_hz,coil_current_rms,coil_turns,coil_current_peak_factor,coil_effective_area_input_geom,coil_effective_area_used_geom,coil_effective_area_used_m2,geom_length_to_m,use_circular_coil_source,use_three_body_layout,use_multi_resolution,dp_air_coarse,air_refinement_levels,use_magnetic_only_air,use_magnetic_only_coil,coil_source_axis_x,coil_source_axis_y,coil_source_axis_z,"
                 << "solver_backend,solver_iterations,solver_estimated_error,"
                 << "used_phi_block_scale,used_diagonal_equilibration_iterations,"
                 << "solved_div_a_l2,solved_current_continuity_l2,solved_phi_equation_residual_l2,"
                 << "solved_physical_current_continuity_l2,solved_relative_current_continuity_l2,"
                 << "solved_linear_residual_norm,solved_total_joule_power,solved_gauge_penalty_energy,solved_magnetic_energy,"
                 << "time_geometry_s,time_particle_generation_s,time_region_setup_s,time_assemble_s,time_solve_s,time_postprocess_s,time_output_s,time_total_s,"
                 << "solved_max_abs_a,solved_max_abs_e,solved_max_abs_j,"
                 << "plate_mean_abs_a,plate_max_abs_a,plate_mean_abs_phi,plate_max_abs_phi,"
                 << "plate_mean_abs_e,plate_max_abs_e,plate_mean_abs_j,plate_max_abs_j,plate_mean_joule,plate_max_joule,"
                 << "coil_mean_abs_a,coil_max_abs_a,coil_mean_abs_phi,coil_max_abs_phi,"
                 << "coil_mean_abs_e,coil_max_abs_e,coil_mean_abs_j,coil_max_abs_j,coil_mean_joule,coil_max_joule,"
                 << "air_mean_abs_a,air_max_abs_a,air_mean_abs_phi,air_max_abs_phi,"
                 << "air_mean_abs_e,air_max_abs_e,air_mean_abs_j,air_max_abs_j,air_mean_joule,air_max_joule,"
                 << "source_mean_abs_a,source_max_abs_a,source_mean_abs_phi,source_max_abs_phi,"
                 << "source_mean_abs_e,source_max_abs_e,source_mean_abs_j,source_max_abs_j,source_mean_joule,source_max_joule\n";

    summary_file << (solve_success ? 1 : 0) << ","
                 << total_particles << ","
                 << plate_particles << ","
                 << coil_particles << ","
                 << air_particles << ","
                 << source_particles << ","
                 << dp_0 << ","
                 << team7_air_box.air_lower[0] << ","
                 << team7_air_box.air_lower[1] << ","
                 << team7_air_box.air_lower[2] << ","
                 << team7_air_box.air_upper[0] << ","
                 << team7_air_box.air_upper[1] << ","
                 << team7_air_box.air_upper[2] << ","
                 << (team7_air_box.used_small_air_box ? 1 : 0) << ","
                 << (team7_air_box.used_tight_stl_air_box ? 1 : 0) << ","
                 << team7_air_box.tight_stl_margin << ","
                 << sigma_air << ","
                 << sigma_plate << ","
                 << sigma_coil << ","
                 << sigma_solve_regularization << ","
                 << nu_air << ","
                 << nu_plate << ","
                 << nu_coil << ","
                 << frequency_hz << ","
                 << coil_current_rms << ","
                 << coil_turns << ","
                 << coil_current_peak_factor << ","
                 << coil_effective_area_input_geom << ","
                 << effective_area_geom << ","
                 << effective_area_si << ","
                 << geom_length_to_m << ","
                 << (use_circular_coil_source ? 1 : 0) << ","
                 << (use_three_body_layout ? 1 : 0) << ","
                 << (use_multi_resolution ? 1 : 0) << ","
                 << dp_air_coarse << ","
                 << air_refinement_levels << ","
                 << (use_magnetic_only_air ? 1 : 0) << ","
                 << (use_magnetic_only_coil ? 1 : 0) << ","
                 << coil_source_axis[0] << ","
                 << coil_source_axis[1] << ","
                 << coil_source_axis[2] << ","
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
                 << time_geometry_s << ","
                 << time_particle_generation_s << ","
                 << time_region_setup_s << ","
                 << time_assemble_s << ","
                 << time_solve_s << ","
                 << time_postprocess_s << ","
                 << time_output_s << ","
                 << time_total_s << ","
                 << solved_diagnostics.max_abs_a << ","
                 << solved_diagnostics.max_abs_e << ","
                 << solved_diagnostics.max_abs_j << ","
                 << plate_summary.mean_abs_a << ","
                 << plate_summary.max_abs_a << ","
                 << plate_summary.mean_abs_phi << ","
                 << plate_summary.max_abs_phi << ","
                 << plate_summary.mean_abs_e << ","
                 << plate_summary.max_abs_e << ","
                 << plate_summary.mean_abs_j << ","
                 << plate_summary.max_abs_j << ","
                 << plate_summary.mean_joule << ","
                 << plate_summary.max_joule << ","
                 << coil_summary.mean_abs_a << ","
                 << coil_summary.max_abs_a << ","
                 << coil_summary.mean_abs_phi << ","
                 << coil_summary.max_abs_phi << ","
                 << coil_summary.mean_abs_e << ","
                 << coil_summary.max_abs_e << ","
                 << coil_summary.mean_abs_j << ","
                 << coil_summary.max_abs_j << ","
                 << coil_summary.mean_joule << ","
                 << coil_summary.max_joule << ","
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
            io_environment.OutputFolder() + "/em_aphi_laplace_team7_stl_em_only_particles.csv",
            std::ios::out | std::ios::trunc);
        particle_file << std::setprecision(12);
        particle_file << "particle_id,x,y,z,is_boundary,is_plate,is_coil,is_source,smoothing_length_ratio,sigma_physical,sigma_em,sigma_solve,nu,source_jx,source_jy,source_jz,"
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
                          << (is_plate[i] ? 1 : 0) << ","
                          << (is_coil[i] ? 1 : 0) << ","
                          << (is_source[i] ? 1 : 0) << ","
                          << vtp_h_ratio[i] << ","
                          << sigma_physical[i] << ","
                          << sigma_em[i] << ","
                          << sigma_solve[i] << ","
                          << nu[i] << ","
                          << source_current_density[i][0] << ","
                          << source_current_density[i][1] << ","
                          << source_current_density[i][2] << ","
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
    std::cout << "[em-aphi-laplace-team7-stl-em-only] " << solver_message
              << ", solver_backend=" << active_solver->LastSolverBackend()
              << ", solver_iterations=" << active_solver->LastSolverIterations()
              << ", solver_estimated_error=" << active_solver->LastSolverEstimatedError()
              << ", particles=" << total_particles
              << ", plate_particles=" << plate_particles
              << ", coil_particles=" << coil_particles
              << ", source_particles=" << source_particles
              << ", solved_total_joule_power=" << solved_diagnostics.total_joule_power
              << ", solved_max_abs_j=" << solved_diagnostics.max_abs_j
              << ", plate_mean_abs_j=" << plate_summary.mean_abs_j
              << ", plate_mean_joule=" << plate_summary.mean_joule
              << ", coil_mean_joule=" << coil_summary.mean_joule
              << ", solved_linear_residual_norm=" << solved_diagnostics.linear_residual_norm
              << ", time_geometry_s=" << time_geometry_s
              << ", time_particle_generation_s=" << time_particle_generation_s
              << ", time_region_setup_s=" << time_region_setup_s
              << ", time_assemble_s=" << time_assemble_s
              << ", time_solve_s=" << time_solve_s
              << ", time_postprocess_s=" << time_postprocess_s
              << ", time_output_s=" << time_output_s
              << ", time_total_s=" << time_total_s
              << std::endl;

    return 0;
}
