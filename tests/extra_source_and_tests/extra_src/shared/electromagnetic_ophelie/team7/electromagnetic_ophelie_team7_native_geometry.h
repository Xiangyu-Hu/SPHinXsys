#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_NATIVE_GEOMETRY_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_NATIVE_GEOMETRY_H

#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_relaxation.h"
#include "sphinxsys.h"

#include "io_environment.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

namespace fs = std::filesystem;

/** TEAM Problem 7 asymmetric conductor + coil (COMSOL / muFEM geometry in mm STL). */
struct OphelieTeam7NativeMesh
{
    /** Multiply STL vertex coordinates (mm) to SI positions (m). */
    Real stl_scale_to_meter_ = 1.0e-3;
    Vec3d air_lower_small_{-50.0, -50.0, -50.0};
    Vec3d air_upper_small_{350.0, 350.0, 200.0};
    Vec3d air_lower_standard_{-1353.0, -1353.0, -300.0};
    Vec3d air_upper_standard_{1647.0, 1647.0, 449.0};
    Real dp_ref_mm_ = 6.0;
    Real dp_coil_mm_ = 6.0;
    Real dp_plate_mm_ = 6.0;
    bool use_small_air_box_ = true;
    /** Relax air SPH particles (TEAM7 conductive bodies still use CoilSourceBody + PlateBody only in EM). */
    bool relax_air_particles_ = true;
    size_t relaxation_steps_ = 1000;
    size_t relaxation_log_every_ = 100;
    size_t relaxation_vtp_every_ = 100;
    /** COMSOL TEAM7 benchmark: 2742 turns, 1 A/turn (see COMSOL blog). */
    Real team7_coil_turns_ = 2742.0;
    Real team7_coil_current_per_turn_ = 1.0;
};

/** COMSOL / muFEM TEAM7 STL mesh bounds in mm (triangle vertices, pre-relax). */
struct Team7NativeStlMeshBBoxMm
{
    Vec3d coil_lower_{94.0, 0.0, 49.0};
    Vec3d coil_upper_{294.0, 200.0, 149.0};
    Vec3d plate_lower_{0.0, 0.0, 0.0};
    Vec3d plate_upper_{294.0, 294.0, 19.0};
};

/** TEAM7 / muFEM A1–B1 Bz probe polyline in COMSOL global mm coordinates. */
struct Team7ReferenceProbeLineMm
{
    static constexpr Real y_mm = 72.0;
    static constexpr Real z_mm = 34.0;
    static constexpr Real x_start_mm = 0.0;
    static constexpr Real x_end_mm = 288.0;
};

struct OphelieTeam7NativeDerivedGeometry
{
    Vecd coil_center_ = Vecd::Zero();
    Vecd plate_center_ = Vecd::Zero();
    BoundingBoxd coil_bbox_;
    BoundingBoxd plate_bbox_;
    Real coil_current_cross_section_m2_ = 0.0;
    Real coil_volume_m3_ = 0.0;
    Real coil_mean_radius_m_ = 0.0;
    Real coil_max_xy_radius_m_ = 0.0;
    Real plate_characteristic_length_m_ = 0.0;
};

inline Vec3d vec3FromVecdMm(const Vecd &v_m, Real inv_mm_to_m)
{
    return Vec3d(v_m[0] * inv_mm_to_m, v_m[1] * inv_mm_to_m, v_m[2] * inv_mm_to_m);
}

/** A_cross = sum_j Vol_j / (2*pi*r_j), r_j = xy distance to coil_center (TEAM7 STL/reload). */
inline Real estimateCoilCurrentCrossSectionFromParticles(BaseParticles &particles, const Vecd &coil_center,
                                                         Real &coil_volume_m3, Real &mean_radius_m)
{
    coil_volume_m3 = 0.0;
    mean_radius_m = 0.0;
    Real cross_section = 0.0;
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        return TinyReal;
    }
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - coil_center[0];
        const Real dy = pos[i][1] - coil_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy) + TinyReal;
        cross_section += vol[i] / (2.0 * Pi * r);
        coil_volume_m3 += vol[i];
        mean_radius_m += r;
    }
    mean_radius_m /= static_cast<Real>(n);
    return std::max(cross_section, TinyReal);
}

/** Max xy distance from coil_center to coil particles (for outer-shell J weighting). */
inline Real estimateCoilMaxXYRadiusFromParticles(BaseParticles &particles, const Vecd &coil_center)
{
    Real r_max = 0.0;
    const size_t n = particles.TotalRealParticles();
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - coil_center[0];
        const Real dy = pos[i][1] - coil_center[1];
        r_max = std::max(r_max, std::sqrt(dx * dx + dy * dy));
    }
    return std::max(r_max, TinyReal);
}

/** A_cross counting only particles with r_xy >= fraction * r_max (current on coil outer pack). */
inline Real estimateCoilCurrentCrossSectionOuterShell(BaseParticles &particles, const Vecd &coil_center, Real r_max_xy,
                                                      Real radius_fraction, Real &coil_volume_shell_m3)
{
    coil_volume_shell_m3 = 0.0;
    Real cross_section = 0.0;
    const Real r_min = radius_fraction * r_max_xy;
    const size_t n = particles.TotalRealParticles();
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - coil_center[0];
        const Real dy = pos[i][1] - coil_center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        if (r >= r_min)
        {
            cross_section += vol[i] / (2.0 * Pi * r);
            coil_volume_shell_m3 += vol[i];
        }
    }
    return std::max(cross_section, TinyReal);
}

inline void printTeam7NativeStlMeshBBoxAudit(const Team7NativeStlMeshBBoxMm &mesh_bbox)
{
    std::cout << "[ophelie] TEAM7 STL mesh bbox (mm, from triangle vertices):" << std::endl;
    std::cout << "[ophelie]   coil_mesh lower=" << mesh_bbox.coil_lower_.transpose()
              << " upper=" << mesh_bbox.coil_upper_.transpose() << std::endl;
    std::cout << "[ophelie]   plate_mesh lower=" << mesh_bbox.plate_lower_.transpose()
              << " upper=" << mesh_bbox.plate_upper_.transpose() << std::endl;
}

inline std::string resolveTeam7NativeStlPath(const std::string &preferred, const StdVec<std::string> &fallbacks)
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

inline std::string team7NativeCoilStlPath()
{
    return resolveTeam7NativeStlPath(
        "./input/coil.stl",
        {
            "./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/input/coil.stl",
            "../data/coil.stl",
            "../../particle_generation_em/data/coil.stl",
            "../particle_generation_em/data/coil.stl",
        });
}

inline std::string team7NativePlateStlPath()
{
    return resolveTeam7NativeStlPath(
        "./input/plate.stl",
        {
            "./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/input/plate.stl",
            "../data/plate.stl",
            "../../particle_generation_em/data/plate.stl",
            "../particle_generation_em/data/plate.stl",
        });
}

class OphelieTeam7NativeCoilShape : public ComplexShape
{
  public:
    OphelieTeam7NativeCoilShape(const std::string &shape_name, const std::string &coil_stl_path, Real stl_scale)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(coil_stl_path, Vec3d::Zero(), stl_scale);
    }
};

class OphelieTeam7NativePlateShape : public ComplexShape
{
  public:
    OphelieTeam7NativePlateShape(const std::string &shape_name, const std::string &plate_stl_path, Real stl_scale)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(plate_stl_path, Vec3d::Zero(), stl_scale);
    }
};

class OphelieTeam7NativeAirShape : public ComplexShape
{
  public:
    OphelieTeam7NativeAirShape(const std::string &shape_name, const Vec3d &box_lower_mm, const Vec3d &box_upper_mm,
                               const std::string &coil_stl_path, const std::string &plate_stl_path, Real stl_scale)
        : ComplexShape(shape_name)
    {
        const Vecd halfsize = 0.5 * stl_scale * (box_upper_mm - box_lower_mm);
        const Vecd center = 0.5 * stl_scale * (box_lower_mm + box_upper_mm);
        add<GeometricShapeBox>(Transform(center), halfsize, "OuterBoundary");
        subtract<TriangleMeshShapeSTL>(coil_stl_path, Vec3d::Zero(), stl_scale, "Coil");
        subtract<TriangleMeshShapeSTL>(plate_stl_path, Vec3d::Zero(), stl_scale, "Plate");
    }
};

inline Vecd particleCentroid(BaseParticles &particles)
{
    const size_t n = particles.TotalRealParticles();
    Vecd sum = Vecd::Zero();
    if (n == 0)
    {
        return sum;
    }
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    for (size_t i = 0; i < n; ++i)
    {
        sum += pos[i];
    }
    return sum / static_cast<Real>(n);
}

inline BoundingBoxd particleBoundingBox(BaseParticles &particles)
{
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        return BoundingBoxd(Vecd::Zero(), Vecd::Zero());
    }
    Vecd lower = Vecd::Constant(std::numeric_limits<Real>::infinity());
    Vecd upper = Vecd::Constant(-std::numeric_limits<Real>::infinity());
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    for (size_t i = 0; i < n; ++i)
    {
        lower = lower.cwiseMin(pos[i]);
        upper = upper.cwiseMax(pos[i]);
    }
    return BoundingBoxd(lower, upper);
}

inline OphelieTeam7NativeDerivedGeometry deriveTeam7NativeGeometry(SolidBody &coil_body, SolidBody &plate_body)
{
    OphelieTeam7NativeDerivedGeometry derived;
    derived.coil_center_ = particleCentroid(coil_body.getBaseParticles());
    derived.plate_center_ = particleCentroid(plate_body.getBaseParticles());
    derived.coil_bbox_ = particleBoundingBox(coil_body.getBaseParticles());
    derived.plate_bbox_ = particleBoundingBox(plate_body.getBaseParticles());

    derived.coil_current_cross_section_m2_ = estimateCoilCurrentCrossSectionFromParticles(
        coil_body.getBaseParticles(), derived.coil_center_, derived.coil_volume_m3_, derived.coil_mean_radius_m_);
    derived.coil_max_xy_radius_m_ =
        estimateCoilMaxXYRadiusFromParticles(coil_body.getBaseParticles(), derived.coil_center_);
    const Vecd plate_extent = derived.plate_bbox_.upper_ - derived.plate_bbox_.lower_;
    derived.plate_characteristic_length_m_ = std::max(plate_extent.maxCoeff(), TinyReal);
    return derived;
}

/** Log coil/plate bbox (mm) and whether the A1–B1 probe line sits in the air gap between bodies. */
inline void printTeam7NativeGeometryProbeAudit(SolidBody &coil_body, SolidBody &plate_body,
                                               const OphelieTeam7NativeDerivedGeometry &derived,
                                               const OphelieParameters &params, const OphelieTeam7NativeMesh &mesh,
                                               Real mm_to_m = 1.0e-3)
{
    const Real inv_mm = 1.0 / mm_to_m;
    const BoundingBoxd coil_bbox = particleBoundingBox(coil_body.getBaseParticles());
    const BoundingBoxd plate_bbox = particleBoundingBox(plate_body.getBaseParticles());
    const Vec3d coil_lo = vec3FromVecdMm(coil_bbox.lower_, inv_mm);
    const Vec3d coil_hi = vec3FromVecdMm(coil_bbox.upper_, inv_mm);
    const Vec3d plate_lo = vec3FromVecdMm(plate_bbox.lower_, inv_mm);
    const Vec3d plate_hi = vec3FromVecdMm(plate_bbox.upper_, inv_mm);
    const Vec3d coil_center_mm = vec3FromVecdMm(derived.coil_center_, inv_mm);
    const Vec3d plate_center_mm = vec3FromVecdMm(derived.plate_center_, inv_mm);

    std::cout << "[ophelie] TEAM7 geometry audit (COMSOL mm frame, stl_scale=" << mesh.stl_scale_to_meter_ << ")"
              << std::endl;
    std::cout << "[ophelie]   coil_bbox_mm lower=" << coil_lo.transpose() << " upper=" << coil_hi.transpose()
              << std::endl;
    std::cout << "[ophelie]   plate_bbox_mm lower=" << plate_lo.transpose() << " upper=" << plate_hi.transpose()
              << std::endl;
    std::cout << "[ophelie]   coil_center_mm=" << coil_center_mm.transpose()
              << " plate_center_mm=" << plate_center_mm.transpose() << std::endl;
    std::cout << "[ophelie]   probe_A1_B1 y_mm=" << Team7ReferenceProbeLineMm::y_mm
              << " z_mm=" << Team7ReferenceProbeLineMm::z_mm << " x_mm=[" << Team7ReferenceProbeLineMm::x_start_mm
              << "," << Team7ReferenceProbeLineMm::x_end_mm << "]" << std::endl;
    std::cout << "[ophelie]   coil_turns=" << mesh.team7_coil_turns_
              << " I_per_turn=" << mesh.team7_coil_current_per_turn_
              << " ampere_turns=" << mesh.team7_coil_turns_ * mesh.team7_coil_current_per_turn_
              << " J0=" << params.coil_j0_override_ << " softening_m=" << params.softening_length_ << std::endl;

    const bool probe_y_in_coil_span =
        Team7ReferenceProbeLineMm::y_mm >= coil_lo[1] && Team7ReferenceProbeLineMm::y_mm <= coil_hi[1];
    const bool probe_y_in_plate_span =
        Team7ReferenceProbeLineMm::y_mm >= plate_lo[1] && Team7ReferenceProbeLineMm::y_mm <= plate_hi[1];
    const bool probe_x_in_domain = Team7ReferenceProbeLineMm::x_end_mm >= coil_lo[0] &&
                                 Team7ReferenceProbeLineMm::x_start_mm <= plate_hi[0];
    const Real dz_above_plate_top_mm = Team7ReferenceProbeLineMm::z_mm - plate_hi[2];
    const Real dz_below_coil_bottom_mm = coil_lo[2] - Team7ReferenceProbeLineMm::z_mm;
    const bool probe_in_air_gap =
        Team7ReferenceProbeLineMm::z_mm > plate_hi[2] && Team7ReferenceProbeLineMm::z_mm < coil_lo[2];
    const bool probe_left_of_coil_x = Team7ReferenceProbeLineMm::x_start_mm < coil_lo[0];
    const bool probe_right_of_coil_x = Team7ReferenceProbeLineMm::x_end_mm > coil_hi[0];

    std::cout << "[ophelie]   probe_y_in_coil_y_span=" << (probe_y_in_coil_span ? 1 : 0)
              << " probe_y_in_plate_y_span=" << (probe_y_in_plate_span ? 1 : 0)
              << " probe_x_in_domain=" << (probe_x_in_domain ? 1 : 0) << std::endl;
    std::cout << "[ophelie]   probe_z_above_plate_top_mm=" << dz_above_plate_top_mm
              << " probe_z_below_coil_bottom_mm=" << dz_below_coil_bottom_mm << std::endl;
    std::cout << "[ophelie]   probe_z_in_air_gap_between_plate_top_and_coil_bottom="
              << (probe_in_air_gap ? 1 : 0) << std::endl;
    std::cout << "[ophelie]   probe_x_left_of_coil_xmin=" << (probe_left_of_coil_x ? 1 : 0)
              << " probe_x_right_of_coil_xmax=" << (probe_right_of_coil_x ? 1 : 0)
              << " coil_x_span_mm=[" << coil_lo[0] << "," << coil_hi[0] << "]" << std::endl;
    if (!probe_in_air_gap)
    {
        std::cout << "[ophelie]   WARNING: A1-B1 z may not lie in the plate-top / coil-bottom air gap."
                  << std::endl;
    }
    if (probe_left_of_coil_x)
    {
        std::cout << "[ophelie]   NOTE: probe line starts left of coil x-min; left-side Bz uses far-field only."
                  << std::endl;
    }
}

inline void applyTeam7NativeParameters(OphelieParameters &params, const OphelieTeam7NativeMesh &mesh,
                                       const OphelieTeam7NativeDerivedGeometry &derived, BaseParticles &coil_particles)
{
    params.coil_center_ = derived.coil_center_;
    params.glass_center_ = derived.plate_center_;
    params.coil_max_xy_radius_m_ = derived.coil_max_xy_radius_m_;
    params.frequency_ = 50.0;
    params.current_amplitude_ = mesh.team7_coil_current_per_turn_;
    params.number_of_turns_ = mesh.team7_coil_turns_;
    Real a_cross = derived.coil_current_cross_section_m2_;
    if (params.coil_j_outer_shell_only_ && params.coil_j_outer_shell_radius_fraction_ > 0.0)
    {
        Real shell_volume = 0.0;
        a_cross = estimateCoilCurrentCrossSectionOuterShell(
            coil_particles, derived.coil_center_, derived.coil_max_xy_radius_m_,
            params.coil_j_outer_shell_radius_fraction_, shell_volume);
        std::cout << "[ophelie] outer-shell J: radius_fraction=" << params.coil_j_outer_shell_radius_fraction_
                  << " A_cross_shell=" << a_cross << " shell_volume=" << shell_volume << std::endl;
    }
    params.coil_j0_override_ = mesh.team7_coil_turns_ * mesh.team7_coil_current_per_turn_ / (a_cross + TinyReal);
    params.target_joule_power_ = 50.0e3;
    params.softening_length_ = 0.25 * mesh.dp_ref_mm_ * mesh.stl_scale_to_meter_;
    const Real mesh_dp_m = mesh.dp_ref_mm_ * mesh.stl_scale_to_meter_;
    constexpr Real french_reduced_edge_flux_dp_m = Real(0.02);
    if (params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux &&
        mesh_dp_m + TinyReal < french_reduced_edge_flux_dp_m)
    {
        const Real dp_ratio = french_reduced_edge_flux_dp_m / mesh_dp_m;
        params.pair_weight_regularization_ *= dp_ratio * dp_ratio;
    }
}

#if SPHINXSYS_USE_SYCL
/**
 * SYCL-CK relaxation for TEAM7 native STL (same pattern as particle_generation_em + mr cylinder case).
 * Air is optional; Reload.xml for EM includes only CoilSourceBody + PlateBody.
 */
inline void relaxTeam7NativeStlBodies(SPHSystem &sph_system, SolidBody &coil_body, SolidBody &plate_body,
                                      LevelSetShape &coil_level_set, LevelSetShape &plate_level_set,
                                      const OphelieTeam7NativeMesh &mesh, RealBody *air_body = nullptr,
                                      LevelSetShape *air_level_set = nullptr)
{
    const bool include_air = mesh.relax_air_particles_ && air_body != nullptr && air_level_set != nullptr;
    OphelieProgressLogger progress("relax:native-stl");
    progress.log("SYCL-CK TEAM7 native relax" + std::string(include_air ? " (coil+plate+air)" : " (coil+plate only)"));

    NearShapeSurface near_coil(coil_body);
    NearShapeSurface near_plate(plate_body);
    Inner<> coil_inner(coil_body);
    Inner<> plate_inner(plate_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    StdVec<RealBody *> relax_bodies = {&coil_body, &plate_body};
    if (include_air)
    {
        relax_bodies.push_back(air_body);
    }

    host_methods.addStateDynamics<RandomizeParticlePositionCK>(relax_bodies).exec();

    ParticleDynamicsGroup update_configuration;
    ParticleDynamicsGroup relaxation_residual;
    ParticleDynamicsGroup update_particle_position;

    auto setup_body_relax = [&](RealBody &body, Inner<> &inner, LevelSetShape &level_set, NearShapeSurface &near_surface)
    {
        update_configuration.add(&main_methods.addCellLinkedListDynamics(body));
        update_configuration.add(&main_methods.addRelationDynamics(inner));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(inner)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(body, level_set));
        update_particle_position.add(&main_methods.addStateDynamics<PositionRelaxationCK>(body));
        update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(near_surface));
    };

    setup_body_relax(coil_body, coil_inner, coil_level_set, near_coil);
    setup_body_relax(plate_body, plate_inner, plate_level_set, near_plate);

    std::unique_ptr<NearShapeSurface> near_air;
    std::unique_ptr<Inner<>> air_inner;
    if (include_air)
    {
        near_air = std::make_unique<NearShapeSurface>(*air_body);
        air_inner = std::make_unique<Inner<>>(*air_body);
        setup_body_relax(*air_body, *air_inner, *air_level_set, *near_air);
    }

    ReduceDynamicsGroup relaxation_scaling = main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(relax_bodies);

    const bool saved_state_recording = sph_system.StateRecording();
    const bool write_relax_vtp = mesh.relaxation_vtp_every_ > 0;
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(true);
    }
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);

    if (write_relax_vtp)
    {
        coil_body.setNewlyUpdated();
        plate_body.setNewlyUpdated();
        if (include_air)
        {
            air_body->setNewlyUpdated();
        }
        body_state_recorder.writeToFile(0);
        progress.log("VTP ite_0 -> " + IO::getEnvironment().OutputFolder() + "/");
    }

    for (size_t step = 0; step < mesh.relaxation_steps_; ++step)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        const Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);

        const size_t step_index = step + 1;
        if (step_index % mesh.relaxation_log_every_ == 0 || step_index == mesh.relaxation_steps_)
        {
            progress.log("step " + std::to_string(step_index) + "/" + std::to_string(mesh.relaxation_steps_));
        }
        if (write_relax_vtp && (step_index % mesh.relaxation_vtp_every_ == 0 || step_index == mesh.relaxation_steps_))
        {
            coil_body.setNewlyUpdated();
            plate_body.setNewlyUpdated();
            if (include_air)
            {
                air_body->setNewlyUpdated();
            }
            body_state_recorder.writeToFile(step_index);
            progress.log("VTP ite_" + std::to_string(step_index));
        }
    }

    coil_body.updateCellLinkedList();
    plate_body.updateCellLinkedList();
    if (include_air)
    {
        air_body->updateCellLinkedList();
    }
    if (write_relax_vtp)
    {
        sph_system.setStateRecording(saved_state_recording);
    }
    progress.finish();
}
#endif

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_NATIVE_GEOMETRY_H
