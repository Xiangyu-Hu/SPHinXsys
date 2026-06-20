/**
 * @file  particle_generation_team7.cpp
 * @brief Particle generation and relaxation for TEAM7 (native STL, SI meters).
 *
 * Coil/plate always uniform at dp_reference. Air is optional (--team7-no-air for EM-only reload).
 * Body names: CoilSourceBody, PlateBody (Ophelie EM reload convention).
 *
 * Examples (published under ./reload_cases/<case_id>/):
 *   Uniform small 3 mm, no air:  --team7-case=uniform-small-dp3 --team7-no-air
 *   Uniform legacy 6 mm:         --team7-case=uniform-legacy-dp6
 */

#include "electromagnetic_ophelie/diagnostics/aphi_team7_native_geometry_config.h"
#include "electromagnetic_ophelie/diagnostics/aphi_team7_native_team7_shapes.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "sphinxsys.h"

#include <filesystem>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using namespace SPH::electromagnetics::test;

namespace
{

class CoilShape : public ComplexShape
{
  public:
    explicit CoilShape(const std::string &shape_name, Real stl_scale) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./input/coil.stl", Vec3d::Zero(), stl_scale);
    }
};

class PlateShape : public ComplexShape
{
  public:
    explicit PlateShape(const std::string &shape_name, Real stl_scale) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./input/plate.stl", Vec3d::Zero(), stl_scale);
    }
};

class AirShape : public ComplexShape
{
  public:
    explicit AirShape(const std::string &shape_name, const Vec3d &box_lower, const Vec3d &box_upper, Real stl_scale)
        : ComplexShape(shape_name)
    {
        Vecd halfsize = 0.5 * (box_upper - box_lower);
        Vecd center = 0.5 * (box_lower + box_upper);
        Transform translation(center);
        add<GeometricShapeBox>(translation, halfsize, "OuterBoundary");
        subtract<TriangleMeshShapeSTL>("./input/coil.stl", Vec3d::Zero(), stl_scale, "Coil");
        subtract<TriangleMeshShapeSTL>("./input/plate.stl", Vec3d::Zero(), stl_scale, "Plate");
    }
};

constexpr Real k_kernel_h_spacing_ratio = 1.15;

} // namespace

int main(int ac, char *av[])
{
    AphiTeam7NativeGeometryConfig config = defaultTeam7NativeGeometryConfig();
    applyTeam7NativeGeometryEnvironment(config);
    if (parseTeam7NativeGeometryCli(ac, av, config))
    {
        return 0;
    }
    printTeam7NativeGeometryConfig(config);

    const bool include_air = config.relax_air_particles;
    const Real sph_system_resolution = team7NativeSphReferenceDp(config);
    const Real solid_global_resolution = team7NativeSolidAdaptiveGlobalResolution(config);
    const int air_refinement_levels = team7NativeAirRefinementLevels(config);
    const int solid_refinement_levels = team7NativeSolidRefinementLevels(config);

    std::cout << "particle_generation_team7 TEAM7 native relax (SI meters)"
              << " discretization=" << team7NativeDiscretizationLabel(config)
              << " reload_case_id=" << team7NativeReloadCaseId(config)
              << " dp_reference_m=" << config.dp_reference
              << " sph_system_resolution_m=" << sph_system_resolution
              << " coil_plate_refinement_levels=" << solid_refinement_levels
              << " air_refinement_levels=" << air_refinement_levels
              << " relax_air_particles=" << (include_air ? 1 : 0)
              << " air_coarsest_spacing_m=" << team7NativeAirAdaptiveGlobalResolution(config) << "\n";

    const BoundingBoxd system_domain_bounds(config.air_box_lower, config.air_box_upper);
    SPHSystem sph_system(system_domain_bounds, sph_system_resolution);

    StdVec<std::string> sph_arg_storage;
    StdVec<char *> sph_argv;
    filterSphinxsysCommandlineArgs(ac, av, sph_arg_storage, sph_argv);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_argv.size()), sph_argv.data());

    auto &inner_boundary_shape =
        sph_system.addShape<AphiTeam7NativeInnerBoundaryShape>("InnerBoundary", config.stl_scale);

    AphiTeam7NativeAdaptiveNearSurface solid_adaptation(solid_global_resolution, k_kernel_h_spacing_ratio, 1.0,
                                                        solid_refinement_levels);
    AphiTeam7NativeAdaptiveNearInnerSurface air_adaptation(
        team7NativeAirAdaptiveGlobalResolution(config), k_kernel_h_spacing_ratio, 1.0, air_refinement_levels,
        &inner_boundary_shape);

    using Team7AdaptiveCoilBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearSurface, SolidBody>;
    using Team7AdaptivePlateBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearSurface, SolidBody>;
    using Team7AdaptiveAirBody = AdaptiveBody<AphiTeam7NativeAdaptiveNearInnerSurface, RealBody>;

    Team7AdaptiveCoilBody &coil_body = sph_system.addAdaptiveBody<SolidBody, AphiTeam7NativeAdaptiveNearSurface>(
        solid_adaptation, makeShared<CoilShape>("CoilSourceBody", config.stl_scale));
    Team7AdaptivePlateBody &plate_body = sph_system.addAdaptiveBody<SolidBody, AphiTeam7NativeAdaptiveNearSurface>(
        solid_adaptation, makeShared<PlateShape>("PlateBody", config.stl_scale));

    Team7AdaptiveAirBody *air_body = nullptr;
    LevelSetShape *air_level_set = nullptr;
    if (include_air)
    {
        air_body = &sph_system.addAdaptiveBody<RealBody, AphiTeam7NativeAdaptiveNearInnerSurface>(
            air_adaptation,
            makeShared<AirShape>("Air", config.air_box_lower, config.air_box_upper, config.stl_scale));
    }

    LevelSetShape &coil_level_set =
        coil_body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet().writeLevelSet();
    LevelSetShape &plate_level_set =
        plate_body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet().writeLevelSet();
    if (include_air && air_body != nullptr)
    {
        air_level_set =
            &air_body->defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet().writeLevelSet();
        sph_system.addShape<LevelSetShape>(*air_body, inner_boundary_shape).writeLevelSet();
    }

    coil_body.generateParticles<BaseParticles, Lattice>();
    plate_body.generateParticles<BaseParticles, Lattice>();
    if (include_air && air_body != nullptr)
    {
        air_body->generateParticles<BaseParticles, Lattice>();
    }

    std::cout << "coil_particles=" << coil_body.getBaseParticles().TotalRealParticles()
              << " plate_particles=" << plate_body.getBaseParticles().TotalRealParticles();
    if (include_air && air_body != nullptr)
    {
        std::cout << " air_particles=" << air_body->getBaseParticles().TotalRealParticles();
    }
    std::cout << std::endl;

    auto &coil_near_surface = coil_body.addBodyPart<NearShapeSurface>();
    auto &plate_near_surface = plate_body.addBodyPart<NearShapeSurface>();
    NearShapeSurface *air_near_surface = nullptr;
    if (include_air && air_body != nullptr)
    {
        air_near_surface = &air_body->addBodyPart<NearShapeSurface>();
    }

    auto &coil_inner = sph_system.addInnerRelation(coil_body);
    auto &plate_inner = sph_system.addInnerRelation(plate_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    StdVec<SPHBody *> all_bodies = {&coil_body, &plate_body};
    if (include_air && air_body != nullptr)
    {
        all_bodies.push_back(air_body);
    }
    host_methods.addStateDynamics<RandomizeParticlePositionCK>(all_bodies).exec();

    auto &update_coil_cell_linked_list = main_methods.addCellLinkedListDynamics(coil_body);
    auto &update_plate_cell_linked_list = main_methods.addCellLinkedListDynamics(plate_body);

    auto &update_coil_inner = main_methods.addRelationDynamics(coil_inner);
    auto &update_plate_inner = main_methods.addRelationDynamics(plate_inner);

    auto &coil_relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(coil_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(coil_body, coil_level_set);
    auto &plate_relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(plate_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(plate_body, plate_level_set);

    auto &coil_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(coil_body);
    auto &plate_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(plate_body);

    auto &coil_update_position = main_methods.addStateDynamics<PositionRelaxationCK>(coil_body);
    auto &plate_update_position = main_methods.addStateDynamics<PositionRelaxationCK>(plate_body);

    auto &coil_levelset_bounding = main_methods.addStateDynamics<LevelsetBounding>(coil_near_surface);
    auto &plate_levelset_bounding = main_methods.addStateDynamics<LevelsetBounding>(plate_near_surface);

    auto &coil_update_smoothing_length_ratio =
        main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(coil_body, coil_level_set);
    auto &plate_update_smoothing_length_ratio =
        main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(plate_body, plate_level_set);

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(coil_body, "SmoothingLengthRatio");
    body_state_recorder.addToWrite<Real>(plate_body, "SmoothingLengthRatio");
    if (include_air && air_body != nullptr)
    {
        body_state_recorder.addToWrite<Real>(*air_body, "SmoothingLengthRatio");
    }

    auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(sph_system);
    write_particle_reload_files.addToReload<Real>(coil_body, "SmoothingLengthRatio");
    write_particle_reload_files.addToReload<Real>(plate_body, "SmoothingLengthRatio");
    registerTeam7CoilPlateNormalsForReload(write_particle_reload_files, coil_body, plate_body);
    if (include_air && air_body != nullptr)
    {
        write_particle_reload_files.addToReload<Real>(*air_body, "SmoothingLengthRatio");
    }

    body_state_recorder.writeToFile(0);

    if (include_air && air_body != nullptr && air_level_set != nullptr && air_near_surface != nullptr)
    {
        auto &air_inner = sph_system.addInnerRelation(*air_body);
        auto &update_air_cell_linked_list = main_methods.addCellLinkedListDynamics(*air_body);
        auto &update_air_inner = main_methods.addRelationDynamics(air_inner);
        auto &air_relaxation_residual =
            main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(air_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(*air_body, *air_level_set);
        auto &air_relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(*air_body);
        auto &air_update_position = main_methods.addStateDynamics<PositionRelaxationCK>(*air_body);
        auto &air_levelset_bounding = main_methods.addStateDynamics<LevelsetBounding>(*air_near_surface);
        auto &air_update_smoothing_length_ratio =
            main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(*air_body, *air_level_set);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            update_coil_cell_linked_list.exec();
            update_plate_cell_linked_list.exec();
            update_air_cell_linked_list.exec();
            update_coil_inner.exec();
            update_plate_inner.exec();
            update_air_inner.exec();

            coil_relaxation_residual.exec();
            plate_relaxation_residual.exec();
            air_relaxation_residual.exec();

            const Real coil_relaxation_step = coil_relaxation_scaling.exec();
            const Real plate_relaxation_step = plate_relaxation_scaling.exec();
            const Real air_relaxation_step = air_relaxation_scaling.exec();
            coil_update_position.exec(coil_relaxation_step);
            plate_update_position.exec(plate_relaxation_step);
            air_update_position.exec(air_relaxation_step);

            coil_levelset_bounding.exec();
            plate_levelset_bounding.exec();
            air_levelset_bounding.exec();
            coil_update_smoothing_length_ratio.exec();
            plate_update_smoothing_length_ratio.exec();
            air_update_smoothing_length_ratio.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
    }
    else
    {
        int ite_p = 0;
        while (ite_p < 1000)
        {
            update_coil_cell_linked_list.exec();
            update_plate_cell_linked_list.exec();
            update_coil_inner.exec();
            update_plate_inner.exec();

            coil_relaxation_residual.exec();
            plate_relaxation_residual.exec();

            const Real coil_relaxation_step = coil_relaxation_scaling.exec();
            const Real plate_relaxation_step = plate_relaxation_scaling.exec();
            coil_update_position.exec(coil_relaxation_step);
            plate_update_position.exec(plate_relaxation_step);

            coil_levelset_bounding.exec();
            plate_levelset_bounding.exec();
            coil_update_smoothing_length_ratio.exec();
            plate_update_smoothing_length_ratio.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(coil_body).exec();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(plate_body).exec();
    std::cout << "NormalFromBodyShapeCK finished (NormalDirection + SignedDistance for coil/plate)." << std::endl;
    write_particle_reload_files.writeToFile();

    std::error_code mkdir_error;
    std::filesystem::create_directories("./reload", mkdir_error);
    const bool metadata_ok = writeTeam7NativeGeometryMetadata("./reload/team7_native_geometry.txt", config);
    const bool case_published = publishTeam7NativeReloadCaseFromWorkingDirectory(config);
    std::cout << "team7_native_geometry_metadata_written=" << (metadata_ok ? 1 : 0)
              << " reload_case_published=" << (case_published ? 1 : 0) << " reload_case_dir="
              << team7NativeReloadCaseDirectory(team7NativeReloadCaseId(config)) << std::endl;
    return metadata_ok && case_published ? 0 : 1;
}
