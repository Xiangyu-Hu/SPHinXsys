/**
 * @file  particle_generation_em.cpp
 * @brief particle geneartion and relaxation for TEAM7 case.
 *
 * NOTE:

 */

#include "sphinxsys.h"


using namespace SPH;

namespace
{
//----------------------------------------------------------------------
//  STL paths
//----------------------------------------------------------------------
const std::string path_coil_stl = "./input/coil.stl";
const std::string path_plate_stl = "./input/plate.stl";

//----------------------------------------------------------------------
//	Geometry parameters (from asymmetric conductor with hole)
//----------------------------------------------------------------------
const Real plate_xy = 294.0;
const Real plate_z = 19.0;
const Real air_margin = 100.0;

/*big air box*/
const Vec3d air_box_lower(-1353.0, -1353.0, -300.0);
const Vec3d air_box_upper(1647.0, 1647.0, 449.0);

/*small air box*/
const Vec3d air_box_lower_small(-50.0, -50.0, -50.0);
const Vec3d air_box_upper_small(350.0, 350.0, 200.0);

const Real dp_0 = 6.0;
/** Single-resolution air for testing (lower resolution to avoid hang) */
const Real dp_air_single = 10.0;
/** Air multi-resolution (commented out when testing single-res air): finest = dp_air_finest, coarsest = dp_air_finest * 2^air_refinement_levels */
const Real dp_air_finest = 3.0;       
const int air_refinement_levels = 4; 
const Real dp_air_coarsest = dp_air_finest * pow(2.0, air_refinement_levels);  /* 24.0 */
BoundingBoxd system_domain_bounds(air_box_lower, air_box_upper);
BoundingBoxd system_domain_bounds_small(air_box_lower_small, air_box_upper_small);
//----------------------------------------------------------------------
//	Coil shape: single STL
//----------------------------------------------------------------------
class CoilShape : public ComplexShape
{
  public:
    explicit CoilShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0);
    }
};

//----------------------------------------------------------------------
//	Coil shape: single STL
//----------------------------------------------------------------------
class PlateShape : public ComplexShape
{
  public:
    explicit PlateShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0);
    }
};

//----------------------------------------------------------------------
//	Air shape: box (OuterBoundary) minus coil minus plate; parts named for component level set.
//----------------------------------------------------------------------
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

//----------------------------------------------------------------------
//	Adaptation that refines only near the inner boundary (coil + plate), not the outer box.
//----------------------------------------------------------------------
class AdaptiveNearInnerSurface : public AdaptiveNearSurface
{
    Shape *inner_shape_;

  public:
    AdaptiveNearInnerSurface(Real global_resolution, Real h_spacing_ratio, Real refinement_to_global,
                            int local_refinement_level, Shape *inner_shape)
        : AdaptiveNearSurface(global_resolution, h_spacing_ratio, refinement_to_global, local_refinement_level),
          inner_shape_(inner_shape)
    {
    }
    virtual Real getLocalSpacing(Shape &shape, const Vecd &position) override
    {
        Real phi = fabs(inner_shape_->findSignedDistance(position));
        return smoothedSpacing(phi, spacing_ref_);
    }
};

} // namespace

//----------------------------------------------------------------------
//	Main program (must be in global namespace, not anonymous — otherwise gtest_main wins at link time)
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    

    SPHSystem sph_system(system_domain_bounds_small, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);


    //------------------------------------------------------------------
    //  Particle generation: coil + plate + adaptive air
    //------------------------------------------------------------------
    SolidBody coil_body(sph_system, makeShared<CoilShape>("Coil"));
    coil_body.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    coil_body.defineMaterial<Solid>();
    LevelSetShape *coil_level_set =
        coil_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet()->writeLevelSet();
    coil_body.generateParticles<BaseParticles, Lattice>();

    SolidBody plate_body(sph_system, makeShared<PlateShape>("Plate"));
    plate_body.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    plate_body.defineMaterial<Solid>();
    LevelSetShape *plate_level_set =
        plate_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet()->writeLevelSet();
    plate_body.generateParticles<BaseParticles, Lattice>();


    //----------------------------------------------------------------------
    //	Air: single-resolution, lower resolution for testing (kept but commented out)
    //	spacing = dp_0 / refinement_to_global => use refinement_to_global = dp_0/dp_air_single so spacing = dp_air_single
    //----------------------------------------------------------------------
    RealBody air_body(sph_system, makeShared<AirShape>("Air", air_box_lower_small, air_box_upper_small));
    LevelSetShape *air_level_set =
        air_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet()->writeLevelSet();
    air_body.generateParticles<BaseParticles, Lattice>();

     //----------------------------------------------------------------------
    //	Air: multi-resolution (active) – only air, refinement near coil/plate inner boundary
    //----------------------------------------------------------------------
    // auto &inner_boundary_shape = sph_system.addShape<InnerBoundaryShape>("InnerBoundary");
    // AdaptiveNearInnerSurface air_adaptation(
    //     dp_air_coarsest, 1.15, 1.0, air_refinement_levels, &inner_boundary_shape);
    // auto &air_body = sph_system.addAdaptiveBody<FluidBody, AdaptiveNearInnerSurface>(
    //     air_adaptation, makeShared<AirShape>("Air", air_box_lower, air_box_upper));
    // LevelSetShape *air_level_set =
    //     air_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet();
    // air_body.generateParticles<BaseParticles, Lattice>();
    // air_body.defineMaterial<Fluid>();
    // auto &near_air_surface = air_body.addBodyPart<NearShapeSurface>();
    // // Level set used for smoothing-length update: only inner boundary (coil + plate), not outer air box
    // LevelSetShape *inner_refinement_level_set =
    //     sph_system.addShape<LevelSetShape>(air_body, inner_boundary_shape).writeLevelSet();

    
    //----------------------------------------------------------------------
    //	Relations: coil, plate, air
    //----------------------------------------------------------------------
    NearShapeSurface coil_near_body_surface(coil_body);
    NearShapeSurface plate_near_body_surface(plate_body);
    NearShapeSurface air_near_body_surface(air_body);

    Inner<> coil_inner(coil_body);
    Inner<> plate_inner(plate_body);
    Inner<> air_inner(air_body);


    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    StdVec<RealBody *> all_bodies = {&coil_body, &plate_body, &air_body};

    host_methods.addStateDynamics<RandomizeParticlePositionCK>(all_bodies).exec();

    ParticleDynamicsGroup update_cell_linked_list = main_methods.addCellLinkedListDynamics(all_bodies);
    ParticleDynamicsGroup update_relation;
    update_relation.add(&main_methods.addRelationDynamics(coil_inner));
    update_relation.add(&main_methods.addRelationDynamics(plate_inner));
    update_relation.add(&main_methods.addRelationDynamics(air_inner));

    ParticleDynamicsGroup update_configuration = update_cell_linked_list + update_relation;
    ParticleDynamicsGroup relaxation_residual;
    relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(coil_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(coil_body, *coil_level_set));
    relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(plate_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(plate_body, *plate_level_set));
    relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(air_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(air_body, *air_level_set));

    ReduceDynamicsGroup relaxation_scaling = main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(all_bodies);

    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(all_bodies);
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(coil_near_body_surface));
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(plate_near_body_surface));
    update_particle_position.add(&main_methods.addStateDynamics<LevelsetBounding>(air_near_body_surface));

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(sph_system);

       
    body_state_recorder.writeToFile(0);

    int ite_p = 0;
    while (ite_p < 1000)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);

        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            body_state_recorder.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;
    write_particle_reload_files.writeToFile();
    return 0;
}
