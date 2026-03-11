/**
 * @file 	particle_generation_em.cpp
 * @brief 	Particle generation and relaxation for a-phi EM case: coil, plate (and optionally air).
 * @details Coil and plate are single-resolution; uses par_ck level set and CK relaxation for SYCL GPU.
 *          Relations must be created via sph_system.addInnerRelation(body) for CK API compatibility.
 * @author 	Generated for SPHinXsys electromagnetic branch
 */

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	File paths to STL data (coil.stl and plate.stl in case input/ folder)
//----------------------------------------------------------------------
const std::string path_coil_stl = "./input/coil.stl";
const std::string path_plate_stl = "./input/plate.stl";

//----------------------------------------------------------------------
//	Geometry parameters (from asymmetric conductor with hole)
//----------------------------------------------------------------------
const Real plate_xy = 294.0;
const Real plate_z = 19.0;
const Real air_margin = 100.0;

const Vec3d air_box_lower(-1353.0, -1353.0, -300.0);
const Vec3d air_box_upper(1647.0, 1647.0, 449.0);

const Real dp_0 = 6.0;
/** Single-resolution air for testing (lower resolution to avoid hang) */
const Real dp_air_single = 10.0;
/** Air multi-resolution (commented out when testing single-res air): finest = dp_air_finest, coarsest = dp_air_finest * 2^air_refinement_levels */
const Real dp_air_finest = 3.0;       /* 最小分辨率（最密）粒子间距 */
const int air_refinement_levels = 4;  /* 向外扩 3 层（最细 1 层 + 再粗 3 层），最粗 = dp_air_finest * 2^3 = 24 */
const Real dp_air_coarsest = dp_air_finest * pow(2.0, air_refinement_levels);  /* 24.0 */
BoundingBoxd system_domain_bounds(air_box_lower, air_box_upper);

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
//	Plate shape: single STL
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
    explicit AirShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize = 0.5 * (air_box_upper - air_box_lower);
        Vecd center = 0.5 * (air_box_lower + air_box_upper);
        Transform translation(center);
        add<GeometricShapeBox>(translation, halfsize, "OuterBoundary");
        subtract<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
        subtract<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
    }
};

//----------------------------------------------------------------------
//	Inner boundary shape: union of coil and plate (for refinement-distance only).
//----------------------------------------------------------------------
class InnerBoundaryShape : public ComplexShape
{
  public:
    explicit InnerBoundaryShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_coil_stl, Vec3d::Zero(), 1.0, "Coil");
        add<TriangleMeshShapeSTL>(path_plate_stl, Vec3d::Zero(), 1.0, "Plate");
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

//----------------------------------------------------------------------
//	Main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    //----------------------------------------------------------------------
    //	Coil (commented out for testing air only)
    //----------------------------------------------------------------------
    // SolidBody coil_body(sph_system, makeShared<CoilShape>("Coil"));
    // coil_body.defineAdaptation<SPHAdaptation>(1.15, 2.0);
    // coil_body.defineMaterial<Solid>();
    // LevelSetShape *coil_level_set = coil_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet();
    // coil_body.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	Plate (commented out for testing air only)
    //----------------------------------------------------------------------
    // SolidBody plate_body(sph_system, makeShared<PlateShape>("Plate"));
    // plate_body.defineAdaptation<SPHAdaptation>(1.15, 2.0);
    // plate_body.defineMaterial<Solid>();
    // LevelSetShape *plate_level_set = plate_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet();
    // plate_body.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	Air: single-resolution, lower resolution for testing (kept but commented out)
    //	spacing = dp_0 / refinement_to_global => use refinement_to_global = dp_0/dp_air_single so spacing = dp_air_single
    //----------------------------------------------------------------------
    // auto &air_body = sph_system.addBody<RealBody>(makeShared<AirShape>("Air"));
    // air_body.defineAdaptation<SPHAdaptation>(1.15, dp_0 / dp_air_single);
    // air_body.defineMaterial<Solid>();
    // LevelSetShape *air_level_set =
    //     air_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet();
    // air_body.generateParticles<BaseParticles, Lattice>();
    // auto &near_air_surface = air_body.addBodyPart<NearShapeSurface>();

    //----------------------------------------------------------------------
    //	Air: multi-resolution (active) – only air, refinement near coil/plate inner boundary
    //----------------------------------------------------------------------
    auto &inner_boundary_shape = sph_system.addShape<InnerBoundaryShape>("InnerBoundary");
    AdaptiveNearInnerSurface air_adaptation(
        dp_air_coarsest, 1.15, 1.0, air_refinement_levels, &inner_boundary_shape);
    auto &air_body = sph_system.addAdaptiveBody<FluidBody, AdaptiveNearInnerSurface>(
        air_adaptation, makeShared<AirShape>("Air"));
    air_body.defineMaterial<Solid>();
    LevelSetShape *air_level_set =
        air_body.defineBodyLevelSetShape(par_ck)->correctLevelSetSign()->cleanLevelSet();
    air_body.generateParticles<BaseParticles, Lattice>();
    auto &near_air_surface = air_body.addBodyPart<NearShapeSurface>();
    // Level set used for smoothing-length update: only inner boundary (coil + plate), not outer air box
    LevelSetShape *inner_refinement_level_set =
        sph_system.addShape<LevelSetShape>(air_body, inner_boundary_shape).writeLevelSet();

    //----------------------------------------------------------------------
    //	Relations: air only (coil/plate commented out)
    //  Use system-managed inner relation (required for CK multi-resolution)
    //----------------------------------------------------------------------
    auto &air_inner = sph_system.addInnerRelation(air_body);
    // auto &coil_inner = sph_system.addInnerRelation(coil_body);
    // auto &plate_inner = sph_system.addInnerRelation(plate_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    host_methods.addStateDynamics<RandomizeParticlePositionCK>(air_body).exec();
    // host_methods.addStateDynamics<RandomizeParticlePositionCK>(coil_body).exec();
    // host_methods.addStateDynamics<RandomizeParticlePositionCK>(plate_body).exec();

    auto &update_cell_linked_list = main_methods.addCellLinkedListDynamics(air_body);
    auto &update_inner_relation = main_methods.addRelationDynamics(air_inner);
    ParticleDynamicsGroup update_configuration;
    update_configuration.add(&update_cell_linked_list);
    update_configuration.add(&update_inner_relation);

    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(air_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(air_body, *air_level_set);
    // ParticleDynamicsGroup relaxation_residual;
    // relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(coil_inner)
    //                             .addPostStateDynamics<LevelsetKernelGradientIntegral>(coil_body, *coil_level_set));
    // relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(plate_inner)
    //                             .addPostStateDynamics<LevelsetKernelGradientIntegral>(plate_body, *plate_level_set));
    // relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(air_inner)
    //                             .addPostStateDynamics<LevelsetKernelGradientIntegral>(air_body, *air_level_set));

    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(air_body);
    // ReduceDynamicsGroup<ReduceMin> relaxation_scaling =
    //     main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(all_bodies);

    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(air_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_air_surface);
    auto &update_smoothing_length_ratio =
        main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(air_body, *inner_refinement_level_set);

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(air_body, "SmoothingLengthRatio");  /* only for adaptive air */

    body_state_recorder.writeToFile(0);
    const int relaxation_steps = 1000;
    int ite_p = 0;
    while (ite_p < relaxation_steps)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();
        update_smoothing_length_ratio.exec();  /* only for multi-resolution air */

        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
            body_state_recorder.writeToFile(ite_p);
        }
    }
    std::cout << "Particle generation and relaxation finished." << std::endl;
    body_state_recorder.writeToFile(ite_p);
    return 0;
}
