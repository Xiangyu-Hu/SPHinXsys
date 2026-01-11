/**
 * @file 	mr_freestream_flow_around_cylinder.cpp
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Define basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;                               /**< Domain length. */
Real DH = 16.0;                               /**< Domain height. */
Real particle_spacing_ref = 0.4;              /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;         /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);   /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;              /**< Radius of the cylinder. */
// Observation locations
Vec2d point_1(3.0, 5.0);
Vec2d point_2(4.0, 5.0);
Vec2d point_3(5.0, 5.0);
StdVec<Vecd> observation_locations = {point_1, point_2, point_3};
//----------------------------------------------------------------------
//	Define global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
Real Re = 100.0;                                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
GeometricShapeBox outer_boundary(BoundingBoxd(Vecd(-DL_sponge, 0.0), Vecd(DL, DH)), "OuterBoundary");
GeometricShapeBall cylinder_shape(insert_circle_center, insert_circle_radius, "Cylinder");

Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
GeometricShapeBox emitter_shape(Transform(emitter_translation), emitter_halfsize);

Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
GeometricShapeBox emitter_buffer_shape(Transform(emitter_buffer_translation), emitter_buffer_halfsize);

Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
GeometricShapeBox disposer_shape(Transform(disposer_translation), disposer_halfsize);
//----------------------------------------------------------------------
//	Define adaptation
//----------------------------------------------------------------------
AdaptiveWithinShape water_body_adaptation(particle_spacing_ref, 1.3, 1.0, 1);

GeometricShapeBox refinement_region(
    BoundingBoxd(Vecd(-DL_sponge - BW, 0.5 * DH - 0.1 * DL), Vecd(DL + BW, 0.5 * DH + 0.1 * DL)),
    "RefinementRegion");
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(true);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    auto &water_body_shape = sph_system.addShape<ComplexShape>("WaterBody");
    water_body_shape.add(&outer_boundary);
    water_body_shape.subtract(&cylinder_shape);
    auto &water_body = sph_system.addAdaptiveBody<FluidBody>(water_body_adaptation, water_body_shape);
    LevelSetShape *outer_boundary_level_set_shape =
        water_body.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet();
    LevelSetShape *refinement_region_level_set_shape =
        sph_system.addShape<LevelSetShape>(water_body, refinement_region).writeLevelSet();

    auto &cylinder = sph_system.addBody<SolidBody>(cylinder_shape);
    cylinder.defineAdaptationRatios(1.15, 4.0);
    LevelSetShape *cylinder_level_set_shape = cylinder.defineBodyLevelSetShape()->writeLevelSet();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        water_body.generateParticles<BaseParticles, Lattice>(*refinement_region_level_set_shape);
        auto &near_water_body_surface = water_body.addBodyPart<NearShapeSurface>(*outer_boundary_level_set_shape);
        
        cylinder.generateParticles<BaseParticles, Lattice>();
        auto &near_cylinder_surface = cylinder.addBodyPart<NearShapeSurface>();

        StdVec<RealBody *> real_bodies = {&water_body, &cylinder};
        StdVec<NearShapeSurface *> near_body_surfaces = {&near_water_body_surface, &near_cylinder_surface};
        //----------------------------------------------------------------------
        //	Define body relation map.
        //	The contact map gives the topological connections between the bodies.
        //	Basically the the range of bodies to build neighbor particle lists.
        //  Generally, we first define all the inner relations, then the contact relations.
        //  At last, we define the complex relaxations by combining previous defined
        //  inner and contact relations.
        //----------------------------------------------------------------------
        auto &water_body_inner = sph_system.addInnerRelation(water_body);
        auto &cylinder_inner = sph_system.addInnerRelation(cylinder);
        auto &water_body_contact = sph_system.addContactRelation(water_body, cylinder);
        //----------------------------------------------------------------------
        // Define SPH solver with particle methods and execution policies.
        // Generally, the host methods should be able to run immediately.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(sph_system);
        //----------------------------------------------------------------------
        // Define the numerical methods used in the simulation.
        // Note that there may be data dependence on the sequence of constructions.
        // Generally, the configuration dynamics, such as update cell linked list,
        // update body relations, are defined first.
        // Then the geometric models or simple objects without data dependencies,
        // such as gravity, initialized normal direction.
        // After that, the major physical particle dynamics model should be introduced.
        // Finally, the auxiliary models such as time step estimator, initial condition,
        // boundary condition and other constraints should be defined.
        //----------------------------------------------------------------------
        auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
        host_methods.addStateDynamics<RandomizeParticlePositionCK>(real_bodies).exec(); // host method able to run immediately
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
        ParticleDynamicsGroup update_cell_linked_list = main_methods.addCellLinkedListDynamics(real_bodies);
        ParticleDynamicsGroup update_relation;
        update_relation.add(&main_methods.addRelationDynamics(water_body_inner, water_body_contact));
        update_relation.add(&main_methods.addRelationDynamics(cylinder_inner));
        ParticleDynamicsGroup update_configuration = update_cell_linked_list + update_relation;

        ParticleDynamicsGroup relaxation_residual;
        relaxation_residual.add(&main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(water_body_inner)
                                     .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_body_contact)
                                     .addPostStateDynamics<LevelsetKernelGradientIntegral>(water_body, *outer_boundary_level_set_shape));
        relaxation_residual.add(&main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(cylinder_inner)
                                     .addPostStateDynamics<LevelsetKernelGradientIntegral>(cylinder, *cylinder_level_set_shape));

        ReduceDynamicsGroup relaxation_scaling = main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(real_bodies);
        ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(real_bodies);
        ParticleDynamicsGroup level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surfaces);
        auto &update_smoothing_length_ratio = main_methods.addStateDynamics<UpdateSmoothingLengthRatio>(water_body, *refinement_region_level_set_shape);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
        body_state_recorder.addToWrite<Real>(water_body, "SmoothingLengthRatio");
        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        body_state_recorder.writeToFile();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            update_configuration.exec();
            relaxation_residual.exec();
            Real relaxation_step = relaxation_scaling.exec();
            update_particle_position.exec(relaxation_step);
            level_set_bounding.exec();
            update_smoothing_length_ratio.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;
        return 0;
    }
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_body.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_body.getName());

    cylinder.defineMaterial<Solid>();
    cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName());

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_locations);

    return 0;
}
