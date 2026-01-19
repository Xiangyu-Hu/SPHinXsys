/**
 * @file 	mr_freestream_flow_around_cylinder.cpp
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Define basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 200.0;                /**< TOTAL SIMULATION TIME*/
Real start_up_time = total_physical_time / 10.0; /**< explicity startup time*/
Real DL = 30.0;                                  /**< Domain length. */
Real DH = 16.0;                                  /**< Domain height. */
Real particle_spacing_ref = 0.4;                 /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0;    /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;            /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);      /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;                 /**< Radius of the cylinder. */
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
AlignedBox emitter_box(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize);

Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, -0.25 * DH) + disposer_halfsize;
AlignedBox disposer_box(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize);
//----------------------------------------------------------------------
//	Define adaptation
//----------------------------------------------------------------------
AdaptiveWithinShape water_body_adaptation(particle_spacing_ref, 1.3, 1.0, 1);

GeometricShapeBox refinement_region(
    BoundingBoxd(Vecd(-DL_sponge - BW, 0.5 * DH - 0.1 * DL), Vecd(DL + BW, 0.5 * DH + 0.1 * DL)),
    "RefinementRegion");
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    FreeStreamVelocity() : u_ref_(U_f), t_ref_(start_up_time) {};
    Real getAxisVelocity(const Vecd &input_position, const Real &input_axis_velocity, Real time)
    {
        return time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * time / t_ref_)) : u_ref_;
    };
};
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
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);
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

    StdVec<RealBody *> real_bodies = {&water_body, &cylinder};
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        water_body.generateParticles<BaseParticles, Lattice>(*refinement_region_level_set_shape);
        cylinder.generateParticles<BaseParticles, Lattice>();

        auto &near_water_body_surface = water_body.addBodyPart<NearShapeSurface>(*outer_boundary_level_set_shape);
        auto &near_cylinder_surface = cylinder.addBodyPart<NearShapeSurface>();
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
        relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(water_body_inner)
                                     .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_body_contact)
                                     .addPostStateDynamics<LevelsetKernelGradientIntegral>(water_body, *outer_boundary_level_set_shape));
        relaxation_residual.add(&main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(cylinder_inner)
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
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&water_body, &cylinder});
        //----------------------------------------------------------------------
        //	Prepare for the time integration loop.
        //----------------------------------------------------------------------
        level_set_bounding.exec();
        update_smoothing_length_ratio.exec();
        //----------------------------------------------------------------------
        //	First output before the simulation.
        body_state_recorder.writeToFile();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 2000)
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

        host_methods.addStateDynamics<NormalFromBodyShapeCK>(cylinder).exec();
        write_particle_reload_files.addToReload<Vecd>(cylinder, "NormalDirection");
        write_particle_reload_files.addToReload<Real>(water_body, "SmoothingLengthRatio");
        write_particle_reload_files.writeToFile();

        std::cout << "The physics relaxation process finish !" << std::endl;
        return 0;
    }
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_body.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_body.getName())
        ->reloadExtraVariable<Real>("SmoothingLengthRatio");
    // //----------------------------------------------------------------------
    // //	Creating body parts.
    // //----------------------------------------------------------------------
    auto &emitter = water_body.addBodyPart<AlignedBoxByParticle>(emitter_box);
    auto &disposer = water_body.addBodyPart<AlignedBoxByCell>(disposer_box);

    cylinder.defineMaterial<Solid>();
    cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    auto &water_body_inner = sph_system.addInnerRelation(water_body);
    auto &water_body_contact = sph_system.addContactRelation(water_body, cylinder);
    auto &fluid_observer_contact = sph_system.addContactRelation(fluid_observer, water_body);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
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
    auto &update_cylinder_cell_linked_list = main_methods.addCellLinkedListDynamics(cylinder);
    ParticleDynamicsGroup update_water_body_configuration;
    update_water_body_configuration.add(&main_methods.addCellLinkedListDynamics(water_body));
    update_water_body_configuration.add(&main_methods.addRelationDynamics(water_body_inner, water_body_contact));
    auto &update_observer_relation = main_methods.addRelationDynamics(fluid_observer_contact);
    auto &particle_sort = main_methods.addSortDynamics(water_body);

    auto &time_dependent_gravity = main_methods.addStateDynamics<GravityForceCK<StartupAcceleration>>(water_body, Vec2d(U_f, 0.0), start_up_time);
    auto &water_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_body);
    auto &water_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_body);

    auto &fluid_boundary_indicator =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::FreeSurfaceIndicationCK>(water_body_inner)
            .addPostContactInteraction(water_body_contact);

    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(water_body_inner)
            .addPostContactInteraction(water_body_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeStream>(water_body);

    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_body_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_body_contact)
            .addPostStateDynamics<fluid_dynamics::FreeStreamCondition<FreeStreamVelocity>>(water_body);
    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_body_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_body_contact);

    auto &transport_correction =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(water_body_inner)
            .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_body_contact)
            .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter, BulkParticles>(water_body)
            .addPostStateDynamics<ConstantConstraintCK, Vecd>(emitter, "Displacement", Vecd::Zero());

    auto &fluid_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_body, U_f);
    auto &fluid_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_body);

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<
                        fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(water_body_inner)
            .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(water_body_contact);

    auto &emitter_injection = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK>(emitter, inlet_particle_buffer);
    auto &inflow_condition = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowConditionCK, FreeStreamVelocity>(emitter);
    auto &disposer_indication = main_methods.addStateDynamics<fluid_dynamics::WithinDisposerIndication>(disposer);
    auto &particle_deletion = main_methods.addStateDynamics<fluid_dynamics::OutflowParticleDeletion>(water_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &write_real_body_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    write_real_body_states.addToWrite<Real>(water_body, "Density");
    write_real_body_states.addToWrite<Real>(water_body, "SmoothingLengthRatio");
    write_real_body_states.addToWrite<int>(water_body, "Indicator");
    write_real_body_states.addToWrite<Vecd>(cylinder, "NormalDirection");
    auto &write_fluid_observation = main_methods.addObserveRecorder<Vecd>("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper time_stepper(sph_system, total_physical_time);
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t advection_steps = 0;
    int screening_interval = 100;
    int observation_interval = screening_interval / 2;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 200.0);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    update_cylinder_cell_linked_list.exec();
    update_water_body_configuration.exec();

    time_dependent_gravity.exec();
    fluid_boundary_indicator.exec();
    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    transport_correction.exec();
    fluid_viscous_force.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    update_observer_relation.exec();
    write_fluid_observation.writeToFile();
    /** statistics for computing time. */
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        inflow_condition.exec();
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            water_update_particle_position.exec();

            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Physical Time = " << time_stepper.getPhysicalTime()
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                update_observer_relation.exec();
                write_fluid_observation.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                write_real_body_states.writeToFile();
            }

            /** Particle creation, deletion, sort and update configuration. */
            time_instance = TickCount::now();
            emitter_injection.exec();
            disposer_indication.exec();
            particle_deletion.exec();

            if (advection_steps % 100)
            {
                particle_sort.exec();
            }

            update_water_body_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            time_dependent_gravity.exec();
            fluid_boundary_indicator.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }

    //----------------------------------------------------------------------
    // Summary for wall time used for real computations.
    //----------------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    return 0;
}
