/**
 * @file 	dam_breach_cks.cpp
 * @brief 	Dam breach simulation.
 * @author  Shuang Li, Xiangyu Hu
 */
#include "dam_breach_ck.h" // case file to setup the test case
#include "sphinxsys_ck.h"
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref, 16);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    SolidBody soil_block(sph_system, makeShared<SoilBlock>("SoilBlock"));
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle, cohesion);
    soil_block.generateParticles<BaseParticles, Lattice>();

    SolidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(350.0);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    AlignedBoxByParticle emitter(water_block, AlignedBox(yAxis, Transform(emitter_translation), emitter_halfsize));
    AlignedBoxByCell outflow_disposer(water_block, AlignedBox(xAxis, Transform(disposer_translation), disposer_halfsize));
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Inner<> soil_block_inner(soil_block);
    Contact<> water_wall_contact(water_block, {&wall_boundary});
    Contact<> water_soil_contact(water_block, {&soil_block});
    Contact<> soil_wall_contact(soil_block, {&wall_boundary});
    Contact<> soil_water_contact(soil_block, {&water_block});
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par);
    auto &host_methods = sph_solver.addParticleMethodContainer(par);
    using MainExecutionPolicy = execution::ParallelPolicy;
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &soil_cell_linked_list = main_methods.addCellLinkedListDynamics(soil_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &water_block_update_complex_relation = main_methods.addRelationDynamics(water_block_inner, water_wall_contact, water_soil_contact);
    auto &soil_block_update_complex_relation = main_methods.addRelationDynamics(soil_block_inner, soil_wall_contact, soil_water_contact);
    auto &particle_sort = main_methods.addSortDynamics(water_block);
    auto &soil_particle_sort = main_methods.addSortDynamics(soil_block);

    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(water_block, gravity);
    auto &soil_constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(soil_block, gravity);
    auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary); // run on CPU

    auto &soil_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(soil_block);
    auto &soil_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(soil_block);
    auto &soil_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(soil_block); // run on CPU
    auto &soil_surface_indicator =  
        main_methods.addInteractionDynamics<
                        fluid_dynamics::FreeSurfaceIndicationCK, WithUpdate>(soil_block_inner)
            .addContactInteraction(soil_wall_contact);
    auto &stress_diffusion = main_methods.addInteractionDynamics<continuum_dynamics::StressDiffusionInnerCK>(soil_block_inner);
    auto &soil_acoustic_step_1st_half =  
        main_methods.addInteractionDynamics<
                        continuum_dynamics::PlasticAcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>(soil_block_inner)
            .addContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrection>(soil_wall_contact);
    auto &soil_acoustic_step_2nd_half = 
        main_methods.addInteractionDynamics<
                        continuum_dynamics::PlasticAcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>(soil_block_inner)
            .addContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrection>(soil_wall_contact)
            .addContactInteraction<Fluid, AcousticRiemannSolverCK, NoKernelCorrection>(soil_water_contact);
    auto &soil_density_regularization = 
        main_methods.addInteractionDynamics<
                        fluid_dynamics::DensityRegularization, WithUpdate, FreeSurface, AllParticles>(soil_block_inner)
            .addContactInteraction(soil_wall_contact);

    auto &water_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);
    auto &fluid_linear_correction_matrix =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(water_block_inner, 0.5)
            .addContactInteraction(water_wall_contact);
    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact)
            .addContactInteraction<Soil, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_soil_contact);
    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_wall_contact)
            .addContactInteraction<AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_soil_contact);
    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::DensityRegularization, WithUpdate, FreeSurface, AllParticles>(water_block_inner)
            .addContactInteraction(water_wall_contact)
            .addContactInteraction(water_soil_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK>
        transport_correction_ck(water_block_inner, water_wall_contact);

    auto &fluid_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_ref);
    auto &fluid_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_state_recorder.addToWrite<Real>(water_block, "Density");
    body_state_recorder.addToWrite<Vecd>(soil_block, "ShearVelocity");
    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);
    /*Inflow condition*/
    auto &inflow_condition = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowConditionCK<AlignedBoxByParticle, InletInflowCondition>>(emitter);
    auto &emitter_injection = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK<AlignedBoxByParticle>>(emitter, inlet_buffer);
    /*Outflow condition*/
    fluid_dynamics::BidirectionalBoundaryCK<MainExecutionPolicy, LinearCorrectionCK, PressurePrescribed<>> outflow_condition(outflow_disposer, inlet_buffer,0.0);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(20.0);
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        time_stepper.setPhysicalTime(restart_io.readRestartFiles(sph_system.RestartStep()));
    }
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep() + 1;
    int screening_interval = 100;
    int restart_output_interval = screening_interval * 10;
    auto &state_recording = time_stepper.addTriggerByInterval(0.1);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec(); // run particle dynamics with host kernels first
    soil_normal_direction.exec(); // run particle dynamics with host kernels first
    constant_gravity.exec();
    soil_constant_gravity.exec();

    water_cell_linked_list.exec();
    soil_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    soil_block_update_complex_relation.exec();

    fluid_density_regularization.exec();
    soil_density_regularization.exec();
    water_advection_step_setup.exec();
    soil_advection_step_setup.exec();
    fluid_linear_correction_matrix.exec();
    outflow_condition.tagBufferParticles();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        inflow_condition.exec();
        outflow_condition.applyBoundaryCondition(acoustic_dt);
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        /*Soil Dynamics*/
        stress_diffusion.exec();
        soil_acoustic_step_1st_half.exec(acoustic_dt);
        soil_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            soil_normal_direction.exec();

            water_update_particle_position.exec();
            soil_update_particle_position.exec();


            emitter_injection.exec();
            outflow_condition.deleteParticles();
            /** Output body state during the simulation according output_interval. */
            time_instance = TickCount::now();
            /** screen output, write body observables and restart files  */
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }


            if (advection_steps % restart_output_interval == 0)
            {
                restart_io.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (advection_steps % 100)
            {
                particle_sort.exec();
                soil_particle_sort.exec();
            }
            water_cell_linked_list.exec();
            soil_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            soil_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            soil_density_regularization.exec();
            water_advection_step_setup.exec();
            soil_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            outflow_condition.tagBufferParticles();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }
    return 0;
};
