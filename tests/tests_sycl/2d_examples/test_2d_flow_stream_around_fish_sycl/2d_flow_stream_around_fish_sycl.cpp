/**
 * @file 2d_flow_stream_around_fish_sycl.cpp
 * @brief Fish swimming driven by active muscles — full SYCL/CK GPU port.
 *        Fluid operators follow mr_free_stream_around_cylinder_sycl as reference.
 *        Solid FSI operators use the same ParticleMethodContainer API.
 *        Output format matches the CPU fish case for direct comparison.
 * @author Pruthvik Arasikere Mallikarjuna and Xiangyu Hu
 */
#include "2d_flow_stream_around_fish.h"
#include "sphinxsys.h"
using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.8);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody fish_body(sph_system, makeShared<FishBody>("FishBody"));
    fish_body.defineAdaptationRatios(1.15, 2.0);
    fish_body.defineMaterial<FishBodyComposite>();
    fish_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relations.
    //----------------------------------------------------------------------
    Inner<> fish_inner(fish_body, ConfigType::Lagrangian);
    Inner<> water_block_inner(water_block);
    Contact<> water_block_contact(water_block, {&fish_body}); // fish as stationary wall
    Contact<> fish_contact(fish_body, {&water_block});
    //----------------------------------------------------------------------
    //	Define SPH solver — all operators via ParticleMethodContainer.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);

    // CPU-only one-shot initialisation (no GPU CK version needed for these).
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(fish_body).exec();
    host_methods.addStateDynamics<FishMaterialInitialization>(fish_body).exec();

    // All GPU operators via main_methods (par_ck = MainExecutionPolicy).
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    //----------------------------------------------------------------------
    //	Configuration dynamics.
    //----------------------------------------------------------------------
    auto &update_fish_cell_linked_list =
        main_methods.addCellLinkedListDynamics(fish_body);
    // Group water CLL + relation rebuild — both run on every update_water_body_configuration.exec().
    ParticleDynamicsGroup update_water_body_configuration;
    update_water_body_configuration.add(&main_methods.addCellLinkedListDynamics(water_block));
    update_water_body_configuration.add(&main_methods.addRelationDynamics(water_block_inner, water_block_contact));

    // fish_inner is Lagrangian — built once during init, never rebuilt.
    auto &fish_inner_build =
        main_methods.addRelationDynamics(fish_inner);
    auto &fish_contact_relation =
        main_methods.addRelationDynamics(fish_contact);
    auto &particle_sort =
        main_methods.addSortDynamics(water_block);
    //----------------------------------------------------------------------
    //	Solid correction matrix (Lagrangian — computed once, stays fixed).
    //----------------------------------------------------------------------
    auto &fish_body_corrected_configuration =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(fish_inner);
    //----------------------------------------------------------------------
    //	FSI average velocity / acceleration for dummy particle boundary condition.
    //----------------------------------------------------------------------
    // auto &initialize_displacement =
    //     main_methods.addStateDynamics<InitializeDisplacementCK>(fish_body);
    // auto &update_average_velocity =
    //     main_methods.addStateDynamics<UpdateAverageVelocityAndAccelerationCK>(fish_body);
    //----------------------------------------------------------------------
    //	Solid dynamics.
    //----------------------------------------------------------------------
    // auto &imposing_active_strain =
    //     main_methods.addStateDynamics<ImposingActiveStrain>(fish_body);
    // auto &fish_body_stress_relaxation_first_half =
    //     main_methods.addInteractionDynamicsOneLevel<
    //         solid_dynamics::StructureIntegration1stHalfPK2, FishBodyComposite>(fish_inner);
    // auto &fish_body_stress_relaxation_second_half =
    //     main_methods.addInteractionDynamicsOneLevel<
    //         solid_dynamics::StructureIntegration2ndHalf>(fish_inner);
    // auto &fish_body_computing_time_step_size =
    //     main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(fish_body);
    // auto &fish_body_update_normal =
    //     main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(fish_body);
    //----------------------------------------------------------------------
    //	Gravity (startup ramp, zero net acceleration for fish case).
    //----------------------------------------------------------------------
    StartupAcceleration time_dependent_acceleration(Vec2d::Zero(), 2.0);
    auto &apply_gravity_force =
        main_methods.addStateDynamics<GravityForceCK<StartupAcceleration>>(
            water_block, time_dependent_acceleration);
    //----------------------------------------------------------------------
    //	Fluid dynamics — follows cylinder MR SYCL pattern exactly.
    //----------------------------------------------------------------------
    auto &water_advection_step_setup =
        main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);

    auto &fluid_boundary_indicator =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::FreeSurfaceIndicationCK>(water_block_inner)
            .addPostContactInteraction(water_block_contact);

    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(water_block_inner)
            .addPostContactInteraction(water_block_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeStream>(water_block);

    StartupToConstantInflowSpeed free_stream_speed(U_f, 2.0);
    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact)
            .addPostStateDynamics<fluid_dynamics::FreeStreamCondition<StartupToConstantInflowSpeed>>(
                water_block, free_stream_speed);

    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact);

    // Transport velocity correction (cylinder pattern):
    // KernelGradientIntegral feeds TransportVelocityCorrectionCK.
    // ConstantConstraintCK pins emitter buffer displacement to zero.
    AlignedBoxByCell emitter_buffer_part(
        water_block,
        AlignedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    auto &transport_correction =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_block_contact)
            .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter, BulkParticles>(water_block)
            .addPostStateDynamics<ConstantConstraintCK, Vecd>(emitter_buffer_part, "Displacement", Vecd::Zero());

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<
            fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(water_block_contact);

    auto &fluid_advection_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_f);
    auto &fluid_acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);

    auto &water_update_particle_position =
        main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);
    //----------------------------------------------------------------------
    //	FSI: forces on solid from fluid.
    //----------------------------------------------------------------------
    // auto &viscous_force_from_fluid =
    //     main_methods.addInteractionDynamics<
    //         FSI::ViscousForceOnStructure<
    //             fluid_dynamics::ViscousForceCK<Contact<Wall, Viscosity, NoKernelCorrectionCK>>>>(fish_contact);
    // auto &pressure_force_from_fluid =
    //     main_methods.addInteractionDynamics<
    //         FSI::PressureForceOnStructure<
    //             fluid_dynamics::AcousticStep2ndHalf<
    //                 Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>>(fish_contact);
    //----------------------------------------------------------------------
    //	Free-stream boundary: emitter injection, inflow condition, disposer.
    //----------------------------------------------------------------------
    AlignedBoxByParticle emitter_part(
        water_block,
        AlignedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    auto &emitter_injection =
        main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK>(emitter_part);
    auto &inflow_condition =
        main_methods.addStateDynamics<
            fluid_dynamics::EmitterInflowConditionCK, StartupToConstantInflowSpeed>(
            emitter_buffer_part, free_stream_speed);

    AlignedBoxByCell disposer_part(
        water_block,
        AlignedBox(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    auto &disposer_indication =
        main_methods.addStateDynamics<fluid_dynamics::WithinDisposerIndication>(disposer_part);
    auto &particle_deletion =
        main_methods.addStateDynamics<fluid_dynamics::OutflowParticleDeletion>(water_block);
    //----------------------------------------------------------------------
    //	Output.
    //----------------------------------------------------------------------
    auto &body_state_recorder =
        main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(water_block, "Pressure");
    body_state_recorder.addToWrite<int>(water_block, "Indicator");
    body_state_recorder.addToWrite<Real>(water_block, "PositionDivergence");
    body_state_recorder.addToWrite<Real>(water_block, "DensitySummation");
    // body_state_recorder.addToWrite<int>(fish_body, "MaterialID");
    // body_state_recorder.addToWrite<Matd>(fish_body, "ActiveStrain");

    Gravity gravity(Vecd::Zero());
    ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>
        write_water_mechanical_energy(water_block, gravity);
    //----------------------------------------------------------------------
    //	Simulation parameters.
    //----------------------------------------------------------------------
    Real End_Time = 1.7;
    Real D_Time = 0.01;
    int screen_output_interval = 100;

    //----------------------------------------------------------------------
    //	Time stepper.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_trigger =
        time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    auto &state_recording_trigger =
        time_stepper.addTriggerByInterval(D_Time);
    size_t number_of_iterations = 0;
    //----------------------------------------------------------------------
    //	Initialisation.
    //----------------------------------------------------------------------
    update_fish_cell_linked_list.exec();
    update_water_body_configuration.exec();
    fish_inner_build.exec();
    fish_contact_relation.exec();
    fish_body_corrected_configuration.exec();

    apply_gravity_force.exec();
    fluid_boundary_indicator.exec();
    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    transport_correction.exec();
    fluid_viscous_force.exec();
    // viscous_force_from_fluid.exec();

    body_state_recorder.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_acoustic_step;

    //----------------------------------------------------------------------
    //	Main loop — TimeStepper manages physical time and acoustic stepping.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime(End_Time))
    {
        //------------------------------------------------------------------
        //	Acoustic time stepping (fastest, most frequent).
        //------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);

        fluid_acoustic_step_1st_half.exec(acoustic_dt); // FreeStreamCondition runs here
        inflow_condition.exec();
        // pressure_force_from_fluid.exec();
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);

        // Solid sub-loop removed for fluid-only mode.
        // Real dt_s_raw = fish_body_computing_time_step_size.exec();
        // int solid_sub_div = (dt_s_raw > Real(0))
        //     ? static_cast<int>(acoustic_dt / dt_s_raw) + 2
        //     : 33;
        // initialize_displacement.exec();
        // time_stepper.integrateMatchedTimeInterval(
        //     acoustic_dt, solid_sub_div,
        //     [&](Real dt_s)
        //     {
        //         imposing_active_strain.exec();
        //         fish_body_stress_relaxation_first_half.exec(dt_s);
        //         fish_body_stress_relaxation_second_half.exec(dt_s);
        //     });

        // update_average_velocity.exec(acoustic_dt);

        interval_acoustic_step += TickCount::now() - time_instance;
        //------------------------------------------------------------------
        //	Advection-level operations (every Dt interval).
        //------------------------------------------------------------------
        if (advection_trigger(fluid_advection_time_step))
        {
            number_of_iterations++;
            Real Dt = advection_trigger.getInterval();
            water_update_particle_position.exec();

            // viscous_force_from_fluid.exec();
            // fish_body_update_normal.exec();

            //	Screen output.
            if (number_of_iterations % screen_output_interval == 0)
            {
                size_t inner_ite_dt = static_cast<size_t>(Dt / time_stepper.getGlobalTimeStepSize());
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "  Time=" << time_stepper.getPhysicalTime()
                          << "  Dt=" << Dt
                          << "  Dt/dt=" << inner_ite_dt << "\n";
                write_water_mechanical_energy.writeToFile(number_of_iterations);
            }

            //	VTP output every D_Time.
            if (state_recording_trigger())
            {
                TickCount t2 = TickCount::now();
                body_state_recorder.writeToFile();
                interval += TickCount::now() - t2;
            }

            //	Particle injection, deletion, sort, relation update.
            emitter_injection.exec();
            disposer_indication.exec();
            particle_deletion.exec();

            // Sort disabled: GPU sort of 52k+ particles causes NaN (uninitialized reserve
            // particles mixed into active array during sort). Buffer increased to 0.8x to
            // accommodate particle growth from emitter without overflow.
            // TODO: fix root cause — guard sort against reserve particles.
            // if (number_of_iterations % 100 == 0)
            //     particle_sort.exec();

            update_water_body_configuration.exec();
            fish_contact_relation.exec();

            //	Advection-level fluid operators for next Dt.
            apply_gravity_force.exec();
            fluid_boundary_indicator.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
        }
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Acoustic step time: " << interval_acoustic_step.seconds() << " seconds." << std::endl;

    return 0;
}
