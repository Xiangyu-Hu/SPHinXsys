/**
 * @file 2d_flow_stream_around_fish.cpp
 * @brief fish swimming driven by active muscles
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
    auto &water_block = sph_system.addBody<FluidBody>(makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    auto &fish_body = sph_system.addBody<SolidBody>(makeShared<FishBody>("FishBody"));
    fish_body.defineAdaptationRatios(1.15, 2.0);
    fish_body.defineMaterial<FishBodyComposite>();
    fish_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    auto &fish_inner = sph_system.addInnerRelation(fish_body, ConfigType::Lagrangian);
    auto &water_block_inner = sph_system.addInnerRelation(water_block);
    auto &water_block_contact = sph_system.addContactRelation(water_block, fish_body);
    auto &fish_contact = sph_system.addContactRelation(fish_body, water_block);
    //----------------------------------------------------------------------
    //	Define SPH solver with host (CPU) and main (GPU) method containers.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    // host_methods: runs on CPU — used for operators with no CK version yet
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    // main_methods: dispatches to GPU via SYCL (CPU fallback if no GPU found)
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    //----------------------------------------------------------------------
    //	Define numerical methods.
    //	Configuration dynamics first, then physics, then auxiliary.
    //----------------------------------------------------------------------

    // --- Cell linked lists and relations (GPU) ---
    auto &particle_sort = main_methods.addSortDynamics(water_block);

    ParticleDynamicsGroup update_water_configuration;
    update_water_configuration.add(&main_methods.addCellLinkedListDynamics(water_block));
    update_water_configuration.add(&main_methods.addRelationDynamics(water_block_inner, water_block_contact));

    ParticleDynamicsGroup update_fish_configuration;
    update_fish_configuration.add(&main_methods.addCellLinkedListDynamics(fish_body));
    update_fish_configuration.add(&main_methods.addRelationDynamics(fish_inner));
    update_fish_configuration.add(&main_methods.addRelationDynamics(fish_contact));

    // --- Solid initialisation (CPU: no CK versions needed for one-shot ops) ---
    auto &fish_body_normal_direction =
        host_methods.addStateDynamics<NormalFromBodyShapeCK>(fish_body);
    auto &composite_material_id =
        host_methods.addStateDynamics<SimpleDynamics<FishMaterialInitialization>>(fish_body);

    // --- Solid correction matrix (GPU, Lagrangian — built once, stays fixed) ---
    auto &fish_body_corrected_configuration =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(fish_inner);

    // --- Free-stream surface indicator (GPU) ---
    auto &free_stream_surface_indicator =
        main_methods.addInteractionDynamics<fluid_dynamics::SpatialTemporalFreeSurfaceIndicationCK>(water_block_inner)
            .addPostContactInteraction(water_block_contact);

    // --- Solid dynamics (GPU) ---
    auto &imposing_active_strain =
        main_methods.addStateDynamics<ImposingActiveStrain>(fish_body);

    auto &fish_body_stress_relaxation_first_half =
        main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration1stHalfPK2, FishBodyComposite>(fish_inner);

    auto &fish_body_stress_relaxation_second_half =
        main_methods.addInteractionDynamicsOneLevel<solid_dynamics::StructureIntegration2ndHalf>(fish_inner);

    auto &fish_body_computing_time_step_size =
        main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(fish_body);

    auto &fish_body_update_normal =
        main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(fish_body);

    // --- Gravity (GPU) ---
    StartupAcceleration time_dependent_acceleration(Vec2d::Zero(), 2.0);
    auto &apply_gravity_force =
        main_methods.addStateDynamics<GravityForceCK<StartupAcceleration>>(water_block, time_dependent_acceleration);

    // --- Fluid dynamics (GPU) ---
    auto &update_fluid_density =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(water_block_inner)
            .addPostContactInteraction(water_block_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(water_block);

    auto &viscous_force =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForce, WithUpdate>(water_block_inner)
            .addPostContactInteraction<Wall>(water_block_contact);

    auto &transport_velocity_correction =
        main_methods.addInteractionDynamics<fluid_dynamics::TransportVelocityCorrectionCK, BulkParticles>(water_block_inner)
            .addPostContactInteraction(water_block_contact);

    auto &pressure_relaxation =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_contact);

    auto &density_relaxation =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_contact);

    auto &get_fluid_advection_time_step_size =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_f);
    auto &get_fluid_time_step_size =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);

    // --- Free-stream boundary conditions (CPU — no CK versions available) ---
    AlignedBoxByParticle emitter(
        water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    auto &emitter_inflow_injection =
        host_methods.addStateDynamics<EmitterInflowInjectionCK>(emitter, inlet_particle_buffer);

    AlignedBoxByCell emitter_buffer(
        water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    auto &emitter_buffer_inflow_condition =
        host_methods.addStateDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>>(emitter_buffer);

    AlignedBoxByCell disposer(
        water_block, AlignedBox(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    auto &disposer_outflow_deletion =
        host_methods.addStateDynamics<fluid_dynamics::DisposerOutflowDeletion>(disposer);

    // FreeStreamVelocityCorrection has no CK version — runs on CPU
    auto &velocity_boundary_condition_constraint =
        host_methods.addStateDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>>(water_block);

    // --- FSI (GPU) ---
    // ViscousForceOnStructure and PressureForceOnStructure are CK versions from shared_ck/
    auto &viscous_force_from_fluid =
        main_methods.addInteractionDynamics<FSI::ViscousForceOnStructure<decltype(viscous_force)>>(fish_contact);

    auto &pressure_force_from_fluid =
        main_methods.addInteractionDynamics<FSI::PressureForceOnStructure<decltype(density_relaxation)>>(fish_contact);

    // AverageVelocityAndAcceleration has no CK version — runs on CPU
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(fish_body);
    //----------------------------------------------------------------------
    //	Define output methods.
    //----------------------------------------------------------------------
    auto &body_state_recorder =
        main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(water_block, "Pressure");
    body_state_recorder.addToWrite<int>(water_block, "Indicator");
    body_state_recorder.addToWrite<int>(fish_body, "MaterialID");
    body_state_recorder.addToWrite<Matd>(fish_body, "ActiveStrain");

    Gravity gravity(Vecd::Zero());
    auto &write_water_mechanical_energy =
        main_methods.addObserveRegression<RegressionTestDynamicTimeWarping, Real>(
            "TotalMechanicalEnergy", water_block);

    auto &write_total_viscous_force =
        main_methods.addObserve<ReducedQuantityRecording<QuantitySummation<Vecd>>>(
            fish_body, "ViscousForceFromFluid");

    auto &write_total_pressure_force =
        main_methods.addObserve<ReducedQuantityRecording<QuantitySummation<Vecd>>>(
            fish_body, "PressureForceFromFluid");
    //----------------------------------------------------------------------
    //	Prepare the simulation.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    fish_body_normal_direction.exec();
    composite_material_id.exec();
    fish_body_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Time-stepping setup.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    Real End_Time = 1.7;
    Real D_Time = 0.01;
    int screen_output_interval = 100;
    size_t number_of_iterations = 0;

    auto &state_recording = time_stepper.addTriggerByInterval(D_Time);

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_acoustic_step;

    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime(End_Time))
    {
        apply_gravity_force.exec();
        Real Dt = get_fluid_advection_time_step_size.exec();

        free_stream_surface_indicator.exec();
        update_fluid_density.exec();
        viscous_force.exec();
        transport_velocity_correction.exec();

        /** FSI: viscous force on fish from fluid */
        viscous_force_from_fluid.exec();
        /** Update fish surface normals after deformation */
        fish_body_update_normal.exec();

        size_t inner_ite_dt = 0;
        size_t inner_ite_dt_s = 0;
        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            TickCount time_instance = TickCount::now();
            Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);

            /** Fluid pressure relaxation, first half */
            pressure_relaxation.exec(dt);
            /** FSI: pressure force on fish from fluid */
            pressure_force_from_fluid.exec(dt);
            /** Fluid pressure relaxation, second half */
            density_relaxation.exec(dt);
            /** Free-stream velocity correction on surface particles (CPU) */
            velocity_boundary_condition_constraint.exec();
            emitter_buffer_inflow_condition.exec(dt);

            /** Solid sub-stepping */
            inner_ite_dt_s = 0;
            Real dt_s_sum = 0.0;
            average_velocity_and_acceleration.initialize_displacement_.exec();
            while (dt_s_sum < dt)
            {
                Real dt_s = SMIN(fish_body_computing_time_step_size.exec(), dt - dt_s_sum);
                imposing_active_strain.exec();
                fish_body_stress_relaxation_first_half.exec(dt_s);
                fish_body_stress_relaxation_second_half.exec(dt_s);
                dt_s_sum += dt_s;
                inner_ite_dt_s++;
            }
            average_velocity_and_acceleration.update_averages_.exec(dt);

            relaxation_time += dt;
            interval_acoustic_step += TickCount::now() - time_instance;
            inner_ite_dt++;
        }

        if (number_of_iterations % screen_output_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9)
                      << "N=" << number_of_iterations
                      << "  Time=" << time_stepper.getPhysicalTime()
                      << "  Dt=" << Dt
                      << "  Dt/dt=" << inner_ite_dt
                      << "  dt/dt_s=" << inner_ite_dt_s << "\n";
            write_total_viscous_force.writeToFile(number_of_iterations);
            write_water_mechanical_energy.writeToFile(number_of_iterations);
            write_total_pressure_force.writeToFile(number_of_iterations);
        }
        number_of_iterations++;

        /** Inject new particles at inlet, delete at outlet */
        emitter_inflow_injection.exec();
        disposer_outflow_deletion.exec();

        /** Periodic particle sorting for cache efficiency */
        if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            particle_sort.exec();

        /** Rebuild cell linked lists and neighbour relations */
        update_water_configuration.exec();
        update_fish_configuration.exec();

        if (state_recording.exec())
        {
            body_state_recorder.writeToFile();
        }
    }

    TimeInterval tt = TickCount::now() - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Acoustic step time: " << interval_acoustic_step.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
        write_water_mechanical_energy.generateDataBase(0.3);
    else if (sph_system.RestartStep() == 0)
        write_water_mechanical_energy.testResult();

    return 0;
}
