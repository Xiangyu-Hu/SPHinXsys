/**
 * @file 2d_flow_stream_around_fish_sycl.cpp
 * @brief Fish swimming driven by active muscles — full SYCL/CK GPU port.
 *        Fluid operators follow mr_free_stream_around_cylinder_sycl as reference.
 *        Solid FSI operators use the same ParticleMethodContainer API.
 *        Output format matches the CPU fish case for direct comparison.
 * @author Pruthvik Arasikere Mallikarjuna and Xiangyu Hu
 */
#include "2d_flow_stream_around_fish.h"
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include "sphinxsys.h"
#include <algorithm>
#include <cmath>
using namespace SPH;

int main(int ac, char *av[])
{
    const paper_bench::BenchmarkDefaults benchmark_defaults{
        0.0025, 1.7, 0.01,
        {{"standard", 0.0025}}};
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av, benchmark_defaults);
    configureParticleSpacing(Real(benchmark_config.dp));
    const Vec2d fish_tail_observation_point = fishTailObservationPoint();
    TickCount benchmark_start = TickCount::now();
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "fish_fsi_2d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    TimeInterval io_interval = TickCount::now() - environment_io_start;
    const double environment_io_seconds = io_interval.seconds();
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.8);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody fish_body(sph_system, makeShared<FishBody>("FishBody"));
    fish_body.defineAdaptationRatios(1.15, 2.0);
    fish_body.defineMatterMaterial<FishBodyComposite>();
    fish_body.generateParticles<BaseParticles, Lattice>();
    fish_body.getBaseParticles().registerStateVariable<Vecd>("Velocity");

    ObserverBody fish_tail_observer(sph_system, "FishTailObserver");
    fish_tail_observer.defineAdaptationRatios(1.15, 2.0);
    fish_tail_observer.generateParticles<ObserverParticles>(
        StdVec<Vecd>{fish_tail_observation_point});
    const size_t initial_particle_count =
        water_block.getBaseParticles().TotalRealParticles() +
        fish_body.getBaseParticles().TotalRealParticles() +
        fish_tail_observer.getBaseParticles().TotalRealParticles();
    //----------------------------------------------------------------------
    //	Define body relations.
    //----------------------------------------------------------------------
    Inner<> fish_inner(fish_body, ConfigType::Lagrangian);
    Inner<> water_block_inner(water_block);
    Contact<> water_block_contact(water_block, {&fish_body}); // fish as stationary wall
    Contact<> fish_contact(fish_body, {&water_block});
    Contact<> fish_tail_observer_contact(fish_tail_observer, {&fish_body});
    //----------------------------------------------------------------------
    //	Define SPH solver — all operators via ParticleMethodContainer.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);

    // CPU-only one-shot initialisation (no GPU CK version needed for these).
    auto &host_methods = sph_solver.getHostMethodContainer();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(fish_body).exec();
    host_methods.addStateDynamics<FishMaterialInitialization>(fish_body).exec();

    // All GPU operators via main_methods (par_ck = MainExecutionPolicy).
    auto &main_methods = sph_solver.getMainMethodContainer();
    //----------------------------------------------------------------------
    //	Configuration dynamics.
    //----------------------------------------------------------------------
    auto &update_fish_cell_linked_list =
        main_methods.addCellLinkedListDynamics(fish_body);
    auto &update_water_cell_linked_list =
        main_methods.addCellLinkedListDynamics(water_block);
    auto &water_body_update_complex_relation =
        main_methods.addRelationDynamics(water_block_inner, water_block_contact);

    // fish_inner is Lagrangian — built once during init, never rebuilt.
    auto &fish_inner_build =
        main_methods.addRelationDynamics(fish_inner);
    auto &fish_contact_relation =
        main_methods.addRelationDynamics(fish_contact);
    auto &fish_tail_observer_contact_relation =
        main_methods.addRelationDynamics(fish_tail_observer_contact);
    //----------------------------------------------------------------------
    //	Solid correction matrix (Lagrangian — computed once, stays fixed).
    //----------------------------------------------------------------------
    auto &fish_body_corrected_configuration =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(fish_inner);
    //----------------------------------------------------------------------
    //	FSI average velocity / acceleration for dummy particle boundary condition.
    //----------------------------------------------------------------------
    auto &initialize_displacement =
        main_methods.addStateDynamics<InitializeDisplacementCK>(fish_body);
    auto &update_average_velocity =
        main_methods.addStateDynamics<UpdateAverageVelocityAndAccelerationCK>(fish_body);
    //----------------------------------------------------------------------
    //	Solid dynamics.
    //----------------------------------------------------------------------
    auto &imposing_active_strain =
        main_methods.addStateDynamics<ImposingActiveStrain>(fish_body);
    auto &fish_body_numerical_damping =
        main_methods.addInteractionDynamicsWithUpdate<
            solid_dynamics::StructureNumericalDamping, FishBodyComposite>(fish_inner);
    auto &fish_body_stress_relaxation_first_half =
        main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration1stHalfPK2, FishBodyComposite>(fish_inner);
    auto &fish_body_stress_relaxation_second_half =
        main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration2ndHalf>(fish_inner);
    auto &fish_body_computing_time_step_size =
        main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(fish_body);
    auto &fish_body_update_normal =
        main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(fish_body);
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

    StartupToConstantInflowSpeed free_stream_speed(0.0, 2.0);
    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact);

    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_contact);

    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(water_block_inner)
            .addPostContactInteraction(water_block_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeStream>(water_block);

    // Transport velocity correction (cylinder pattern):
    // KernelGradientIntegral feeds TransportVelocityCorrectionCK.
    // ConstantConstraintCK pins emitter buffer displacement to zero.
    OrientedBoxByCell emitter_buffer_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
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
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(water_block);

    auto &water_update_particle_position =
        main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);
    //----------------------------------------------------------------------
    //	FSI: forces on solid from fluid.
    //----------------------------------------------------------------------
    auto &viscous_force_from_fluid =
        main_methods.addInteractionDynamics<
            FSI::ViscousForceOnStructure<
                fluid_dynamics::ViscousForceCK<Contact<Wall, Viscosity, NoKernelCorrectionCK>>>>(fish_contact);
    auto &pressure_force_from_fluid =
        main_methods.addInteractionDynamics<
            FSI::PressureForceOnStructure<
                fluid_dynamics::AcousticStep2ndHalf<
                    Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>>>(fish_contact);
    //----------------------------------------------------------------------
    //	Free-stream boundary: emitter injection, inflow condition, disposer.
    //----------------------------------------------------------------------
    OrientedBoxByParticle emitter_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    // Sort is constructed AFTER emitter_part so body_parts_by_particle_ includes it.
    // This mirrors the cylinder pattern and ensures sort.exec() updates the emitter's particle list.
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);
    auto &emitter_injection = main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK>(emitter_part);
    auto &inflow_condition =
        main_methods.addStateDynamics<fluid_dynamics::EmitterInflowConditionCK, StartupToConstantInflowSpeed>(
            emitter_buffer_part, free_stream_speed);
    auto &free_stream_condition = main_methods.addStateDynamics<
        fluid_dynamics::FreeStreamCondition<StartupToConstantInflowSpeed>>(water_block, free_stream_speed);

    Vec2d correct_disposer_translation = Vec2d(DL, -0.25 * DH) + disposer_halfsize;
    OrientedBoxByCell disposer_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(correct_disposer_translation)), disposer_halfsize));
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
    body_state_recorder.addToWrite<int>(fish_body, "MaterialID");
    body_state_recorder.addToWrite<Matd>(fish_body, "ActiveStrain");

    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);

    Gravity gravity(Vecd::Zero());
    ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>
        write_water_mechanical_energy(water_block, gravity);
    auto &write_fish_tail_position =
        main_methods.addObserveRecorder<Vecd>(
            fish_tail_observer_contact, "Position");
    //----------------------------------------------------------------------
    //	Simulation parameters.
    //----------------------------------------------------------------------
    Real End_Time = Real(benchmark_config.end_time);
    Real D_Time = Real(benchmark_config.output_interval);
    int screen_output_interval = 100;
    const bool verification_mode = !benchmark_config.benchmark_mode;
    const bool write_heavy_output = benchmark_config.output_enabled;
    size_t fish_tail_samples = 0;

    //----------------------------------------------------------------------
    //	Time stepper.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_trigger =
        time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    auto &state_recording_trigger =
        time_stepper.addTriggerByInterval(D_Time);
    size_t number_of_iterations = 0;
    size_t number_of_advection_steps = 0;
    size_t number_of_acoustic_steps = 0;
    size_t number_of_solid_steps = 0;
    size_t number_of_output_steps = 0;
    //----------------------------------------------------------------------
    //	Initialisation.
    //----------------------------------------------------------------------
    update_fish_cell_linked_list.exec();
    update_water_cell_linked_list.exec();
    water_body_update_complex_relation.exec();
    fish_inner_build.exec();
    fish_contact_relation.exec();
    fish_tail_observer_contact_relation.exec();
    fish_body_corrected_configuration.exec();

    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        TickCount restart_io_start = TickCount::now();
        restart_io.readRestartFiles(sph_system.RestartStep());
        io_interval += TickCount::now() - restart_io_start;
        // rebuild spatial hash from restored positions
        update_fish_cell_linked_list.exec();
        update_water_cell_linked_list.exec();
        // re-sort water particles for cache locality, then rebuild CLL
        particle_sort.exec();
        update_water_cell_linked_list.exec();
        // reconstruct neighbor lists
        water_body_update_complex_relation.exec();
        fish_contact_relation.exec();
        fish_tail_observer_contact_relation.exec();
        // rebuild B-matrix from restored fish positions
        fish_body_corrected_configuration.exec();
        // repopulate density from restored positions
        fluid_density_regularization.exec();
        number_of_iterations = sph_system.RestartStep();
    }

    apply_gravity_force.exec();
    fluid_boundary_indicator.exec();
    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    transport_correction.exec();
    fluid_viscous_force.exec();
    viscous_force_from_fluid.exec();

    const double init_seconds = std::max(
        0.0, (TickCount::now() - initialization_start).seconds() -
                 (io_interval.seconds() - environment_io_seconds));

    {
        TickCount initial_io_start = TickCount::now();
        fish_tail_observer_contact_relation.exec();
        write_fish_tail_position.writeToFile(number_of_iterations);
        fish_tail_samples++;
        if (write_heavy_output)
        {
            body_state_recorder.writeToFile();
            write_water_mechanical_energy.writeToFile(number_of_iterations);
        }
        io_interval += TickCount::now() - initial_io_start;
    }

    TimeInterval interval_acoustic_step;
    bool simulation_failed = false;
    std::string failure_reason;

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
        number_of_acoustic_steps++;

        fluid_acoustic_step_1st_half.exec(acoustic_dt); // FreeStreamCondition runs here
        inflow_condition.exec();
        free_stream_condition.exec();
        pressure_force_from_fluid.exec();
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);

        initialize_displacement.exec();

        Real dt_s_sum = 0.0;
        int guard = 0;
        int solid_sub_count = 0;

        while (dt_s_sum < acoustic_dt && guard++ < 10000)
        {
            Real dt_s_candidate =
                fish_body_computing_time_step_size.exec();

            if (!std::isfinite(dt_s_candidate) || dt_s_candidate <= 0.0)
            {
                std::cerr << "ERROR: invalid solid time step\n";
                simulation_failed = true;
                failure_reason = "non-positive or non-finite solid time step";
                break;
            }

            Real dt_s =
                std::min(dt_s_candidate,
                         acoustic_dt - dt_s_sum);

            imposing_active_strain.exec();

            fish_body_numerical_damping.exec(dt_s);
            fish_body_stress_relaxation_first_half.exec(dt_s);
            fish_body_stress_relaxation_second_half.exec(dt_s);

            dt_s_sum += dt_s;
            solid_sub_count++;
            number_of_solid_steps++;
        }

        if (!simulation_failed && dt_s_sum < acoustic_dt)
        {
            std::cerr << "ERROR: solid substep guard exhausted\n";
            simulation_failed = true;
            failure_reason = "solid substep guard exhausted";
        }
        if (simulation_failed)
        {
            interval_acoustic_step += TickCount::now() - time_instance;
            break;
        }

        update_average_velocity.exec(acoustic_dt);

        interval_acoustic_step += TickCount::now() - time_instance;
        //------------------------------------------------------------------
        //	Advection-level operations (every Dt interval).
        //------------------------------------------------------------------
        if (advection_trigger(fluid_advection_time_step))
        {
            number_of_iterations++;
            number_of_advection_steps++;
            Real Dt = advection_trigger.getInterval();
            water_update_particle_position.exec();

            viscous_force_from_fluid.exec();
            fish_body_update_normal.exec();

            //	Screen output.
            if (number_of_iterations % screen_output_interval == 0)
            {
                TickCount console_io_start = TickCount::now();
                size_t inner_ite_dt = static_cast<size_t>(Dt / time_stepper.getGlobalTimeStepSize());
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "  Time=" << time_stepper.getPhysicalTime()
                          << "  Dt=" << Dt
                          << "  Dt/dt=" << inner_ite_dt
                          << "  dt/dt_s=" << solid_sub_count << "\n";
                if (verification_mode && write_heavy_output)
                {
                    write_water_mechanical_energy.writeToFile(
                        number_of_iterations);
                }
                io_interval += TickCount::now() - console_io_start;
            }

            //	VTP output every D_Time.
            if (state_recording_trigger())
            {
                TickCount periodic_io_start = TickCount::now();
                fish_tail_observer_contact_relation.exec();
                write_fish_tail_position.writeToFile(number_of_iterations);
                fish_tail_samples++;
                if (write_heavy_output)
                {
                    body_state_recorder.writeToFile();
                }
                io_interval += TickCount::now() - periodic_io_start;
                number_of_output_steps++;
            }

            if (verification_mode && write_heavy_output &&
                number_of_iterations % 500 == 0)
            {
                TickCount restart_io_start = TickCount::now();
                restart_io.writeToFile(number_of_iterations);
                io_interval += TickCount::now() - restart_io_start;
            }

            //	Particle injection, deletion, sort, relation update.
            emitter_injection.exec();
            disposer_indication.exec();
            particle_deletion.exec();

            if (number_of_iterations % 100 == 0)
                particle_sort.exec();

            update_fish_cell_linked_list.exec();
            update_water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fish_contact_relation.exec();
            fish_tail_observer_contact_relation.exec();

            //	Advection-level fluid operators for next Dt.
            apply_gravity_force.exec();
            fluid_boundary_indicator.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
        }
    }

    const TimeInterval total_wall = TickCount::now() - benchmark_start;
    const double compute_seconds = std::max(
        0.0, total_wall.seconds() - init_seconds -
                 io_interval.seconds());
    std::cout << "Total wall time: " << total_wall.seconds()
              << " seconds (compute: " << compute_seconds
              << ", initialization: " << init_seconds
              << ", I/O: " << io_interval.seconds() << ")." << std::endl;
    std::cout << "Acoustic step time: " << interval_acoustic_step.seconds() << " seconds." << std::endl;

    paper_bench::BenchmarkSummary summary;
    summary.dp = particle_spacing_ref;
    summary.initial_particle_count = initial_particle_count;
    summary.particle_count =
        water_block.getBaseParticles().TotalRealParticles() +
        fish_body.getBaseParticles().TotalRealParticles() +
        fish_tail_observer.getBaseParticles().TotalRealParticles();
    summary.physical_time = time_stepper.getPhysicalTime();
    summary.outer_steps = number_of_output_steps;
    summary.advection_steps = number_of_advection_steps;
    summary.acoustic_steps = number_of_acoustic_steps;
    summary.wall_seconds = total_wall.seconds();
    summary.compute_seconds = compute_seconds;
    summary.io_seconds = io_interval.seconds();
    summary.init_seconds = init_seconds;
    summary.time_per_outer_step =
        number_of_output_steps == 0
            ? 0.0
            : compute_seconds / static_cast<double>(number_of_output_steps);
    summary.status = simulation_failed ? "failed" : "completed";
    summary.solid_steps = number_of_solid_steps;
    summary.setParticleUpdates(
        summary.particle_count * summary.advection_steps,
        "particle_count*advection_steps");
    summary.acoustic_component_seconds =
        paper_bench::formatMetric(interval_acoustic_step.seconds());
    summary.extra_fields.push_back(
        {"fish_tail_reference_x",
         std::to_string(fish_tail_observation_point[0])});
    summary.extra_fields.push_back(
        {"fish_tail_reference_y",
         std::to_string(fish_tail_observation_point[1])});
    summary.extra_fields.push_back(
        {"fish_tail_samples", std::to_string(fish_tail_samples)});
    summary.extra_fields.push_back(
        {"failure_reason", failure_reason});
    summary.extra_fields.push_back(
        {"time_per_advection_step_seconds",
         number_of_advection_steps == 0
             ? "0"
             : std::to_string(
                   compute_seconds /
                   static_cast<double>(number_of_advection_steps))});
    summary.extra_fields.push_back(
        {"outer_step_definition", "periodic feature/output trigger"});
    summary.extra_fields.push_back(
        {"feature_definition",
         "interpolated Position at the dp-adjusted tail-center reference point"});
    benchmark_recorder.writeSummary(summary);

    return simulation_failed ? 1 : 0;
}
