/**
 * @file    2d_flow_stream_around_fish.cpp
 * @brief   fish swimming driven by active muscles
 * @author  Yaru Ren and Xiangyu Hu
 */
#include "2d_flow_stream_around_fish.h"
#include "sphinxsys.h"
using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    system.setRunParticleRelaxation(true);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    system.setReloadParticles(false);
    // handle command line arguments
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief   Particles and body creation for fish.
     */
    SolidBody fish_body(system, makeShared<FishBody>("FishBody"));
    fish_body.defineAdaptationRatios(1.15, 2.0);
    fish_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    fish_body.defineParticlesAndMaterial<ElasticSolidParticles, FishBodyComposite>();
    //  Using relaxed particle distribution if needed
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? fish_body.generateParticles<ParticleGeneratorReload>(io_environment, fish_body.getName())
        : fish_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation fish_inner(fish_body);
    InnerRelation water_block_inner(water_block);
    ComplexRelation water_block_complex(water_block, {&fish_body});
    ContactRelation fish_contact(fish_body, {&water_block});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        /**
         * @brief 	Methods used for particle relaxation.
         */
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_fish_body_particles(fish_body);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_fish_body(io_environment, fish_body);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&fish_body});

        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(fish_inner);
        /**
         * @brief 	Particle relaxation starts here.
         */
        random_fish_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_fish_body.writeToFile();

        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_fish_body.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;

        /** Output results. */
        write_particle_reload_files.writeToFile();
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero()));
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    /** Emitter buffer inflow condition. */
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer);
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);
    /** time-space method to detect surface particles. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
        free_stream_surface_indicator(water_block_complex);
    /** Evaluation of density by freestream approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density(water_block_complex);
    /** We can output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** modify the velocity of boundary particles with free-stream velocity. */
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
    /** Pressure relaxation using verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Impose transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<IndicatedParticles<0>>> transport_velocity_correction(water_block_complex);
    /** Computing vorticity in the flow. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> fish_body_normal_direction(fish_body);
    /** Corrected configuration for the elastic insert body. */
    InteractionWithUpdate<CorrectedConfigurationInner> fish_body_corrected_configuration(fish_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(fish_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluidRiemann>
        fluid_force_on_fish_update(fish_contact, viscous_force_on_solid);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(fish_body);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    SimpleDynamics<FishMaterialInitialization> composite_material_id(fish_body);
    SimpleDynamics<ImposingActiveStrain> imposing_active_strain(fish_body);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> fish_body_computing_time_step_size(fish_body);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> fish_body_stress_relaxation_first_half(fish_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> fish_body_stress_relaxation_second_half(fish_inner);
    /** Update norm .*/
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> fish_body_update_normal(fish_body);
    fish_body.addBodyStateForRecording<Real>("Density");
    fish_body.addBodyStateForRecording<int>("MaterialID");
    fish_body.addBodyStateForRecording<Matd>("ActiveStrain");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    RestartIO restart_io(io_environment, system.real_bodies_);
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_viscous_force_on_insert_body(io_environment, viscous_force_on_solid, "TotalViscousForceOnSolid");
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_force_on_insert_body(io_environment, fluid_force_on_fish_update, "TotalForceOnSolid");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    system.initializeSystemConfigurations();
    /** computing surface normal direction for the fish. */
    fish_body_normal_direction.exec();
    /** computing linear reproducing configuration for the fish. */
    fish_body_corrected_configuration.exec();
    /** initialize material ids for the fish. */
    composite_material_id.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real End_Time = 2.5; /**< End time. */
    Real D_Time = 0.01;  /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;

        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_solid.exec();
            /** Update normal direction on elastic body.*/
            fish_body_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = get_fluid_time_step_size.exec();
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** FSI for fluid force on solid body. */
                fluid_force_on_fish_update.exec();
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
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
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                emitter_buffer_inflow_condition.exec(dt);
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            fish_body.updateCellLinkedList();
            /** one need update configuration after periodic condition. */
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            fish_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_on_insert_body.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
