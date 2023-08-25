/**
 * @file 	mr_freestream_flow_around_cylinder.cpp
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "mr_free_stream_around_cylinder.h"
#include "sphinxsys.h"

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<ParticleRefinementWithinShape>(1.3, 1.0, 1);
    water_block.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(io_environment);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    MultiPolygonShape refinement_region(MultiPolygon(initial_refinement_region), "RefinementRegion");
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
        : water_block.generateParticles<ParticleGeneratorMultiResolution>(refinement_region);
    water_block.addBodyStateForRecording<Real>("SmoothingLengthRatio");

    SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 4.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
        : cylinder.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    AdaptiveInnerRelation water_block_inner(water_block);
    AdaptiveContactRelation water_contact(water_block, {&cylinder});
    ComplexRelation water_block_complex(water_block_inner, water_contact);
    AdaptiveContactRelation cylinder_contact(cylinder, {&water_block});
    AdaptiveContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
        BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
        relax_dynamics::RelaxationStepComplex relaxation_step_complex(water_block_complex, "OuterBoundary", true);
        SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(water_block, refinement_region);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            update_smoothing_length_ratio.exec();
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finished !" << std::endl;
        /** Output results. */
        write_real_body_particle_reload_files.writeToFile(0);
        return 0;
    }

    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero()));
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    /** Emitter buffer inflow condition. */
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer, 0.1);
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);
    /** time-space method to detect surface particles. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
        free_stream_surface_indicator(water_block_complex);
    /** Evaluation of density by freestream approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplexAdaptive> update_fluid_density(water_block_complex);
    /** We can output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("VolumetricMeasure");
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** modify the velocity of boundary particles with free-stream velocity. */
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
    /** Pressure relaxation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
    /** Density relaxation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Apply transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplexAdaptive>
        transport_velocity_correction(water_block_complex);
    /** compute the vorticity. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_pressure_force_on_inserted_body(cylinder_contact);
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> fluid_viscous_force_on_inserted_body(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    ObservedQuantityRecording<Vecd>
        write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_viscous_force_on_inserted_body(io_environment, fluid_viscous_force_on_inserted_body, "TotalViscousForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>
        write_total_force_on_inserted_body(io_environment, fluid_pressure_force_on_inserted_body, "TotalPressureForceOnSolid");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the insert body. */
    cylinder_normal_direction.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 200.0;
    Real output_interval = end_time / 400.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                emitter_buffer_inflow_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            /** write run-time observation into file */
            cylinder_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
        write_total_force_on_inserted_body.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
