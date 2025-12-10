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
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<AdaptiveWithinShape>(1.3, 1.0, 1);
    water_block.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(sph_system);
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    MultiPolygonShape refinement_region(MultiPolygon(initial_refinement_region), "RefinementRegion");
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice, AdaptiveByShape>(refinement_region);

    SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 4.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();

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
    AdaptiveInnerRelation water_block_inner(water_block);
    AdaptiveContactRelation water_contact(water_block, {&cylinder});
    AdaptiveContactRelation cylinder_contact(cylinder, {&water_block});
    AdaptiveContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which are only use for now for updating configurations.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_contact);
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
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_real_body_particle_reload_files({&water_block, &cylinder});
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(cylinder_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            DynamicsArgs(water_block_inner, std::string("OuterBoundary")), water_contact);
        SimpleDynamics<UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(water_block, refinement_region);
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
    StartupAcceleration time_dependent_acceleration(Vec2d(U_f, 0.0), 2.0);
    SimpleDynamics<GravityForce<StartupAcceleration>> apply_gravity_force(water_block, time_dependent_acceleration);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplexAdaptive> update_fluid_density(water_block_inner, water_contact);

    AlignedBoxByParticle emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);
    AlignedBoxByCell emitter_buffer(water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer, 0.1);
    AlignedBoxByCell disposer(water_block, AlignedBox(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer);

    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    /** modify the velocity of boundary particles with free-stream velocity. */
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
    /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_contact);
    /** Apply transport velocity formulation. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplexAdaptive<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_contact);
    /** compute the vorticity. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    write_real_body_states.addToWrite<Real>(water_block, "VolumetricMeasure");
    write_real_body_states.addToWrite<Real>(water_block, "SmoothingLengthRatio");
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_from_fluid(cylinder, "PressureForceFromFluid");
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
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            apply_gravity_force.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();
            viscous_force.exec();
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
                physical_time += dt;
                emitter_buffer_inflow_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";

                viscous_force_from_fluid.exec();
                pressure_force_from_fluid.exec();
                write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
                write_total_pressure_force_from_fluid.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            /** write run-time observation into file */
            cylinder_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        // The lift force at the cylinder is very small and not important in this case.
        write_total_viscous_force_from_fluid.generateDataBase(1.0e-3);
    }
    else
    {
        write_total_viscous_force_from_fluid.testResult();
    }

    return 0;
}
