/**
 * @file 	freestream_flow_around_square.cpp
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "case.h"
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
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody square(sph_system, makeShared<Square>("square"));
    square.defineAdaptationRatios(1.15, 2.0);
    square.defineBodyLevelSetShape();
    square.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? square.generateParticles<ParticleGeneratorReload>(io_environment, square.getName())
        : square.generateParticles<ParticleGeneratorLattice>();

    SolidBody beam(sph_system, makeShared<DefaultShape>("beam"));
    beam.defineAdaptationRatios(1.15, particle_spacing_ref / particle_spacing_beam);
    beam.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam.generateParticles<BeamParticleGenerator>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation square_inner(square);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(square);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, {&square});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&square});
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(square_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation beam_inner(beam);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> beam_curvature(beam_inner);
    InnerRelation water_block_inner(water_block);
    ComplexRelation water_block_complex(water_block_inner, RealBodyVector{&square, &beam});
    ContactRelation beam_water_contact(beam, {&water_block});
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
    /** Pressure relaxation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
    /** Density relaxation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Apply transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_complex);
    /** compute the vorticity. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> square_normal_direction(square);
    /** Corrected configuration for the elastic insert body. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> beam_corrected_configuration(beam_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_solid(beam_water_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdate>
        fluid_force_on_solid_update(beam_water_contact, viscous_force_on_solid);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(beam);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> beam_computing_time_step_size(beam);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> beam_stress_relaxation_first_half(beam_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> beam_stress_relaxation_second_half(beam_inner);
    /** Constrain region of the inserted body. */
    BoundaryGeometry beam_boundary_geometry(beam, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_beam_base(beam_boundary_geometry);
    /** damping */
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        beam_position_damping(0.2, beam_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        beam_rotation_damping(0.2, beam_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    beam.addBodyStateForRecording<Real>("MeanCurvature");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the insert body. */
    square_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    beam_corrected_configuration.exec();
    beam_curvature.compute_initial_curvature();
    water_block_complex.getContactRelation().updateConfiguration();
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

            /** FSI for viscous force. */
            viscous_force_on_solid.exec();

            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                fluid_force_on_solid_update.exec();
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(beam_computing_time_step_size.exec(), dt - dt_s_sum);
                    beam_stress_relaxation_first_half.exec(dt_s);
                    constraint_beam_base.exec();
                    beam_position_damping.exec();
                    beam_rotation_damping.exec();
                    beam_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

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
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt
                          << "	dt / dt_s = " << inner_ite_dt_s
                          << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            beam_curvature.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            beam.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            beam_water_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
