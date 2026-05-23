/**
 * @file 2d_flow_stream_around_fish_sycl.cpp
 * @brief fish swimming driven by active muscles — SYCL CK direct API (filling_tank pattern)
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
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody fish_body(sph_system, makeShared<FishBody>("FishBody"));
    fish_body.defineAdaptationRatios(1.15, 2.0);
    fish_body.defineMaterial<FishBodyComposite>();
    fish_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    Inner<> fish_inner(fish_body, ConfigType::Lagrangian);
    Inner<> water_block_inner(water_block);
    Contact<> water_block_contact(water_block, {&fish_body});
    Contact<> fish_contact(fish_body, {&water_block});
    Contact<> fish_fsi_contact(fish_body, {&water_block});
    //----------------------------------------------------------------------
    //	Configuration dynamics: cell linked lists, relations, sorting.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> fish_cell_linked_list(fish_body);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>>
        water_body_update_complex_relation(water_block_inner, water_block_contact);
    // fish_inner is Lagrangian — built once during init, never rebuilt in the loop.
    UpdateRelation<MainExecutionPolicy, Inner<>> fish_inner_build(fish_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> fish_contact_relation(fish_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fish_fsi_contact_relation(fish_fsi_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);
    //----------------------------------------------------------------------
    //	Solid initialisation (CPU — no CK version needed for one-shot ops).
    //----------------------------------------------------------------------
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK>
        fish_body_normal_direction(fish_body);
    StateDynamics<execution::ParallelPolicy, FishMaterialInitialization>
        composite_material_id(fish_body);
    //----------------------------------------------------------------------
    //	Solid correction matrix (GPU, Lagrangian — built once, stays fixed).
    //----------------------------------------------------------------------
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        fish_body_corrected_configuration(fish_inner);
    //----------------------------------------------------------------------
    //	GPU replacements for AverageVelocityAndAcceleration.
    //----------------------------------------------------------------------
    StateDynamics<MainExecutionPolicy, InitializeDisplacementCK>
        initialize_displacement(fish_body);
    StateDynamics<MainExecutionPolicy, UpdateAverageVelocityAndAccelerationCK>
        update_average_velocity(fish_body);
    //----------------------------------------------------------------------
    //	Free-stream surface indicator (GPU).
    //----------------------------------------------------------------------
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationCK<Inner<WithUpdate>>>
        free_stream_surface_indicator(water_block_inner);
    free_stream_surface_indicator.addPostContactInteraction(water_block_contact);
    //----------------------------------------------------------------------
    //	Solid dynamics (GPU).
    //----------------------------------------------------------------------
    StateDynamics<MainExecutionPolicy, ImposingActiveStrain>
        imposing_active_strain(fish_body);

    InteractionDynamicsCK<MainExecutionPolicy,
                          solid_dynamics::StructureIntegration1stHalfPK2<Inner<OneLevel, FishBodyComposite>>>
        fish_body_stress_relaxation_first_half(fish_inner);

    InteractionDynamicsCK<MainExecutionPolicy,
                          solid_dynamics::StructureIntegration2ndHalf<Inner<OneLevel>>>
        fish_body_stress_relaxation_second_half(fish_inner);

    ReduceDynamicsCK<MainExecutionPolicy, solid_dynamics::AcousticTimeStepCK>
        fish_body_computing_time_step_size(fish_body);

    StateDynamics<MainExecutionPolicy, ZeroForceCK> zero_force(fish_body);

    StateDynamics<MainExecutionPolicy, solid_dynamics::UpdateElasticNormalDirectionCK>
        fish_body_update_normal(fish_body);
    //----------------------------------------------------------------------
    //	Gravity (GPU).
    //----------------------------------------------------------------------
    StartupAcceleration time_dependent_acceleration(Vec2d::Zero(), 2.0);
    StateDynamics<MainExecutionPolicy, GravityForceCK<StartupAcceleration>>
        apply_gravity_force(water_block, time_dependent_acceleration);
    //----------------------------------------------------------------------
    //	Fluid dynamics (GPU).
    //	Construction order matters for variable registration — see comments in
    //	the original CPU version.
    //----------------------------------------------------------------------
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>>
        fluid_linear_correction_matrix(water_block_inner);
    fluid_linear_correction_matrix.addPostContactInteraction(water_block_contact);

    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup>
        water_advection_step_setup(water_block);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
    pressure_relaxation(water_block_inner, water_block_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
    density_relaxation(water_block_inner, water_block_contact);

    // Wall-contact operator for FSI type extraction only — not called in time loop
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalf<Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>>>
        density_relaxation_wall(water_block_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensitySummationCK<Inner<>>>
        update_fluid_density(water_block_inner);
    update_fluid_density.addPostContactInteraction(water_block_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, FreeSurface>>
        fluid_density_regularization(water_block);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        viscous_force(water_block_inner, water_block_contact);

    // Wall-contact operators for FSI type extraction only — not called in time loop
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::ViscousForceCK<Contact<Wall, Viscosity, NoKernelCorrectionCK>>>
        viscous_force_wall(water_block_contact);

    InteractionDynamicsCK<MainExecutionPolicy,
                          KernelGradientIntegral<Inner<NoKernelCorrectionCK>>>
        water_kernel_gradient_integral(water_block_inner);
    water_kernel_gradient_integral.addPostContactInteraction<Boundary, NoKernelCorrectionCK>(water_block_contact);

    StateDynamics<MainExecutionPolicy,
                  fluid_dynamics::TransportVelocityCorrectionCK<SPHBody, TruncatedLinear, BulkParticles>>
        transport_velocity_correction(water_block);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK>
        get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>>
        get_fluid_time_step_size(water_block);

    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition>
        water_update_particle_position(water_block);
    //----------------------------------------------------------------------
    //	Free-stream boundary conditions.
    //----------------------------------------------------------------------
    AlignedBoxByParticle emitter(
        water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    StateDynamics<MainExecutionPolicy, fluid_dynamics::EmitterInflowInjectionCK<AlignedBoxByParticle>>
        emitter_inflow_injection(emitter);

    AlignedBoxByCell emitter_buffer(
        water_block, AlignedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    StateDynamics<execution::ParallelPolicy,
                  fluid_dynamics::EmitterInflowConditionCK<AlignedBoxByCell, FreeStreamVelocity>>
        emitter_buffer_inflow_condition(emitter_buffer);

    // DisposerOutflowDeletion has no UpdateKernel — use CPU SimpleDynamics directly.
    AlignedBoxByCell disposer(
        water_block, AlignedBox(xAxis, Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer);

    // FreeStreamVelocityCorrection has no CK version — use CPU SimpleDynamics directly.
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>>
        velocity_boundary_condition_constraint(water_block);
    //----------------------------------------------------------------------
    //	FSI (GPU).
    //	Use fish_fsi_contact (plain Contact<>) because ViscousForceOnStructure
    //	and PressureForceOnStructure constructors require Contact<>.
    //----------------------------------------------------------------------
    using VFWallType = std::remove_reference_t<decltype(viscous_force_wall)>;
    InteractionDynamicsCK<MainExecutionPolicy, FSI::ViscousForceOnStructure<VFWallType>>
        viscous_force_from_fluid(fish_fsi_contact);

    using DRWallType = std::remove_reference_t<decltype(density_relaxation_wall)>;
    InteractionDynamicsCK<MainExecutionPolicy, FSI::PressureForceOnStructure<DRWallType>>
        pressure_force_from_fluid(fish_fsi_contact);
    //----------------------------------------------------------------------
    //	Define output methods.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_state_recorder(sph_system);
    body_state_recorder.addToWrite<Real>(water_block, "Pressure");
    body_state_recorder.addToWrite<int>(water_block, "Indicator");
    body_state_recorder.addToWrite<int>(fish_body, "MaterialID");
    body_state_recorder.addToWrite<Matd>(fish_body, "ActiveStrain");

    Gravity gravity(Vecd::Zero());
    ReducedQuantityRecording<TotalMechanicalEnergy>
        write_water_mechanical_energy(water_block, gravity);
    ReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_viscous_force(fish_body, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_pressure_force(fish_body, "PressureForceFromFluid");
    //----------------------------------------------------------------------
    //	Prepare the simulation.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");

    fish_body_normal_direction.exec();
    composite_material_id.exec();

    water_cell_linked_list.exec();
    fish_cell_linked_list.exec();
    water_body_update_complex_relation.exec();
    fish_inner_build.exec();   // Lagrangian — built once, never rebuilt
    fish_contact_relation.exec();
    fish_fsi_contact_relation.exec();

    fish_body_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Time-stepping setup.
    //----------------------------------------------------------------------
    Real End_Time = 1.7;
    Real D_Time = 0.01;
    int screen_output_interval = 100;
    size_t number_of_iterations = 0;

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_acoustic_step;

    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < D_Time)
        {
            apply_gravity_force.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();

            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();
            fluid_density_regularization.exec();
            viscous_force.exec();
            water_kernel_gradient_integral.exec();
            transport_velocity_correction.exec();
            viscous_force_from_fluid.exec();
            fish_body_update_normal.exec();

            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                TickCount time_instance = TickCount::now();
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);

                pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec(dt);
                density_relaxation.exec(dt);
                velocity_boundary_condition_constraint.exec();
                emitter_buffer_inflow_condition.exec(dt);

                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                initialize_displacement.exec();
                while (dt_s_sum < dt)
                {
                    zero_force.exec();
                    Real dt_s = SMIN(fish_body_computing_time_step_size.exec(), dt - dt_s_sum);
                    if (dt_s <= Real(0))
                        break;
                    imposing_active_strain.exec();
                    fish_body_stress_relaxation_first_half.exec(dt_s);
                    fish_body_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                update_average_velocity.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
                interval_acoustic_step += TickCount::now() - time_instance;
                inner_ite_dt++;
            }

            water_update_particle_position.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "  Time=" << sv_physical_time->getValue()
                          << "  Dt=" << Dt
                          << "  Dt/dt=" << inner_ite_dt
                          << "  dt/dt_s=" << inner_ite_dt_s << "\n";
                write_total_viscous_force.writeToFile(number_of_iterations);
                write_water_mechanical_energy.writeToFile(number_of_iterations);
                write_total_pressure_force.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                particle_sort.exec();

            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fish_cell_linked_list.exec();
            fish_contact_relation.exec();
            fish_fsi_contact_relation.exec();
        }

        TickCount t2 = TickCount::now();
        body_state_recorder.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Acoustic step time: " << interval_acoustic_step.seconds() << " seconds." << std::endl;

    return 0;
}
