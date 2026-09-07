#include "test_2d_turbulent_wavy_channel.h"
using namespace SPH;

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRestartStep(0);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.Name())
        : water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineMatterMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.Name())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer_center_point(sph_system, "ObserverCenterPoint");
    observer_center_point.generateParticles<ObserverParticles>(observer_location_center_point);

    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation observer_centerpoint_contact(observer_center_point, {&water_block});

    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);

    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;

        InnerRelation wall_boundary_inner(wall_boundary);

        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles_water(water_block);

        BodyStatesRecordingToVtp write_inserted_body_to_vtp(wall_boundary);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block);

        ReloadParticleIO write_particle_reload_files(SPHBodyVector{&water_block, &wall_boundary});

        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);

        random_inserted_body_particles.exec(0.25);
        random_inserted_body_particles_water.exec(0.25);

        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_inner_water.SurfaceBounding().exec();

        write_inserted_body_to_vtp.writeToFile(0);
        write_inserted_body_to_vtp_water.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            relaxation_step_inner_water.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
                write_inserted_body_to_vtp_water.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of the wall_boundary finish !" << std::endl;
        std::cout << "The physics relaxation process of the water_block finish !" << std::endl;

        write_particle_reload_files.writeToFile(0);

        return 0;
    }

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::udf::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact, y_p_constant);
    InteractionWithUpdate<fluid_dynamics::udf::TurbulentLinearGradientCorrectionMatrixComplex> corrected_configuration_fluid_separated_inner_wall(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerNoRiemann> density_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::udf::Integration2ndHalfOnlyWallAcousticRiemannAdjusted> density_relaxation_wall(water_wall_contact);
    density_relaxation.post_processes_.push_back(&density_relaxation_wall);
    InteractionWithUpdate<fluid_dynamics::udf::kOmega_GetVelocityGradientComplex> get_velocity_gradient(water_block_inner, water_wall_contact);
    SimpleDynamics<fluid_dynamics::udf::kOmega_kTransportEquationInner> k_equation_relaxation(water_block_inner, initial_turbu_values, is_AMRD, is_blended);
    InteractionDynamics<fluid_dynamics::udf::kOmega_TKE_Diffusion> compute_TKE_diffusion(water_block_inner);
    SimpleDynamics<fluid_dynamics::udf::kOmega_omegaTransportEquationInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::udf::kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner> compute_TSDR_diffusion_and_gradient_k_omega(water_block_inner);
    InteractionDynamics<fluid_dynamics::udf::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::udf::kOmega_WallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::udf::P_refinement_GetVelocityGradientInner> get_velocity_gradient_inner_only_for_P(water_block_inner);
    SimpleDynamics<fluid_dynamics::udf::P_refinement<num_node_sublayer_model,type_tdma_sublayer_model>> get_friction_velocity_from_sublayer(water_block, y_p_constant);
    InteractionWithUpdate<fluid_dynamics::udf::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::udf::TVC_ModifiedLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    IncreaseToFullGravity time_dependent_acceleration(external_acc, external_acc_gradually_impose_t);
    SimpleDynamics<GravityForce<Gravity>> apply_gravity_force(water_block, time_dependent_acceleration);
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, periodic_along_x);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_fluid_density(water_block_inner, water_wall_contact);
    SimpleDynamics<UpdateVolume> update_volume(water_block);
    ReduceDynamics<fluid_dynamics::udf::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    SimpleDynamics<fluid_dynamics::udf::kOmegaTurbulentEddyViscosity> update_eddy_viscosity(water_block);
    ParticleSorting particle_sorting(water_block);

    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    ObservedQuantityRecording<Real> write_centerpoint_quantity("TurbulentViscosity", observer_centerpoint_contact);

    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 150.0;
    Real num_output_files = 4.0;
    Real Output_Time = end_time / num_output_files;
    Real dt = 0.0;

    TickCount t1 = TickCount::now();
    TimeInterval interval;

    wall_boundary_normal_direction.exec();
    inlet_outlet_surface_particle_indicator.exec();
    update_near_wall_status.exec();
    corrected_configuration_fluid.exec();
    corrected_configuration_fluid_separated_inner_wall.exec();
    get_velocity_gradient.exec();
    get_velocity_gradient_inner_only_for_P.exec();
    update_eddy_viscosity.exec();

    body_states_recording.writeToFile();
    write_centerpoint_quantity.writeToFile(number_of_iterations);

    int num_output_file = 0;
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            apply_gravity_force.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            update_volume.exec();
            corrected_configuration_fluid.exec();
            corrected_configuration_fluid_separated_inner_wall.exec();
            if (physical_time > turbulent_module_activate_time)
            {
                update_eddy_viscosity.exec();
                update_near_wall_status.exec();
                standard_wall_function_correction.exec();
                get_velocity_gradient_inner_only_for_P.exec();
                get_friction_velocity_from_sublayer.exec();
            }
            turbulent_viscous_force.exec();
            if (physical_time > turbulent_module_activate_time)
            {
                get_velocity_gradient.exec();
                compute_TKE_diffusion.exec();
                compute_TSDR_diffusion_and_gradient_k_omega.exec();
            }
            transport_velocity_correction.exec();
            kernel_summation.exec();
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                if (physical_time > turbulent_module_activate_time)
                {
                    turbulent_kinetic_energy_force.exec();
                }
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                if (physical_time > turbulent_module_activate_time)
                {
                    k_equation_relaxation.exec(dt);
                    epsilon_equation_relaxation.exec(dt);
                }
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_itr++;
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerpoint_quantity.writeToFile(number_of_iterations);
                }
            }

            number_of_iterations++;
            periodic_condition_x.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();
            inlet_outlet_surface_particle_indicator.exec();
        }
        body_states_recording.writeToFile();
        observer_centerpoint_contact.updateConfiguration();
        num_output_file++;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    // if (sph_system.GenerateRegressionData())
    // {
    //     write_centerpoint_quantity.generateDataBase(1.0e-3);
    // }
    // else
    // {
    //     write_centerpoint_quantity.testResult();
    // }
    return 0;
}
