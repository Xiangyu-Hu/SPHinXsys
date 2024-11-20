/**
 * @file 	cohesive_soil_failure.cpp
 * @brief 	2D cohesive soil failure.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#include "cohesive_soil_failure.h" //SPHinXsys Library.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    RealBody soil_block(sph_system, makeShared<Soil>("GranularBody"));
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle, cohesion);
    soil_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation soil_block_inner(soil_block);
    ContactRelation soil_block_contact(soil_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation soil_block_complex(soil_block_inner, soil_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(soil_block, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<SoilInitialCondition> soil_initial_condition(soil_block);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> correction_matrix(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration1stHalfWithWallRiemann> granular_stress_relaxation(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration2ndHalfWithWallRiemann> granular_density_relaxation(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> soil_density_by_summation(soil_block_inner, soil_block_contact);
    InteractionDynamics<continuum_dynamics::StressDiffusion> stress_diffusion(soil_block_inner);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<FreeSurfaceNormalComplex> free_surface_normal(soil_block_inner, soil_block_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> soil_acoustic_time_step(soil_block, 0.4);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(soil_block, "Pressure");
    body_states_recording.addToWrite<Real>(soil_block, "Density");
    SimpleDynamics<continuum_dynamics::VerticalStress> vertical_stress(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "VerticalStress");
    SimpleDynamics<continuum_dynamics::AccDeviatoricPlasticStrain> accumulated_deviatoric_plastic_strain(soil_block);
    body_states_recording.addToWrite<Real>(soil_block, "AccDeviatoricPlasticStrain");
    body_states_recording.addToWrite<int>(soil_block, "Indicator");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_mechanical_energy(soil_block, gravity);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    soil_initial_condition.exec();
    correction_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real End_Time = 2.0;         /**< End time. */
    Real D_Time = End_Time / 50; /**< Time stamps for output of body states. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_soil_stress_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();

            soil_density_by_summation.exec();
            surface_indicator.exec();
            free_surface_normal.exec();
            transport_velocity_correction.exec();
            Real dt = soil_acoustic_time_step.exec();
            stress_diffusion.exec();
            granular_stress_relaxation.exec(dt);
            granular_density_relaxation.exec(dt);
            integration_time += dt;
            physical_time += dt;

            interval_computing_soil_stress_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4) << "	Time = "
                          << physical_time
                          << std::scientific << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_mechanical_energy.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            time_instance = TickCount::now();
            /** Update cell linked list and configuration. */
            soil_block.updateCellLinkedList();
            soil_block_complex.updateConfiguration();
            correction_matrix.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        vertical_stress.exec();
        accumulated_deviatoric_plastic_strain.exec();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << std::fixed << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_soil_stress_relaxation = "
              << interval_computing_soil_stress_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_mechanical_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_mechanical_energy.testResult();
    }

    return 0;
};
