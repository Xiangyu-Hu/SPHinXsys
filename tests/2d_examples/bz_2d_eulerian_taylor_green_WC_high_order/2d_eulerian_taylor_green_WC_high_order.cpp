/**
 * @file 	eulerian_taylor_green.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_taylor_green vortex flow example.
 * @author 	Chi Zhang, Zhentong Wang and Xiangyu Hu
 */
#include "2d_eulerian_taylor_green_WC_high_order.h"
#include "sphinxsys.h"
using namespace SPH; //	Namespace cite here.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); //Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);       //Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);// handle command line arguemnts.
#endif
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    EulerianFluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
        : water_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map.
        //----------------------------------------------------------------------
        InnerRelation water_body_inner(water_body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_body);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_water_body_to_vtp(io_environment, { &water_body });
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, { &water_body });

        /* Relaxation method: including based on the 0th and 1st order consistency. */
        relax_dynamics::RelaxationStepInner relaxation_0th_inner(water_body_inner, false);
        relax_dynamics::RelaxationStepImplicitInner relaxation_0th_implicit_inner(water_body_inner, false);
        InteractionDynamics<relax_dynamics::CalculateCorrectionMatrix> calculate_correction_matrix(water_body_inner, false);
        relax_dynamics::RelaxationStepByCMInner relaxation_1st_inner(water_body_inner, false);
        relax_dynamics::RelaxationStepByCMImplicitInner relaxation_1st_implicit_inner(water_body_inner, false);

        PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
        PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);

        InteractionDynamics<relax_dynamics::CheckCorrectedZeroOrderConsistency> check_corrected_zero_order_consistency(water_body_inner);
        ReduceAverage<QuantitySummation<Real>> calculate_particle_average_zero_error(water_body, "CorrectedZeroOrderErrorNorm");
        ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_zero_error(water_body, "CorrectedZeroOrderErrorNorm");
        InteractionDynamics<relax_dynamics::CheckCorrectedFirstOrderConsistency> check_corrected_first_order_consistency(water_body_inner);
        ReduceAverage<QuantitySummation<Real>> calculate_particle_average_first_error(water_body, "CorrectedFirstOrderErrorNorm");
        ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_first_error(water_body, "CorrectedFirstOrderErrorNorm");

        //----------------------------------------------------------------------  
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_water_body_particles.exec(0.15);
        sph_system.initializeSystemCellLinkedLists();
        periodic_condition_x.update_cell_linked_list_.exec();
        periodic_condition_y.update_cell_linked_list_.exec();
        sph_system.initializeSystemConfigurations();
        write_water_body_to_vtp.writeToFile(0);

        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        TickCount t1 = TickCount::now();

        int ite = 0; //iteration step for the total relaxation step.

        Real last_zero_maximum_residual = 1;
        Real last_zero_average_residual = 1;
        Real last_first_maximum_residual = 1;
        Real last_first_average_residual = 1;

        Real current_zero_maximum_residual = 1;
        Real current_zero_average_residual = 1;
        Real current_first_maximum_residaul = 1;
        Real current_first_average_residaul = 1;

        GlobalStaticVariables::physical_time_ = ite;
        /* The procedure to obtain uniform particle distribution that satisfies the 0ht order consistency. */
        while (current_zero_average_residual > 0.0001)
        {
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            water_body.updateCellLinkedList();
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            water_body_inner.updateConfiguration();

            relaxation_0th_implicit_inner.exec(0.1);
            //calculate_correction_matrix.exec();
            //relaxation_1st_implicit_inner.exec(0.1);

            ite++;

            if (ite % 100 == 0)
            {
                check_corrected_zero_order_consistency.exec();
                current_zero_average_residual = calculate_particle_average_zero_error.exec();
                current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();

                std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
                std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << ": average = " << current_zero_average_residual << std::endl;
                write_water_body_to_vtp.writeToFile(ite);
            }
        }

        check_corrected_zero_order_consistency.exec();
        current_zero_average_residual = calculate_particle_average_zero_error.exec();
        current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
        std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << ": average = " << current_zero_average_residual << std::endl;

        ite++;
        write_water_body_to_vtp.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);

        TickCount t2 = TickCount::now();
        TickCount::interval_t tt;
        tt = t2 - t1;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_body);
    InteractionWithUpdate<KernelGradientWithCorrectionInner> kernel_gradient_update(water_body_inner);
    InteractionWithUpdate<ConfigurationInner> configuration_fluid(water_body_inner);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Initial condition with momentum and energy field */
    SimpleDynamics<TaylorGreenInitialCondition> initial_condition(water_body);
    /** Initialize a time step. */
    SimpleDynamics<EulerianWCTimeStepInitialization> time_step_initialization(water_body);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
    /** Periodic BCs in y direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<EulerianWCAcousticTimeStepSize> get_fluid_time_step_size(water_body);
    /** Pressure relaxation algorithm by using verlet time stepping. */
    /** Here, we can use HLLC with Limiter Riemann solver for pressure relaxation and density and energy relaxation  */
    InteractionWithUpdate<Integration1stHalfAcousticRiemann> pressure_relaxation(water_body_inner, 0);
    InteractionWithUpdate<Integration2ndHalfAcousticRiemann> density_and_energy_relaxation(water_body_inner, 0);
    /** Computing viscous acceleration. */
    InteractionDynamics<WCEulerianViscousAccelerationInner> viscous_acceleration(water_body_inner);
    water_body.addBodyStateForRecording<Real>("Pressure");
    water_body.addBodyStateForRecording<Matd>("CorrectionMatrix");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    /** Output the mechanical energy of fluid body. */
    RegressionTestEnsembleAverage<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_total_mechanical_energy(io_environment, water_body);
    /** Output the maximum speed of the fluid body. */
    RegressionTestEnsembleAverage<ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>>>
        write_maximum_speed(io_environment, water_body);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    //kernel_gradient_update.exec();
    configuration_fluid.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 5.0;
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    /** Output the mechanical energy of fluid. */
    write_total_mechanical_energy.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force. */
            time_step_initialization.exec();
            Real dt = get_fluid_time_step_size.exec();
            viscous_acceleration.exec();
            /** Dynamics including pressure relaxation. */
            integration_time += dt;
            pressure_relaxation.exec(dt);
            density_and_energy_relaxation.exec(dt);
            GlobalStaticVariables::physical_time_ += dt;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_total_mechanical_energy.writeToFile(number_of_iterations);
        write_maximum_speed.writeToFile(number_of_iterations);
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    write_total_mechanical_energy.testResult();
    write_maximum_speed.testResult();

    return 0;
}
