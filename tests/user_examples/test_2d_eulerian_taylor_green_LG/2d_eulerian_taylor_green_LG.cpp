/**
 * @file 	eulerian_taylor_green_LG.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation with the Laguerre Gauss kernel.
 * @details 2D eulerian_taylor_green vortex flow example.
 * @author 	Chi Zhang, Zhentong Wang and Xiangyu Hu
 */
#include "2d_eulerian_taylor_green_LG.h"
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
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    IOEnvironment io_environment(sph_system);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    EulerianFluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.sph_adaptation_->resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    water_body.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_f, heat_capacity_ratio, mu_f);
    water_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_body);
    InteractionWithUpdate<KernelGradientWithCorrectionInner> kernel_gradient_update(water_body_inner);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Initial condition with momentum and energy field */
    SimpleDynamics<TaylorGreenInitialCondition> initial_condition(water_body);
    /** Initialize a time step. */
    SimpleDynamics<EulerianCompressibleTimeStepInitialization> time_step_initialization(water_body);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
    /** Periodic BCs in y direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<EulerianCompressibleAcousticTimeStepSize> get_fluid_time_step_size(water_body);
    /** Pressure relaxation algorithm by using verlet time stepping. */
    /** Here, we can use HLLC with Limiter Riemann solver for pressure relaxation and density and energy relaxation  */
    InteractionWithUpdate<Integration1stHalfHLLCWithLimiterRiemann> pressure_relaxation(water_body_inner);
    InteractionWithUpdate<Integration2ndHalfHLLCWithLimiterRiemann> density_and_energy_relaxation(water_body_inner);
    /** Computing viscous acceleration. */
    InteractionDynamics<EulerianCompressibleViscousAccelerationInner> viscous_acceleration(water_body_inner);
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
    kernel_gradient_update.exec();
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
