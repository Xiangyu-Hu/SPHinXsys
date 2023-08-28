/**
 * @file 	shock_tube.cpp
 * @brief 	This is a test to show the standard Sod shock tube case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "shock_tube.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Create body, materials and particles.
    //----------------------------------------------------------------------
    EulerianFluidBody wave_body(sph_system, makeShared<WaveBlock>("WaveBody"));
    wave_body.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_l, heat_capacity_ratio);
    wave_body.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The inner relation defines the particle configuration for particles within a body.
    //	The contact relation defines the particle configuration between the bodies.
    //----------------------------------------------------------------------
    InnerRelation wave_body_inner(wave_body);
    // using correction Matrix to correct kernel gradient.
    InteractionWithUpdate<KernelGradientWithCorrectionInner> kernel_gradient_update(wave_body_inner);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<ShockTubeInitialCondition> waves_initial_condition(wave_body);
    wave_body.addBodyStateForRecording<Real>("TotalEnergy");
    wave_body.addBodyStateForRecording<Real>("Density");
    // Initialize particle acceleration.
    SimpleDynamics<EulerianCompressibleTimeStepInitialization> initialize_wave_step(wave_body);
    // Periodic BCs in y direction.
    PeriodicConditionUsingCellLinkedList periodic_condition_y(wave_body, wave_body.getBodyShapeBounds(), yAxis);
    // Time step size with considering sound wave speed.
    ReduceDynamics<EulerianCompressibleAcousticTimeStepSize> get_wave_time_step_size(wave_body);
    // Pressure, density and energy relaxation algorithm by use HLLC Riemann solver.
    InteractionWithUpdate<Integration1stHalfHLLCRiemann> pressure_relaxation(wave_body_inner);
    InteractionWithUpdate<Integration2ndHalfHLLCRiemann> density_and_energy_relaxation(wave_body_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations of the simulation.
    //	Regression tests are also defined here.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt body_states_recording(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    waves_initial_condition.exec();
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.2;
    Real output_interval = 0.01; /**< Time stamps for output of body states. */
    //----------------------------------------------------------------------
    // Output the start states of bodies.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Statistics for computing CPU time.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        //	Integrate time (loop) until the next output time.
        while (integration_time < output_interval)
        {
            initialize_wave_step.exec();
            Real dt = get_wave_time_step_size.exec();
            // Dynamics including pressure and density and energy relaxation.
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
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
