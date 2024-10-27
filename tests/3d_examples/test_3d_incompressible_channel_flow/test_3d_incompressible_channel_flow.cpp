/**
 * @file 	test_3d_incompressible_channel_flow.cpp
 * @brief 	This is the inviscid incompressible channel flow test for the realization of FVM in the SPHinXsys.
 * @author 	Yash Mandaokar, Zhentong Wang and Xiangyu Hu
 */
#include "test_3d_incompressible_channel_flow.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh read_mesh_data(mesh_fullpath);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, read_mesh_data.MinMeshEdge());
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody air_block(sph_system, makeShared<AirBody>("AirBody"));
    air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    air_block.generateParticlesWithReserve<BaseParticles, UnstructuredMesh>(ghost_boundary, read_mesh_data);
    GhostCreationFromMesh ghost_creation(air_block, read_mesh_data, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM air_block_inner(air_block, read_mesh_data);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<InvCFInitialCondition> initial_condition(air_block);
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
     * the value is larger, the numerical dissipation larger. */
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfInnerRiemann> pressure_relaxation(air_block_inner, 500.0);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfInnerRiemann> density_relaxation(air_block_inner, 8000.0);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(air_block, read_mesh_data.MinMeshEdge(), 0.6);
    /** Boundary conditions set up */
    InvCFBoundaryConditionSetup boundary_condition_setup(air_block_inner, ghost_creation);
    //----------------------------------------------------------------------
    BodyStatesRecordingInMeshToVtu write_real_body_states(air_block, read_mesh_data);
    write_real_body_states.addToWrite<Real>(air_block, "Density");
    write_real_body_states.addToWrite<Real>(air_block, "Pressure");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(air_block);

    air_block_inner.updateConfiguration();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 30;
    Real output_interval = 2; /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            boundary_condition_setup.resetBoundaryConditions();
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);

            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
                write_maximum_speed.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            // write_real_body_states.writeToFile();
        }
        TickCount t2 = TickCount::now();
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
