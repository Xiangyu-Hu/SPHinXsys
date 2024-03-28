/**
 * @file 	test_3d_incompressible_channel_flow.cpp
 * @brief 	This is the inviscid incompressible channel flow test for the realization of FVM in the SPHinXsys.
 * @author 	Yash Mandaokar, Zhentong Wang and Xiangyu Hu
 */
#include "test_3d_incompressible_channel_flow.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh_3d read_mesh_data(mesh_fullpath);
    //----------------------------------------------------------------------
   //	Build up the environment of a SPHSystem.
   //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody air_block(sph_system, makeShared<WaveBody>("AirBody"));
    air_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    air_block.generateParticlesWithReserve<UnstructuredMesh_3d>(ghost_boundary, read_mesh_data);
    air_block.addBodyStateForRecording<Real>("Density");
    air_block.addBodyStateForRecording<Real>("Pressure");
    SimpleDynamics<InvCFInitialCondition> initial_condition(air_block);
    GhostCreationFromMesh_3d ghost_creation(air_block, read_mesh_data, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM_3d air_block_inner(air_block, read_mesh_data);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
        /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfInnerRiemann> pressure_relaxation(air_block_inner, 500.0);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfInnerRiemann> density_relaxation(air_block_inner, 8000.0);
    /** Boundary conditions set up */
    InvCFBoundaryConditionSetup boundary_condition_setup(air_block_inner, ghost_creation);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(air_block, read_mesh_data.min_distance_between_nodes_,0.6);
    //----------------------------------------------------------------------
    // Visualization in FVM with date in cell.
    BodyStatesRecordingInMeshToVtu write_real_body_states(air_block, read_mesh_data);
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(air_block);

    air_block_inner.updateConfiguration();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
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
    while (GlobalStaticVariables::physical_time_ < end_time)
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
            GlobalStaticVariables::physical_time_ += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
                    << GlobalStaticVariables::physical_time_
                    << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            //write_real_body_states.writeToFile();
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
    return 0;
}