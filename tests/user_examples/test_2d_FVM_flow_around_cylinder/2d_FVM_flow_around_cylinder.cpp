/**
 * @file 	2d_FVM_flow_around_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder in FVM.
 * @details We consider a flow passing by a cylinder in 2D in FVM framework.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#include "2d_FVM_flow_around_cylinder.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANASYS mesh.file
    readMeshFile read_mesh_data(zero_three_flow_around_cylinder_mesh_file_fullpath);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    EulerianFluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorInFVM>(read_mesh_data.elements_center_coordinates_, read_mesh_data.elements_volumes_);
    water_block.addBodyStateForRecording<Real>("Density");
    /** Initial condition */
    SimpleDynamics<WeaklyCompressibleFluidInitialCondition> initial_condition(water_block);
    GhostCreationFromMesh ghost_creation(water_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(water_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
    water_block_inner.updateConfiguration();
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Boundary conditions set up */
    FACBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation.each_boundary_type_with_all_ghosts_index_,
                                                       ghost_creation.each_boundary_type_with_all_ghosts_eij_, ghost_creation.each_boundary_type_contact_real_index_);
    SimpleDynamics<EulerianWCTimeStepInitialization> initialize_a_fluid_step(water_block);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(water_block, read_mesh_data.max_distance_between_nodes_);
    InteractionDynamics<WCEulerianViscousAccelerationInner> viscous_acceleration(water_block_inner);
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<Integration1stHalfAcousticRiemann> pressure_relaxation(water_block_inner, 200.0);
    InteractionWithUpdate<Integration2ndHalfAcousticRiemann> density_relaxation(water_block_inner, 200.0);
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionDynamics<ViscousForceFromFluidInFVM> viscous_force_on_solid(water_block_inner, ghost_creation.each_boundary_type_contact_real_index_);
    InteractionDynamics<AllForceAccelerationFromFluid> fluid_force_on_solid_update(water_block_inner, viscous_force_on_solid, ghost_creation.each_boundary_type_contact_real_index_);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_viscous_force_on_inserted_body(io_environment, viscous_force_on_solid, "TotalViscousForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>
        write_total_force_on_inserted_body(io_environment, fluid_force_on_solid_update, "TotalPressureForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>> write_maximum_speed(io_environment, water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with case specified initial condition if necessary.
    //----------------------------------------------------------------------
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 100.0;
    Real output_interval = 5.0; /**< time stamps for output. */
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
            initialize_a_fluid_step.exec();
            Real dt = get_fluid_time_step_size.exec();
            boundary_condition_setup.resetBoundaryConditions();
            viscous_acceleration.exec();
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);

            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
                     << GlobalStaticVariables::physical_time_
                     << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
        write_total_force_on_inserted_body.writeToFile(number_of_iterations);
        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
    write_total_viscous_force_on_inserted_body.testResult();

    return 0;
}
