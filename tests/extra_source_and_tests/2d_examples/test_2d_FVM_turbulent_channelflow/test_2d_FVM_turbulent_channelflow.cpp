/**
 * @file 	test_2d_FVM_turbulent_channelflow.h
 * @brief 	This is a test to show high reynolds number turbulent channel flow in FVM.
 * @author 	Yash Mandaokar, Feng Wang and Xiangyu Hu
 */

#include "test_2d_FVM_turbulent_channelflow.h"
#include "common_weakly_compressible_FVM_classes.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh read_mesh_data(mesh_file_path);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, read_mesh_data.MinMeshEdge());
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("TurbulentChannelFlow"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, UnstructuredMesh>(ghost_boundary, read_mesh_data);
    GhostCreationFromMesh ghost_creation(water_block, read_mesh_data, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(water_block, read_mesh_data);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfInnerRiemann> pressure_relaxation(water_block_inner, 1.0);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfInnerRiemann> density_relaxation(water_block_inner, 10000.0);
    SimpleDynamics<TurbulentChannelFlowInitialCondition> initial_condition(water_block);
    SimpleDynamics<fluid_dynamics::WallAdjacentCells> wall_adj_cell(water_block_inner, ghost_creation);

    InteractionWithUpdate<fluid_dynamics::KEpsilonStd1stHalfExtendedHLLCRiemannSolver> tke(water_block_inner, ghost_creation, 0.0);
    InteractionWithUpdate<fluid_dynamics::KEpsilonStd2ndHalfExtendedHLLCRiemannSolver> dissipationrate(water_block_inner, ghost_creation, 0.0);

    TurbulentChannelFlowBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(water_block, read_mesh_data.MinMeshEdge(), 1.0);
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceInner> turbulent_viscous_force(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::TkeGradientForceInner> tke_gradient_force(water_block_inner);

    // visualization in FVM with data in cell
    BodyStatesRecordingToMeshVtu write_real_body_states(water_block, read_mesh_data);
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(water_block);

    initial_condition.exec();
    wall_adj_cell.exec();
    water_block_inner.updateConfiguration();
    write_real_body_states.addToWrite<Real>(water_block, "Density");
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states.addToWrite<Vecd>(water_block, "Velocity");
    write_real_body_states.addToWrite<Real>(water_block, "TKE");
    write_real_body_states.addToWrite<Real>(water_block, "Dissipation");
    write_real_body_states.addToWrite<Real>(water_block, "TurblunetViscosity");
    write_real_body_states.addToWrite<Real>(water_block, "TKEProduction");
    write_real_body_states.addToWrite<Real>(water_block, "WallShearStress");
    write_real_body_states.addToWrite<Real>(water_block, "Ystar");
    write_real_body_states.addToWrite<Real>(water_block, "StrainRate");
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            boundary_condition_setup.resetBoundaryConditions();
            tke.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            dissipationrate.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            viscous_force.exec();
            turbulent_viscous_force.exec();
            tke_gradient_force.exec();
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);
            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
                     << physical_time
                     << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
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
