/**
 * @file 	2d_FVM_double_mach_reflection.cpp
 * @brief 	This is the compressible test for the realization of FVM in the SPHinXsys.
 * @details We consider a double mach reflection case.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "2d_FVM_double_mach_reflection.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh ansys_mesh(double_mach_reflection_mesh1_fullpath);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, ansys_mesh.MinMeshEdge());
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody wave_block(sph_system, makeShared<WaveBody>("WaveBody"));
    wave_block.defineMaterial<CompressibleFluid>(rho0_another, heat_capacity_ratio);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    wave_block.generateParticlesWithReserve<BaseParticles, UnstructuredMesh>(ghost_boundary, ansys_mesh);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(wave_block, ansys_mesh);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfHLLCRiemann> pressure_relaxation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfHLLCRiemann> density_relaxation(water_block_inner);

    SimpleDynamics<DMFInitialCondition> initial_condition(wave_block);
    GhostCreationFromMesh ghost_creation(wave_block, ansys_mesh, ghost_boundary);
    DMFBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation);
    ReduceDynamics<CompressibleAcousticTimeStepSizeInFVM> get_fluid_time_step_size(wave_block, ansys_mesh.MinMeshEdge(), 0.2);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToMeshVtp write_real_body_states(wave_block, ansys_mesh);
    write_real_body_states.addToWrite<Real>(wave_block, "Density");
    write_real_body_states.addToWrite<Real>(wave_block, "Pressure");
    RegressionTestEnsembleAverage<ReducedQuantityRecording<MaximumSpeed>> write_maximum_speed(wave_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_block_inner.updateConfiguration();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 0.2;
    Real output_interval = 0.01; /**< time stamps for output. */
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
            }
            number_of_iterations++;
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

    if (sph_system.GenerateRegressionData())
    {
        write_maximum_speed.generateDataBase(1.0e-3, 1.0e-3);
    }
    else
    {
        write_maximum_speed.testResult();
    }
    return 0;
}
