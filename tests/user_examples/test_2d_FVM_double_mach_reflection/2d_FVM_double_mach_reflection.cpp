/**
 * @file 	2d_FVM_double_mach_reflection.cpp
 * @brief 	This is the compressible test for the realizaiton of FVM in the SPHinXsys.
 * @details We consider a double mach reflection case.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "2d_FVM_double_mach_reflection.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANASYS mesh.file
    readMeshFile read_mesh_data(double_mach_reflection_mesh_fullpath);
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
    EulerianFluidBody wave_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    wave_block.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_another, heat_capacity_ratio);
    wave_block.generateParticles<ParticleGeneratorInFVM>(read_mesh_data.elements_center_coordinates_, read_mesh_data.elements_volumes_);
    wave_block.addBodyStateForRecording<Real>("Density");
    /** Initial condition and register variables*/
    SimpleDynamics<DMFInitialCondition> initial_condition(wave_block);
    GhostCreationFromMesh ghost_creation(wave_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(wave_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
    water_block_inner.updateConfiguration();
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Boundary conditions set up */
    DMFBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation.each_boundary_type_with_all_ghosts_index_,
                                                       ghost_creation.each_boundary_type_with_all_ghosts_eij_, ghost_creation.each_boundary_type_contact_real_index_);
    SimpleDynamics<EulerianCompressibleTimeStepInitialization> initialize_a_fluid_step(wave_block);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<CompressibleAcousticTimeStepSizeInFVM> get_fluid_time_step_size(wave_block, read_mesh_data.max_distance_between_nodes_, 0.08);
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<Integration1stHalfHLLCRiemann> pressure_relaxation(water_block_inner);
    InteractionWithUpdate<Integration2ndHalfHLLCRiemann> density_relaxation(water_block_inner);
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with case specified initial condition if necessary.
    //----------------------------------------------------------------------
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real dt = get_fluid_time_step_size.exec();
            boundary_condition_setup.resetBoundaryConditions();
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
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

    return 0;
}
