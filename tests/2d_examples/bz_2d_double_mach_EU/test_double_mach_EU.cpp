/**
 * @file 	Lid_driven_cavity.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_Lid_driven_cavity example.
 * @author 	Zhentong Wang
 */
#include "sphinxsys.h" //	SPHinXsys Library.
#include "test_lid_driven_cavity_FVM.h"
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//read data from ANASYS mesh.file
	readMeshFile read_mesh_data(lid_driven_mesh_file_fullpath);
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
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorInFVM>(read_mesh_data.elements_center_coordinates_, read_mesh_data.elements_volumes_);
    water_block.addBodyStateForRecording<Real>("Density");
    /** Initial condition */
    GhostCreationFromMesh ghost_creation(water_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelationInFVM water_block_inner(water_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_2D_);
    water_block_inner.updateConfiguration();
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	FACBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation.each_boundary_type_with_all_ghosts_index_,
                                                       ghost_creation.each_boundary_type_with_all_ghosts_eij_, ghost_creation.each_boundary_type_contact_real_index_);
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(water_block, read_mesh_data.min_distance_between_nodes_);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_block_inner);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfAcousticRiemann> pressure_relaxation(water_block_inner);
	InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfAcousticRiemann> density_relaxation(water_block_inner);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingInMeshToVtp body_states_recording(io_environment, sph_system.real_bodies_,read_mesh_data.elements_nodes_connection_,read_mesh_data.point_coordinates_2D_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	boundary_condition_setup.resetBoundaryConditions();
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0; /**< End time. */
	Real D_Time = 1.0;	 /**< Time stamps for output of body states. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
    TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
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
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
        interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
