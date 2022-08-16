/**
 * @file 	eulerian_taylor_green.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_taylor_green vortex flow example.
 * @author 	Chi Zhang, Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 50.0; /**< Global reference resolution. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 100;						/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
Real heat_capacity_ratio = 1.4;		/**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_body_shape;
		water_body_shape.push_back(Vecd(0.0, 0.0));
		water_body_shape.push_back(Vecd(0.0, DH));
		water_body_shape.push_back(Vecd(DL, DH));
		water_body_shape.push_back(Vecd(DL, 0.0));
		water_body_shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Application dependent initial condition
//----------------------------------------------------------------------
class TaylorGreenInitialCondition
	: public eulerian_compressible_fluid_dynamics::CompressibleFluidInitialCondition
{
public:
	explicit TaylorGreenInitialCondition(EulerianFluidBody &water)
		: eulerian_compressible_fluid_dynamics::CompressibleFluidInitialCondition(water){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial momentum and energy profile */
		rho_[index_i] = rho0_f;
		p_[index_i] = pow(c_f, 2) * rho_[index_i] / gamma_;
		vel_[index_i][0] = -cos(2.0 * Pi * pos_[index_i][0]) *
							 sin(2.0 * Pi * pos_[index_i][1]);
		vel_[index_i][1] = sin(2.0 * Pi * pos_[index_i][0]) *
							 cos(2.0 * Pi * pos_[index_i][1]);
		mom_[index_i] = rho_[index_i] * vel_[index_i];
		Real rho_e = p_[index_i] / (gamma_ - 1.0);
		E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].normSqr();
	}
};
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
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/** output environment. */
	InOutput in_output(sph_system);
	//handle command line arguments
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	EulerianFluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_body.defineParticlesAndMaterial<CompressibleFluidParticles, CompressibleFluid>(rho0_f, heat_capacity_ratio, mu_f);
	water_body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner water_body_inner(water_body);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initial condition with momentum and energy field */
	TaylorGreenInitialCondition initial_condition(water_body);
	/** Initialize a time step. */
	eulerian_compressible_fluid_dynamics::CompressibleFlowTimeStepInitialization time_step_initialization(water_body);
	/** Periodic BCs in x direction. */
	PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
	/** Periodic BCs in y direction. */
	PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
	/** Time step size with considering sound wave speed. */
	eulerian_compressible_fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_body);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	/** Here, we can use HLLC with Limiter Riemann solver for pressure relaxation and density and energy relaxation  */
	eulerian_compressible_fluid_dynamics::PressureRelaxationHLLCWithLimiterRiemannInner pressure_relaxation(water_body_inner);
	eulerian_compressible_fluid_dynamics::DensityAndEnergyRelaxationHLLCWithLimiterRiemannInner density_and_energy_relaxation(water_body_inner);
	/** Computing viscous acceleration. */
	eulerian_compressible_fluid_dynamics::ViscousAccelerationInner viscous_acceleration(water_body_inner);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	RegressionTestEnsembleAveraged<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
		write_total_mechanical_energy(in_output, water_body);
	/** Output the maximum speed of the fluid body. */
	RegressionTestEnsembleAveraged<BodyReducedQuantityRecording<MaximumSpeed>>
		write_maximum_speed(in_output, water_body);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	initial_condition.exec();
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_x.update_cell_linked_list_.parallel_exec();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_body.updateCellLinkedList();
		periodic_condition_x.update_cell_linked_list_.parallel_exec();
		periodic_condition_y.update_cell_linked_list_.parallel_exec();
		water_body_inner.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 5.0; /**< End time. */
	Real D_Time = 0.1;	 /**< Time stamps for output of body states. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
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
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force. */
			time_step_initialization.parallel_exec();
			Real dt = get_fluid_time_step_size.parallel_exec();
			viscous_acceleration.parallel_exec();
			/** Dynamics including pressure relaxation. */
			integration_time += dt;
			pressure_relaxation.parallel_exec(dt);
			density_and_energy_relaxation.parallel_exec(dt);
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
				{
					restart_io.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;
		}

		tick_count t2 = tick_count::now();
		write_total_mechanical_energy.writeToFile(number_of_iterations);
		write_maximum_speed.writeToFile(number_of_iterations);
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		 << " seconds." << endl;

	write_total_mechanical_energy.newResultTest();
	write_maximum_speed.newResultTest();

	return 0;
}
