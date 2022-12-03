/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;						  /**< Reference length. */
Real DH = 3.0;						  /**< Reference and the height of main channel. */
Real DL1 = 0.7 * DL;				  /**< The length of the main channel. */
Real resolution_ref = 0.15;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(1.0, DH / (2.0 * (DL - DL1)));
Real Re = 100.0;					/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the water block in T shape polygon. */
std::vector<Vecd> water_block_shape{
	Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH),
	Vecd(DL, 2.0 * DH), Vecd(DL, -DH), Vecd(DL1, -DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0)};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
	Vecd(-DL_sponge - BW, -BW), Vecd(-DL_sponge - BW, DH + BW), Vecd(DL1 - BW, DH + BW), Vecd(DL1 - BW, 2.0 * DH + BW),
	Vecd(DL + BW, 2.0 * DH + BW), Vecd(DL + BW, -DH - BW), Vecd(DL1 - BW, -DH - BW), Vecd(DL1 - BW, -BW), Vecd(-DL_sponge - BW, -BW)};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
	Vecd(-DL_sponge - BW, 0.0), Vecd(-DL_sponge - BW, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH + BW),
	Vecd(DL, 2.0 * DH + BW), Vecd(DL, -DH - BW), Vecd(DL1, -DH - BW), Vecd(DL1, 0.0), Vecd(-DL_sponge - BW, 0.0)};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};

class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Define emitter buffer inflow boundary condition
//----------------------------------------------------------------------
class EmitterBufferInflowCondition : public fluid_dynamics::InflowVelocityCondition
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowCondition(BodyAlignedBoxByCell &aligned_box_part)
		: InflowVelocityCondition(aligned_box_part),
		  u_ave_(0), u_ref_(U_f), t_ref_(4.0) {}

	// here every argument parameters and return value are in frame (local) coordinate
	Vecd getPrescribedVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];

		if (position[0] < halfsize_[0])
		{
			u = 1.5 * u_ave_ * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
			v = 0.0;
		}
		return Vec2d(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -DH - BW), Vec2d(DL + BW, 2.0 * DH + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex_relation(water_block_inner, {&wall_boundary});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex_relation);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex_relation);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_relation);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
		inlet_outlet_surface_particle_indicator(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_complex_relation);
	water_block.addBodyStateForRecording<Real>("Pressure");		   // output for debug
	water_block.addBodyStateForRecording<int>("SurfaceIndicator"); // output for debug

	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

	Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
	Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
	BodyAlignedBoxByParticle emitter(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(emitter_translation)), emitter_halfsize));
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection, BodyAlignedBoxByParticle> emitter_inflow_injection(emitter, 10, 0);

	Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
	Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;
	BodyAlignedBoxByCell emitter_buffer(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
	SimpleDynamics<EmitterBufferInflowCondition, BodyAlignedBoxByCell> emitter_buffer_inflow_condition(emitter_buffer);

	Vec2d disposer_up_halfsize = Vec2d(0.3 * DH, 0.5 * BW);
	Vec2d disposer_up_translation = Vec2d(DL + 0.05 * DH, 2.0 * DH) - disposer_up_halfsize;
	BodyAlignedBoxByCell disposer_up(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(disposer_up_translation)), disposer_up_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion, BodyAlignedBoxByCell> disposer_up_outflow_deletion(disposer_up, yAxis);

	Vec2d disposer_down_halfsize = disposer_up_halfsize;
	Vec2d disposer_down_translation = Vec2d(DL1 - 0.05 * DH, - DH) + disposer_down_halfsize;
	BodyAlignedBoxByCell disposer_down(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Rotation2d(Pi), Vec2d(disposer_down_translation)), disposer_down_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion, BodyAlignedBoxByCell> disposer_down_outflow_deletion(disposer_down, yAxis);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);
	RestartIO restart_io(io_environment, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 100.0;
	Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;							 /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_body_states.writeToFile();
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			inlet_outlet_surface_particle_indicator.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				emitter_buffer_inflow_condition.exec();
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** inflow injection*/
			emitter_inflow_injection.exec();
			disposer_up_outflow_deletion.parallel_exec();
			disposer_down_outflow_deletion.parallel_exec();

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_complex_relation.updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;

	return 0;
}
