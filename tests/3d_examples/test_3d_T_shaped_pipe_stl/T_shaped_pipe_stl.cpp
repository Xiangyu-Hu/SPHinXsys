/* ---------------------------------------------------------------------------*
 *                       SPHinXsys: 3D T-shaped-pipe example                        *
 * ----------------------------------------------------------------------------*
 * This is the test case of 3D T shaped	pipe by boundary face		  				  *
 * ---------------------------------------------------------------------------*/
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"

using namespace SPH;

// for geometry
Real resolution_ref = 0.025; // particle spacing

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-0.14459, -1.1172, -0.32998), Vecd(1.5073, 0.62912, 0.86797));
// for material properties of the fluid
Real rho0_f = 1.0;
Real U_f = 1.0;
Real L_f = 0.5; // reference length
Real c_f = 10.0 * U_f;
Real Re = 100.0;
Real mu_f = rho0_f * U_f * L_f / Re;

// boundary faces
Vecd inlet_center(0);	// center of face
Vecd inlet_direction(0.91068, -0.24402, 0.33333);	// direction pointing from face to computational domain
// points describing the border of face
StdVec<Vecd> inlet_bound_points{
	Vecd(-0.11956391412724343, 0.03929345153227706, 0.3554194966507927),
	Vecd(0.12122779105891612, -0.04292526559422116, -0.362623966253959),
	Vecd(-0.10867079399866064, -0.3718001341002144, 0.02471754554198849),
	Vecd(0.10667774485754022, 0.36451415483264793, -0.024606140397148306)};
SegmentFace inlet_face(inlet_bound_points, inlet_direction, inlet_center);

Vecd outlet_center_1(1.1774, 0.48453, 0.13812);
Vecd outlet_direction_1(-0.33333, -0.91068, 0.24402); // direction pointing from face to computational domain
StdVec<Vecd> outlet_bound_points_1{
	Vecd(1.177682940165783, 0.5829555094646648, 0.5059034703281549),
	Vecd(0.8236588125639355, 0.6100340957001089, 0.12335607178208308),
	Vecd(1.1676159033619389, 0.38927654392002875, -0.23066805581976407),
	Vecd(1.5400962651041668, 0.3464510438660897, 0.11832255338016111)};
SegmentFace outlet_face_1(outlet_bound_points_1, outlet_direction_1, outlet_center_1);

Vecd outlet_center_2(0.64402, -0.97256, 0.52855);
Vecd outlet_direction_2(0.33333, 0.91068, -0.24402);
StdVec<Vecd> outlet_bound_points_2{
	Vecd(0.6476893177675176, -0.8778725585931972, 0.8869559411054205),
	Vecd(0.28561573376649535, -0.847517281860303, 0.5056417318676291),
	Vecd(0.6774248294970703, -1.0804477908307577, 0.17155451184736226),
	Vecd(0.9800273900389873, -1.0820049567715868, 0.5791059373171119)};
SegmentFace outlet_face_2(outlet_bound_points_2, outlet_direction_2, outlet_center_2);

//	import the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		std::string file_name = "./input/T_shaped_tube_fluid.stl"; // in milimeter
		TriangleMeshShapeSTL stl_shape(file_name, Vecd(0), 0.001); // convert unit from milimeter to meter
		body_shape_.add<LevelSetShape>(this, stl_shape, true);
	}
};

//	import the static solid wall boundary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		std::string file_name = "./input/T_shaped_tube_solid.stl";
		TriangleMeshShapeSTL stlshape(file_name, Vecd(0), 0.001);
		body_shape_.add<LevelSetShape>(this, stlshape, true);
	}
};

/**
 * @brief define velocity buffer
 */
class EmitterBufferInflowConditionWithFace : public VelocityInflowConditionWithFace
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowConditionWithFace(FluidBody &body, BodyRegionByCellsWithFace &body_part)
		: VelocityInflowConditionWithFace(body, body_part),
		  u_ave_(0), u_ref_(U_f), t_ref_(4.0) {}

	Vecd defineVelocityProfile(Vecd &position, Vecd &velocity) override
	{
		Vecd vel = u_ave_ * inlet_direction.normalize();
		return vel;
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}
};

// the main program
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0;
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock water_block(system, "WaterBody");
	// water_block.setBodyDomainBounds(fluid_body_domain_bounds);
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	//---------------------------------------------
	// Relaxing solid wall particles.
	//--------------------------------------------
	/* body topology only for particle relation */
	BodyRelationInner wall_bound_model_inner(wall_boundary);
	/* Random reset the solid body particle position */
	RandomizePartilePosition random_wall_particles(wall_boundary);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_bound_model_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_wall_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	//----------------------------------------------------------------------
	//	Relax particles of the solid body.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		relaxation_step_inner.parallel_exec();
		ite_p++;
		if (ite_p % 200 == 0)
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
	}
	cout << "The physics relaxation process of wall_boundary body finish ! \n";
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner water_block_inner(water_block);
	ComplexBodyRelation water_block_complex_relation(water_block_inner, {&wall_boundary});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block);

	/* Define particle emitter at inlet */
	BodyRegionByParticleWithFace emitter(water_block, inlet_face, 2);
	InflowInjectingWithFace emitter_inflow_injecting(water_block, emitter, 10);

	/** Define inflow condition. */
	BodyRegionByCellsWithFace emitter_buffer(water_block, inlet_face, 2);
	EmitterBufferInflowConditionWithFace emitter_buffer_inflow_condition(water_block, emitter_buffer);

	/** time-space method to detect surface particles. */
	fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex
		inlet_outlet_surface_particle_indicator(water_block_complex_relation);
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(water_block_complex_relation);
	/** We can output a method-specific particle data for debug */
	fluid_particles.addAVariableToWrite<Real>("Pressure");
	fluid_particles.addAVariableToWrite<int>("SurfaceIndicator");
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex_relation);
	/** Density relaxation. */
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex_relation);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex_relation);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex_relation);
	/** Define outlet face to */
	BodyRegionByCellsWithFace face1(water_block, outlet_face_1);
	BodyRegionByCellsWithFace face2(water_block, outlet_face_2);
	/** recycle real fluid particle to buffer particles at outlet. */
	DoNothingConditionWithFace outflow1(water_block, face1);
	DoNothingConditionWithFace outflow2(water_block, face2);

	/**/
	// ModifiedDoNothingConditionWithFace outflow1(water_block, face1);
	// ModifiedDoNothingConditionWithFace outflow2(water_block, face2);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	inlet_outlet_surface_particle_indicator.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex_relation.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 100.0;			/**< End time. */
	Real D_Time = End_Time / 400.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;					/**< Default acoustic time step sizes. */
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
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			int iii = 0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				emitter_buffer_inflow_condition.parallel_exec();
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				++iii;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "  Dt/dt= " << iii << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** inflow injecting */
			emitter_inflow_injecting.exec();

			/* check particle at outlet */
			outflow1.parallel_exec();
			outflow2.parallel_exec();

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedList();
			water_block_complex_relation.updateConfiguration();
			inlet_outlet_surface_particle_indicator.parallel_exec();
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