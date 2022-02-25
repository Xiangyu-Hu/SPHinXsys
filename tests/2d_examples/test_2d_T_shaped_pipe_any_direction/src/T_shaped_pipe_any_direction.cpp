/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file T_shaped_pipe.h
 * @brief 	This is the benchmark test of inlet and outlet condition in any direction.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author	Huiqiang Yue and Anyong Zhang
 */

#include "sphinxsys.h"
#include "utilities.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;						  /**< Reference length. */
Real DH = 3.0;						  /**< Reference height. */
Real DL1 = 0.7 * DL;				  /**< The length of the main channel. */
Real resolution_ref = 0.15;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */

// angle to rotate around the center
Real rotate_angle = 135.0 / 180.0 * Pi;
Vecd rotate_center((DL1 - DL_sponge) * 0.5, DH * 0.5);

/* system domain bounds */
BoundingBox new_bounding_box = calculateNewBoundingBox({Vec2d(-DL_sponge, -DH), Vec2d(DL + BW, 2.0 * DH)}, rotate_angle, rotate_center);
BoundingBox system_domain_bounds(new_bounding_box);
/** Prescribed fluid body domain bounds*/
BoundingBox fluid_body_domain_bounds(new_bounding_box);
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
/** the water block in T shape polygen. */
std::vector<Vecd> water_block_shape{
	rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge, DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, 2.0 * DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL, 2.0 * DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL, -DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, -DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, 0.0), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center)};
/** the outer wall polygen. */
std::vector<Vecd> outer_wall_shape{
	rotatePointAroundPoint(Vecd(-DL_sponge, -BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge, DH + BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1 - BW, DH + BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1 - BW, 2.0 * DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL + BW, 2.0 * DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL + BW, -DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1 - BW, -DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1 - BW, -BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge, -BW), rotate_angle, rotate_center)};
/** the inner wall polygen. */
std::vector<Vecd> inner_wall_shape{
	rotatePointAroundPoint(Vecd(-DL_sponge - BW, 0.0), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge - BW, DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, 2.0 * DH + BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL, 2.0 * DH + BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL, -DH - BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, -DH - BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1, 0.0), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge - BW, 0.0), rotate_angle, rotate_center)};
//-------------------------------------------------

// define boundary faces
// points describing the border of face(line in 2D)
StdVec<Vecd> outlet_points1{
	rotatePointAroundPoint(Vecd(DL1 - BW, 2.0 * DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL + BW, 2.0 * DH), rotate_angle, rotate_center)};
// direction pointing from face to computational domain
Vecd outlet_direction = rotatePointAroundPoint(Vecd(0, -1), rotate_angle);
SegmentFace outflow_face1(outlet_points1, outlet_direction);

StdVec<Vecd> outlet_points2{
	rotatePointAroundPoint(Vecd(DL + BW, -DH), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(DL1 - BW, -DH), rotate_angle, rotate_center),
};
SegmentFace outflow_face2(outlet_points2, rotatePointAroundPoint(Vecd(0, 1), rotate_angle));

StdVec<Vecd> inlet_points{
	rotatePointAroundPoint(Vecd(-DL_sponge, -BW), rotate_angle, rotate_center),
	rotatePointAroundPoint(Vecd(-DL_sponge, DH + BW), rotate_angle, rotate_center)};
SegmentFace inflow_face(inlet_points, rotatePointAroundPoint(Vecd(1, 0), rotate_angle));

//----------------------------------------------------------------------
//	Define case dependent SPH bodies.
//----------------------------------------------------------------------
/** Water block body definition. */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/** Wall boundary body definition. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Define case dependent SPH body parts.
//----------------------------------------------------------------------
/** create the emitter shape. */
MultiPolygon creatEmitterShape()
{
	std::vector<Vecd> emmiter_shape{
		rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge, DH), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge + BW, DH), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge + BW, 0.0), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center)};

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(emmiter_shape, ShapeBooleanOps::add);
	return multi_polygon;
}

/** create the emitter buffer shape . */
MultiPolygon createEmitterBufferShape()
{
	std::vector<Vecd> emitter_buffer_shape{
		rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge, DH), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(0.0, DH), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(0.0, 0.0), rotate_angle, rotate_center),
		rotatePointAroundPoint(Vecd(-DL_sponge, 0.0), rotate_angle, rotate_center)};

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(emitter_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Define emitter buffer inflow boundary condition
//----------------------------------------------------------------------
class EmitterBufferInflowConditionWithFace : public VelocityInflowConditionWithFace
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowConditionWithFace(FluidBody &body, BodyRegionByCellsWithFace &body_part)
		: VelocityInflowConditionWithFace(body, body_part),
		  u_ave_(0), u_ref_(U_f), t_ref_(4.0) {}

	Vecd defineVelocityProfile(Vecd &position, Vecd &velocity) override
	{
		Vecd old_pos = rotatePointAroundPoint(position, -rotate_angle, rotate_center);
		Real vv = u_ave_ * 6.0 * old_pos[1] * (DH - old_pos[1]) / DH / DH;
		Real u = vv * std::cos(rotate_angle);
		Real v = vv * std::sin(rotate_angle);
		return Vecd(u, v);
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
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	WaterBlock water_block(system, "WaterBody");
	water_block.setBodyDomainBounds(fluid_body_domain_bounds);
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
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
	/** Emmiter. */
	BodyRegionByParticleWithFace emitter(water_block, inflow_face, 4);
	InflowInjectingWithFace emitter_inflow_injecting(water_block, emitter, 10);
	/** Emitter condition. */
	BodyRegionByCellsWithFace emitter_buffer(water_block, inflow_face, 8);
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
	/** recycle real fluid particle to buffer particles at outlet. */
	BodyRegionByCellsWithFace face1(water_block, outflow_face1);
	BodyRegionByCellsWithFace face2(water_block, outflow_face2);

	DoNothingConditionWithFace outflow1(water_block, face1);
	DoNothingConditionWithFace outflow2(water_block, face2);

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
	Real D_Time = End_Time / 200.0; /**< Time stamps for output of body states. */
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

			/** inflow injecting*/
			emitter_inflow_injecting.exec();

			/*check outflow boundary*/
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
