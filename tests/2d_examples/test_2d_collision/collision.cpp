/**
 * @file 	collision.cpp
 * @brief 	two soft balls with and without internal damping bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   //Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0;				  /**< box length. */
Real DH = 4.0;				  /**< box height. */
Real resolution_ref = 0.025;  /**< reference resolution. */
Real BW = resolution_ref * 4; /**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center_1(2.0, 2.0);
Vec2d ball_center_2(6.0, 2.0);
Real ball_radius = 0.5;
Real gravity_g = 1.0;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
Real physical_viscosity = 10000.0;
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name)
	{
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
class FreeBall : public SolidBody
{
public:
	FreeBall(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
class DampingBall : public SolidBody
{
public:
	DampingBall(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(ball_center_2, ball_radius, 100, ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
//----------------------------------------------------------------------
//	Observer particle generator.
//----------------------------------------------------------------------

class FreeBallObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	FreeBallObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(ball_center_1, 0.0));
	}
};
class DampingObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	DampingObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(ball_center_2, 0.0));
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FreeBall free_ball(sph_system, "FreeBall");
	SharedPtr<ParticleGenerator> free_ball_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		free_ball_particle_generator = makeShared<ParticleGeneratorReload>(in_output, free_ball.getBodyName());
	SharedPtr<NeoHookeanSolid> free_ball_material = makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	ElasticSolidParticles free_ball_particles(free_ball, free_ball_material, free_ball_particle_generator);

	DampingBall damping_ball(sph_system, "DampingBall");
	SharedPtr<ParticleGenerator> damping_ball_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		damping_ball_particle_generator = makeShared<ParticleGeneratorReload>(in_output, damping_ball.getBodyName());
	SharedPtr<NeoHookeanSolid> damping_ball_material = makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	ElasticSolidParticles damping_ball_particles(damping_ball, damping_ball_material, damping_ball_particle_generator);

	WallBoundary wall_boundary(sph_system, "Wall");
	SharedPtr<Solid> wall_material = makeShared<Solid>(rho0_s, free_ball_material->ContactStiffness());
	SolidParticles solid_particles(wall_boundary, wall_material);

	ObserverBody free_ball_observer(sph_system, "FreeBallObserver");
	ObserverParticles free_ball_observer_particles(free_ball_observer, makeShared<FreeBallObserverParticleGenerator>());

	ObserverBody damping_ball_observer(sph_system, "DampingBallObserver");
	ObserverParticles damping_ball_observer_particles(damping_ball_observer, makeShared<DampingObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner free_ball_inner(free_ball);
		BodyRelationInner damping_ball_inner(damping_ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		RandomizePartilePosition free_ball_random_particles(free_ball);
		RandomizePartilePosition damping_ball_random_particles(damping_ball);
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		relax_dynamics::RelaxationStepInner damping_ball_relaxation_step_inner(damping_ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_ball_state(in_output, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(in_output, {&free_ball, &damping_ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		free_ball_random_particles.parallel_exec(0.25);
		damping_ball_random_particles.parallel_exec(0.25);
		write_ball_state.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			free_ball_relaxation_step_inner.exec();
			damping_ball_relaxation_step_inner.exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_ball_state.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner free_ball_inner(free_ball);
	SolidBodyRelationContact free_ball_contact(free_ball, {&wall_boundary});
	BodyRelationInner damping_ball_inner(damping_ball);
	SolidBodyRelationContact damping_ball_contact(damping_ball, {&wall_boundary});
	BodyRelationContact free_ball_observer_contact(free_ball_observer, {&free_ball});
	BodyRelationContact damping_all_observer_contact(damping_ball_observer, {&damping_ball});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	TimeStepInitialization free_ball_initialize_timestep(free_ball, gravity);
	TimeStepInitialization damping_ball_initialize_timestep(damping_ball, gravity);
	solid_dynamics::CorrectConfiguration free_ball_corrected_configuration(free_ball_inner);
	solid_dynamics::CorrectConfiguration damping_ball_corrected_configuration(damping_ball_inner);
	solid_dynamics::AcousticTimeStepSize free_ball_get_time_step_size(free_ball);
	solid_dynamics::AcousticTimeStepSize damping_ball_get_time_step_size(damping_ball);
	/** stress relaxation for the balls. */
	solid_dynamics::StressRelaxationFirstHalf free_ball_stress_relaxation_first_half(free_ball_inner);
	solid_dynamics::StressRelaxationSecondHalf free_ball_stress_relaxation_second_half(free_ball_inner);
	solid_dynamics::StressRelaxationFirstHalf damping_ball_stress_relaxation_first_half(damping_ball_inner);
	solid_dynamics::StressRelaxationSecondHalf damping_ball_stress_relaxation_second_half(damping_ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation free_ball_update_contact_density(free_ball_contact);
	solid_dynamics::ContactForce free_ball_compute_solid_contact_forces(free_ball_contact);
	solid_dynamics::ContactDensitySummation damping_ball_update_contact_density(damping_ball_contact);
	solid_dynamics::ContactForce damping_ball_compute_solid_contact_forces(damping_ball_contact);
	/** Damping for one ball */
	DampingWithRandomChoice<DampingPairwiseInner<Vec2d>>
		damping(damping_ball_inner, 0.5, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		free_ball_displacement_recording("Position", in_output, free_ball_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		damping_ball_displacement_recording("Position", in_output, damping_all_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	free_ball_corrected_configuration.parallel_exec();
	damping_ball_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Initial states output.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile(0);
	free_ball_displacement_recording.writeToFile(0);
	damping_ball_displacement_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 10.0;
	Real End_Time = T0;
	Real D_Time = 0.01 * T0;
	Real Dt = 0.1 * D_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				free_ball_initialize_timestep.parallel_exec();
				damping_ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				free_ball_update_contact_density.parallel_exec();
				free_ball_compute_solid_contact_forces.parallel_exec();
				free_ball_stress_relaxation_first_half.parallel_exec(dt);
				free_ball_stress_relaxation_second_half.parallel_exec(dt);

				free_ball.updateCellLinkedList();
				free_ball_contact.updateConfiguration();

				damping_ball_update_contact_density.parallel_exec();
				damping_ball_compute_solid_contact_forces.parallel_exec();
				damping_ball_stress_relaxation_first_half.parallel_exec(dt);
				damping.parallel_exec(dt);
				damping_ball_stress_relaxation_second_half.parallel_exec(dt);

				damping_ball.updateCellLinkedList();
				damping_ball_contact.updateConfiguration();

				ite++;
				Real dt_free = free_ball_get_time_step_size.parallel_exec();
				Real dt_damping = damping_ball_get_time_step_size.parallel_exec();
				dt = SMIN(dt_free, dt_damping);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				free_ball_displacement_recording.writeToFile(ite);
				damping_ball_displacement_recording.writeToFile(ite);
			}
		}
		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	free_ball_displacement_recording.newResultTest();
	damping_ball_displacement_recording.newResultTest();

	return 0;
}
