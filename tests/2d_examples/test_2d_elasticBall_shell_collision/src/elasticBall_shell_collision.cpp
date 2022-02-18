/**
 * @file 	elasticBall_shell_collision.cpp
 * @brief 	an elastic ball bouncing within a confined shell boundary
 * @details This is a case to test elasticSolid -> shell impact/collision.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h"	//SPHinXsys Library.
using namespace SPH;	//Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0; 					/**< box length. */
Real DH = 4.0; 					/**< box height. */
Real resolution_ref = 0.025; 	/**< reference resolution. */
Real BW = resolution_ref * 1.; 	/**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center(2.0, 2.0);
Real ball_radius = 0.5;			
Real gravity_g = 1.0;
Real initial_ball_speed = 4.0;
Vec2d initial_velocity = initial_ball_speed * Vec2d(0.0, -1.);
//----------------------------------------------------------------------
//	Global paramters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45; 			
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class WallBoundary : public ThinStructure
{
public:
	WallBoundary(SPHSystem &sph_system, std::string body_name) 
	: ThinStructure(sph_system, body_name)
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
class BallBody : public SolidBody
{
public:
	BallBody(SPHSystem& system, std::string body_name) : SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(ball_center, ball_radius, 100, ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
/**
 * application dependent initial condition
 */
class BallInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	BallInitialCondition(SolidBody &body)
		: solid_dynamics::ElasticDynamicsInitialCondition(body) {};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		vel_n_[index_i] = initial_velocity;
	};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	In_Output 	in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	BallBody ball(sph_system, "BallBody");
	SharedPtr<ParticleGenerator> ball_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		ball_particle_generator = makeShared<ParticleGeneratorReload>(in_output, ball.getBodyName());
	SharedPtr<NeoHookeanSolid> ball_material = makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	ElasticSolidParticles ball_particles(ball, ball_material, ball_particle_generator);

	WallBoundary wall_boundary(sph_system, "Wall");
	SharedPtr<ParticleGenerator> wall_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		wall_particle_generator = makeShared<ParticleGeneratorReload>(in_output, wall_boundary.getBodyName());
	SharedPtr<LinearElasticSolid> wall_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	ShellParticles solid_particles(wall_boundary, wall_material, wall_particle_generator, BW);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner ball_inner(ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		RandomizePartilePosition  			ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_ball_state(in_output, sph_system.real_bodies_);
				ReloadParticleIO write_particle_reload_files(in_output, {&ball});

		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		write_ball_state.writeToFile(0.0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begines.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_ball_state.writeToFile(Real(ite) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0.0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner ball_inner(ball);
	SolidBodyRelationContact ball_contact(ball, {&wall_boundary});
	SolidBodyRelationContact wall_ball_contact(wall_boundary, {&ball});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	TimeStepInitialization 	ball_initialize_timestep(ball, gravity);
	solid_dynamics::CorrectConfiguration ball_corrected_configuration(ball_inner);
	solid_dynamics::AcousticTimeStepSize ball_get_time_step_size(ball);
	/** stress relaxation for the balls. */
	solid_dynamics::StressRelaxationFirstHalf ball_stress_relaxation_first_half(ball_inner);
	solid_dynamics::StressRelaxationSecondHalf ball_stress_relaxation_second_half(ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity ball_update_contact_density(ball_contact);
	solid_dynamics::ContactDensitySummation wall_ball_update_contact_density(wall_ball_contact);
	solid_dynamics::ContactForce ball_compute_solid_contact_forces(ball_contact);
	/** initial condition */
	BallInitialCondition ball_initial_velocity(ball);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp	body_states_recording(in_output, sph_system.real_bodies_);
	BodyStatesRecordingToVtp 	write_ball_state(in_output, { ball });
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	solid_particles.initializeNormalDirectionFromBodyShape();
	ball_corrected_configuration.parallel_exec();
	ball_initial_velocity.exec();
	
	/** Initial states output. */
	body_states_recording.writeToFile(0);
	/** Main loop. */
	int ite 		= 0;
	Real T0 		= 10.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.01*T0;
	Real Dt 		= 0.1*D_Time;			
	Real dt 		= 0.0; 	
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
				ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0) 
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				ball_update_contact_density.parallel_exec();
				wall_ball_update_contact_density.parallel_exec();
				ball_compute_solid_contact_forces.parallel_exec();
				ball_stress_relaxation_first_half.parallel_exec(dt);
				ball_stress_relaxation_second_half.parallel_exec(dt);

				ball.updateCellLinkedList();
				ball_contact.updateConfiguration();
				wall_ball_contact.updateConfiguration();

				ite++;
				Real dt_ball = ball_get_time_step_size.parallel_exec();
				dt = dt_ball;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}
		tick_count t2 = tick_count::now();
		write_ball_state.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
