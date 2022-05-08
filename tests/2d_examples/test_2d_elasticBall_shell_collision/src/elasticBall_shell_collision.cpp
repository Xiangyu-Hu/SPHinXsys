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
Real thickness = resolution_ref * 1.; 	/**< wall width for BCs. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
BoundingBox system_domain_bounds(Vec2d(-thickness, -thickness), Vec2d(DL + thickness, DH + thickness));
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
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-thickness, -thickness));
		outer_wall_shape.push_back(Vecd(-thickness, DH + thickness));
		outer_wall_shape.push_back(Vecd(DL + thickness, DH + thickness));
		outer_wall_shape.push_back(Vecd(DL + thickness, -thickness));
		outer_wall_shape.push_back(Vecd(-thickness, -thickness));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
class BallBody : public MultiPolygonShape
{
public:
	explicit BallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(ball_center, ball_radius, 100, ShapeBooleanOps::add);
	}
};
/**
 * application dependent initial condition
 */
class BallInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit BallInitialCondition(SolidBody &body)
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
	sph_system.run_particle_relaxation_ = true;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	InOutput 	in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody ball(sph_system, makeShared<BallBody>("BallBody"));
	ball.defineBodyLevelSetShape();
	ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? ball.generateParticles<ParticleGeneratorReload>(in_output, ball.getBodyName())
		: ball.generateParticles<ParticleGeneratorLattice>();

	// Note the wall boundary here has sharp corner, and is a numerical invalid elastic shell structure,
	// and its dynamics is not able to be modeled by the shell dynamics in SPHinXsys in the current version.
	// Here, we use it simply as a rigid shell.
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineAdaptation<SPHAdaptation>(1.15, 1.0);
	wall_boundary.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(wall_boundary);
	//here dummy linear elastic solid is use because no solid dynamics in particle relaxation
	wall_boundary.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(1.0, 1.0, 0.0);
	wall_boundary.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness);
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner ball_inner(ball);
		BodyRelationInner wall_boundary_inner(wall_boundary);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation for ball.
		//----------------------------------------------------------------------
		RandomizePartilePosition  			ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation for wall boundary.
		//----------------------------------------------------------------------
		RandomizePartilePosition  			wall_boundary_random_particles(wall_boundary);
		relax_dynamics::ShellRelaxationStepInner
		relaxation_step_wall_boundary_inner(wall_boundary_inner, thickness, level_set_refinement_ratio);
		relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(wall_boundary_inner, thickness, cos(Pi / 3.75));
		wall_boundary.addBodyStateForRecording<int>("UpdatedIndicator");
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(in_output, sph_system.real_bodies_);
		MeshRecordingToPlt write_mesh_cell_linked_list(in_output, wall_boundary, wall_boundary.cell_linked_list_);
		ReloadParticleIO write_particle_reload_files(in_output, {&ball, &wall_boundary});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		wall_boundary_random_particles.parallel_exec(0.25);

		relaxation_step_wall_boundary_inner.mid_surface_bounding_.parallel_exec();
		write_relaxed_particles.writeToFile(0);
		wall_boundary.updateCellLinkedList();
		write_mesh_cell_linked_list.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begines.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.parallel_exec();
			relaxation_step_wall_boundary_inner.parallel_exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relaxed_particles.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		shell_normal_prediction.exec();
		write_relaxed_particles.writeToFile(ite);
		write_particle_reload_files.writeToFile(0);
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
	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
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
