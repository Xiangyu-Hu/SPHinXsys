/**
 * @file 	shell_collision.cpp
 * @brief 	an elastic ball bouncing within a confined shell boundary
 * @details This is a case to test shell-shell modeling.
 * @details Both te ball and box are thin shell structures.
 * @author 	Massoud Rezavand and Xiangyu Hu
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
Real initial_ball_speed = 1.0;
Vec2d initial_velocity = initial_ball_speed*Vec2d(0.0, -1.);
//----------------------------------------------------------------------
//	Global paramters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0;				   /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;				   /** Poisson ratio. */
Real physical_viscosity = 200.0;   /** physical damping, here we choose the same value as numerical viscosity. */
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
/** Define application dependent particle generator for thin structure. */
class BallParticleGenerator : public ParticleGeneratorDirect
{
public:
	BallParticleGenerator() : ParticleGeneratorDirect()
	{
		for (Real t = 0; t <= 4*Pi*ball_radius; t += resolution_ref / 2.) 
		{
            Real x = ball_radius * cos(t) + ball_center[0];
            Real y = ball_radius * sin(t) + ball_center[1];
			positions_volumes_.push_back(std::make_pair(Vecd(x, y), resolution_ref));
		}
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
	ThinStructure free_ball(sph_system, "FreeBall", makeShared<SPHAdaptation>(1.15, 1.0));
	ShellParticles free_ball_particles(free_ball,
										makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
										makeShared<BallParticleGenerator>(), BW);

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
		BodyRelationInner free_ball_inner(free_ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		RandomizePartilePosition  			free_ball_random_particles(free_ball);
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_ball_state(in_output, sph_system.real_bodies_);
				ReloadParticleIO write_particle_reload_files(in_output, {&free_ball});

		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		free_ball_random_particles.parallel_exec(0.25);
		write_ball_state.writeToFile(0.0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begines.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			free_ball_relaxation_step_inner.exec();
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
	BodyRelationInner free_ball_inner(free_ball);
	SolidBodyRelationContact free_ball_contact(free_ball, {&wall_boundary});
	SolidBodyRelationContact wall_ball_contact(wall_boundary, {&free_ball});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	TimeStepInitialization 	free_ball_initialize_timestep(free_ball, gravity);
	thin_structure_dynamics::ShellCorrectConfiguration free_ball_corrected_configuration(free_ball_inner);
	thin_structure_dynamics::ShellAcousticTimeStepSize free_ball_get_time_step_size(free_ball);
	/** stress relaxation for the balls. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf free_ball_stress_relaxation_first_half(free_ball_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf free_ball_stress_relaxation_second_half(free_ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity free_ball_update_contact_density(free_ball_contact);
	solid_dynamics::ShellContactDensity wall_ball_update_contact_density(wall_ball_contact);
	solid_dynamics::ContactForce free_ball_compute_solid_contact_forces(free_ball_contact);
	solid_dynamics::ContactForce wall_compute_solid_contact_forces(wall_ball_contact);
	/** initial condition */
	BallInitialCondition ball_initial_velocity(free_ball);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp	body_states_recording(in_output, sph_system.real_bodies_);
	BodyStatesRecordingToVtp 	write_ball_state(in_output, { free_ball });
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	solid_particles.initializeNormalDirectionFromBodyShape();
	free_ball_corrected_configuration.parallel_exec();
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
				free_ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0) 
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				free_ball_update_contact_density.parallel_exec();
				wall_ball_update_contact_density.parallel_exec();
				free_ball_compute_solid_contact_forces.parallel_exec();
				wall_compute_solid_contact_forces.parallel_exec();
				free_ball_stress_relaxation_first_half.parallel_exec(dt);
				free_ball_stress_relaxation_second_half.parallel_exec(dt);

				free_ball.updateCellLinkedList();
				free_ball_contact.updateConfiguration();
				wall_boundary.updateCellLinkedList();
				wall_ball_contact.updateConfiguration();

				ite++;
				Real dt_free = free_ball_get_time_step_size.parallel_exec();
				dt = dt_free;
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
