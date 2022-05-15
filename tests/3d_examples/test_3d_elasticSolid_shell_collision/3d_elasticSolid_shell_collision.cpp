/**
 * @file 	3d_elasticSolid_shell_collision.cpp
 * @brief 	This is a benchmark test of the 3D elastic solid->shell contact/impact formulations.
 * @details  We consider the collisoin of an elastic ball bouncing in a spherical shell box.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;						  /**< plate length. */
Real DW = 4.0;						  /**< plate width. */
Real DH = 4.0;						  /**< plate height. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);	  /** Pseudo-normal. */
Real resolution_ref = 0.05;			  /**< reference resolution. */
Real thickness = resolution_ref * 1.; /**< plate thickness. */
size_t particle_number = 80;		  /** Initial reference particle spacing. */
int BWD = 1;						  /** Width of the boundary layer measured by number of particles. */
Real BW = resolution_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
BoundingBox system_domain_bounds(Vec3d(- thickness, - thickness, -DH / 2.0 - thickness),
								 Vec3d(DL  + thickness, DW + thickness, DH / 2.0 + thickness));
Real ball_radius = 0.5;
Real initial_ball_speed = 4.0;
Vec3d initial_velocity = initial_ball_speed * Vec3d(0.0, 0.0, -1.0);
Real gravity_g = 0.1;
//----------------------------------------------------------------------
//	Global paramters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries.
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit PlateParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			for (int j = 0; j < (particle_number + 2 * BWD); j++)
			{
				Real x = resolution_ref * i - BW + resolution_ref * 0.5;
				Real y = resolution_ref * j - BW + resolution_ref * 0.5;
				initializePositionAndVolume(Vecd(x, y, -DL / 2.0), resolution_ref * resolution_ref);
				initializeSurfaceProperties(n_0, thickness);
			}
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
	explicit BallInitialCondition(SolidBody &body)
		: solid_dynamics::ElasticDynamicsInitialCondition(body){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		vel_n_[index_i] = initial_velocity;
	};
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
	sph_system.run_particle_relaxation_ = true;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	InOutput in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	/** create a plate body. */
	SolidBody plate(sph_system, makeShared<DefaultShape>("Plate"));
	plate.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	plate.generateParticles<PlateParticleGenerator>();

	SolidBody ball(sph_system, makeShared<GeometricShapeSphere>(Vecd(DL/2.0, DW/2.0, 0.0), ball_radius, "BallBody"));
	ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
	{
		ball.generateParticles<ParticleGeneratorReload>(in_output, ball.getBodyName());
	}
	else
	{
		ball.defineBodyLevelSetShape()->writeLevelSet(ball);
		ball.generateParticles<ParticleGeneratorLattice>();
	}
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
		//	Define the methods for particle relaxation for ball.
		//----------------------------------------------------------------------
		RandomizePartilePosition ball_random_particles(ball);
		relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(in_output, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(in_output, {&ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		ball_random_particles.parallel_exec(0.25);
		write_relaxed_particles.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begines.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			ball_relaxation_step_inner.parallel_exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relaxed_particles.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner ball_inner(ball);
	SolidBodyRelationContact ball_contact(ball, {&plate});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	TimeStepInitialization ball_initialize_timestep(ball, gravity);
	solid_dynamics::CorrectConfiguration ball_corrected_configuration(ball_inner);
	solid_dynamics::AcousticTimeStepSize ball_get_time_step_size(ball);
	/** stress relaxation for the balls. */
	solid_dynamics::StressRelaxationFirstHalf ball_stress_relaxation_first_half(ball_inner);
	solid_dynamics::StressRelaxationSecondHalf ball_stress_relaxation_second_half(ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity ball_update_contact_density(ball_contact);
	solid_dynamics::ContactForce ball_compute_solid_contact_forces(ball_contact);
	/** initial condition */
	BallInitialCondition ball_initial_velocity(ball);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	BodyStatesRecordingToVtp write_ball_state(in_output, {ball});
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
	int ite = 0;
	Real T0 = 25.0;
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
				ball_initialize_timestep.parallel_exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				ball_update_contact_density.parallel_exec();
				ball_compute_solid_contact_forces.parallel_exec();
				ball_stress_relaxation_first_half.parallel_exec(dt);
				ball_stress_relaxation_second_half.parallel_exec(dt);

				ball.updateCellLinkedList();
				ball_contact.updateConfiguration();

				ite++;
				Real dt_free = ball_get_time_step_size.parallel_exec();
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
