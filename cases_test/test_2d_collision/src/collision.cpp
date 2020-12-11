/**
 * @file 	collision.cpp
 * @brief 	a soft ball with initial velocity bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Xiangyu Hu
 * @version 0.1
 * @version 0.3.0
 *			Here, I will try to implement a contact model for collision. 
 *			-- Chi ZHANG
 */
 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 8.0; 						/**< box length. */
Real DH = 4.0; 						/**< box height. */
Real particle_spacing_ref = 0.025; 	/**< reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< wall width for BCs. */
/**
 * @brief Material properties of the sphere.
 */
Real rho0_s = 1.0e3; 		
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45; 			
Vec2d ball_center_1(2.0, 2.0);
Vec2d ball_center_2(6.0, 2.0);
Real ball_radius = 0.5;			
Real initial_ball_speed = 1.0;
Real initial_direction = 0.25 * M_PI;
Vec2d initial_velocity = initial_ball_speed*Vec2d(cos(initial_direction), sin(initial_direction));
Real gravity_g = 1.0;
Real physical_viscosity = 10000.0; 
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name, int refinement_level)
		: SolidBody(sph_system, body_name, refinement_level)
	{
		/** Geometry definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(0.0, 0.0));
		inner_wall_shape.push_back(Point(0.0, DH));
		inner_wall_shape.push_back(Point(DL, DH));
		inner_wall_shape.push_back(Point(DL, 0.0));
		inner_wall_shape.push_back(Point(0.0, 0.0));

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/** Definition of the ball as a elastic structure. */
class FreeBall : public SolidBody
{
public:
	FreeBall(SPHSystem& system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/** Definition of the ball as a elastic structure. */
class DampingBall : public SolidBody
{
public:
	DampingBall(SPHSystem& system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addACircle(ball_center_2, ball_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/**
 * @brief Define ball material.
 */
class BallMaterial : public LinearElasticSolid
{
public:
	BallMaterial() : LinearElasticSolid()
	{
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;
		eta_0_ = physical_viscosity;

		assignDerivedMaterialParameters();
	}
};
/**
 * Setup material properties of myocardium
 */
class Material : public NeoHookeanSolid
{
public:
	Material() : NeoHookeanSolid()
	{
		rho_0_ 	= rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;
		eta_0_ = physical_viscosity;

		assignDerivedMaterialParameters();
	}
};
/**
 * application dependent initial condition
 */
class BallInitialCondition
	: public solid_dynamics::ElasticSolidDynamicsInitialCondition
{
public:
	BallInitialCondition(SolidBody* beam)
		: solid_dynamics::ElasticSolidDynamicsInitialCondition(beam) {};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		vel_n_[index_i] = initial_velocity;
	};
};
/** fluid observer body */
class FreeBallObserver : public FictitiousBody
{
public:
	FreeBallObserver(SPHSystem& system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(make_pair(ball_center_1, 0.0));
	}
};
/** fluid observer body */
class DampingBallObserver : public FictitiousBody
{
public:
	DampingBallObserver(SPHSystem& system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(make_pair(ball_center_2, 0.0));
	}
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref, 2);
	/** tag for run particle relaxation for the initial body fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** tag for computation start with relaxed body fitted particles distribution */
	sph_system.reload_particles_ = true;
	/** tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** output environment. */
	In_Output 	in_output(sph_system);
	/** handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary = new WallBoundary(sph_system, "Wall",	0);
	BallMaterial* wall_material = new BallMaterial();
	SolidParticles 					solid_particles(wall_boundary, wall_material);
	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	FreeBall* free_ball = new FreeBall(sph_system, "FreeBall", 0);
	Material* free_ball_material = new Material();
	ElasticSolidParticles 	free_ball_particles(free_ball, free_ball_material);
	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	DampingBall* damping_ball = new DampingBall(sph_system, "DampingBall", 0);
	Material* damping_ball_material = new Material();
	ElasticSolidParticles 	damping_ball_particles(damping_ball, damping_ball_material);
	/** Observer. */
	FreeBallObserver* free_ball_observer = new FreeBallObserver(sph_system, "FreeBallObserver", 0);
	BaseParticles 			free_ball_observer_particles(free_ball_observer);
	
	DampingBallObserver* damping_ball_observer = new DampingBallObserver(sph_system, "DampingBallObserver", 0);
	BaseParticles 			damping_ball_observer_particles(damping_ball_observer);
	/** Output the body states. */
	WriteBodyStatesToVtu 		write_body_states(in_output, sph_system.real_bodies_);
	/** check whether run particle relaxation for body fiited particle distribution. */
	if (sph_system.run_particle_relaxation_)
	{
		/** topology */
		SPHBodyInnerRelation* free_ball_inner = new SPHBodyInnerRelation(free_ball);
		SPHBodyInnerRelation* damping_ball_inner = new SPHBodyInnerRelation(damping_ball);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToPlt 		write_ball_state(in_output, { free_ball, damping_ball });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { free_ball, damping_ball});
		/** Random reset the relax solid particle position. */
		RandomizePartilePosition  			free_ball_random_particles(free_ball);
		RandomizePartilePosition  			damping_ball_random_particles(damping_ball);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		relax_dynamics::RelaxationStepInner damping_ball_relaxation_step_inner(damping_ball_inner);
		/**Particle relaxation starts here.*/
		free_ball_random_particles.parallel_exec(0.25);
		damping_ball_random_particles.parallel_exec(0.25);
		write_ball_state.WriteToFile(0.0);
		/**
		 * From here the time stepping begines.
		 * Set the starting time.
		 */
		int ite = 0;
		int relax_step = 1000;
		Real dt = 0.0;
		while (ite < relax_step)
		{
			free_ball_relaxation_step_inner.exec();
			damping_ball_relaxation_step_inner.exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_ball_state.WriteToFile(Real(ite) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.WriteToFile(0.0);
		return 0;
	}
	/** Algorithms. */
	SPHBodyInnerRelation*   free_ball_inner = new SPHBodyInnerRelation(free_ball);
	SolidBodyContactRelation* free_ball_contact = new SolidBodyContactRelation(free_ball, {wall_boundary});
	SPHBodyInnerRelation*   damping_ball_inner = new SPHBodyInnerRelation(damping_ball);
	SolidBodyContactRelation* damping_ball_contact = new SolidBodyContactRelation(damping_ball, {wall_boundary});
	SPHBodyContactRelation* free_ball_observer_contact = new SPHBodyContactRelation(free_ball_observer, { free_ball });
	SPHBodyContactRelation* damping_all_observer_contact = new SPHBodyContactRelation(damping_ball_observer, { damping_ball });
	/** Dynamics. */
	InitializeATimeStep 	free_ball_initialize_timestep(free_ball, &gravity);
	InitializeATimeStep 	damping_ball_initialize_timestep(damping_ball, &gravity);
	/** Kernel correction. */
	solid_dynamics::CorrectConfiguration free_ball_corrected_configuration(free_ball_inner);
	solid_dynamics::CorrectConfiguration damping_ball_corrected_configuration(damping_ball_inner);
	/** Time step size. */
	solid_dynamics::AcousticTimeStepSize free_ball_get_time_step_size(free_ball);
	solid_dynamics::AcousticTimeStepSize damping_ball_get_time_step_size(damping_ball);
	/** stress relaxation for the beam. */
	solid_dynamics::StressRelaxationFirstHalf free_ball_stress_relaxation_first_half(free_ball_inner);
	solid_dynamics::StressRelaxationSecondHalf free_ball_stress_relaxation_second_half(free_ball_inner);

	solid_dynamics::StressRelaxationFirstHalf damping_ball_stress_relaxation_first_half(damping_ball_inner);
	solid_dynamics::StressRelaxationSecondHalf damping_ball_stress_relaxation_second_half(damping_ball_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::SummationContactDensity free_ball_update_contact_density(free_ball_contact);
	solid_dynamics::ContactForce free_ball_compute_solid_contact_forces(free_ball_contact);

	solid_dynamics::SummationContactDensity damping_ball_update_contact_density(damping_ball_contact);
	solid_dynamics::ContactForce damping_ball_compute_solid_contact_forces(damping_ball_contact);
	/** Damping*/
	DampingBySplittingWithRandomChoice<SPHBodyInnerRelation, DampingBySplittingPairwise<Vec2d>, Vec2d>
		damping(damping_ball_inner, 0.5, damping_ball_particles.vel_n_, physical_viscosity);
	/** Observer and output. */
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_free_ball_displacement("Displacement", in_output, free_ball_observer_contact);
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_damping_ball_displacement("Displacement", in_output, damping_all_observer_contact);
	/** Now, pre-simulation. */
	if (sph_system.reload_particles_) 
	{
		unique_ptr<ReadReloadParticle>	
			free_ball_reload_particles(new ReadReloadParticle(in_output, {free_ball}, { "FreeBall"}));
		unique_ptr<ReadReloadParticle>	
			damping_ball_reload_particles(new ReadReloadParticle(in_output, { damping_ball}, { "DampingBall"}));
		free_ball_reload_particles->ReadFromFile();
		damping_ball_reload_particles->ReadFromFile();
	}
	GlobalStaticVariables::physical_time_ = 0.0;
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	free_ball_corrected_configuration.parallel_exec();
	damping_ball_corrected_configuration.parallel_exec();
	/** Initial states output. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_free_ball_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	write_damping_ball_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Main loop. */
	int ite 		= 0;
	Real T0 		= 10.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.01*T0;
	Real Dt 		= 0.1*D_Time;			
	Real dt 		= 0.0; 	
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** computation loop starts. */
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
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				free_ball_update_contact_density.parallel_exec();
				free_ball_compute_solid_contact_forces.parallel_exec();
				free_ball_stress_relaxation_first_half.parallel_exec(dt);
				free_ball_stress_relaxation_second_half.parallel_exec(dt);

				free_ball->updateCellLinkedList();
				free_ball_contact->updateConfiguration();

				damping_ball_update_contact_density.parallel_exec();
				damping_ball_compute_solid_contact_forces.parallel_exec();
				damping_ball_stress_relaxation_first_half.parallel_exec(dt);
				damping.parallel_exec(dt);
				damping_ball_stress_relaxation_second_half.parallel_exec(dt);

				damping_ball->updateCellLinkedList();
				damping_ball_contact->updateConfiguration();

				ite++;
				Real dt_free = free_ball_get_time_step_size.parallel_exec();
				Real dt_damping = damping_ball_get_time_step_size.parallel_exec();
				dt = SMIN(dt_free, dt_damping);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				write_free_ball_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
				write_damping_ball_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
			}
		}
		tick_count t2 = tick_count::now();
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
	return 0;
}
