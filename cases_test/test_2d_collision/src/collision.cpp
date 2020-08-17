/**
 * @file 	collision.cpp
 * @brief 	a soft ball with initial velocity bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Xiangyu Hu
 * @version 0.1
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
Real DL = 5.366; 						/**< box length. */
Real DH = 5.366; 						/**< box height. */
Real particle_spacing_ref = 0.05; 		/**< reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< wall width for BCs. */
/**
 * @brief Material properties of the sphere.
 */
Real rho0_s = 1.0e3; 			//reference density
Real Youngs_modulus = 2.0e6;	//reference Youngs modulus
Real poisson = 0.3975; 			//Poisson ratio
Vec2d ball_center(2.0, 2.0);		/**< Location of the cylinder center. */
Real ball_radius = 0.5;			/**< Radius of the cylinder. */
Real initial_ball_speed = 1.0;
Real initial_direction = 0.25 * M_PI;
Vec2d initial_velocity = initial_ball_speed*Vec2d(cos(initial_direction), sin(initial_direction));
Real gravity_g = 1.0;					/**< gravity */

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(sph_system, body_name, refinement_level, op)
	{
		/** Geometry definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));
		body_shape_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(0.0, 0.0));
		inner_wall_shape.push_back(Point(0.0, DH));
		inner_wall_shape.push_back(Point(DL, DH));
		inner_wall_shape.push_back(Point(DL, 0.0));
		inner_wall_shape.push_back(Point(0.0, 0.0));
		body_shape_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/** Definition of the ball as a elastic structure. */
class Ball : public SolidBody
{
public:
	Ball(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		body_shape_.addACircle(ball_center, ball_radius, 100, ShapeBooleanOps::add);
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
	void Update(size_t index_particle_i, Real dt) override {
		/** initial velocity profile */
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		base_particle_data_i.vel_n_ = initial_velocity;
	};
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** tag for run particle relaxation for the initial body fitted distribution */
	sph_system.run_particle_relaxation_ = true;
	/** tag for computation start with relaxed body fitted particles distribution */
	sph_system.reload_particles_ = false;
	/** tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** output environment. */
	In_Output 	in_output(sph_system);
	//handle command line arguments
	sph_system.handleCommandlineOptions(ac, av);

	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary 
		= new WallBoundary(sph_system, "Wall",	0, ParticlesGeneratorOps::lattice);
	SolidParticles 					solid_particles(wall_boundary);
	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	Ball* ball = new Ball(sph_system, "Ball", 0, ParticlesGeneratorOps::lattice);
	BallMaterial* ball_material = new BallMaterial();
	ElasticSolidParticles 	inserted_body_particles(ball, ball_material);

	/** Output the body states. */
	WriteBodyStatesToVtu 		write_body_states(in_output, sph_system.real_bodies_);
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	/** check whether run particle relaxation for body fiited particle distribution. */
	if (sph_system.run_particle_relaxation_)
	{
		/** add background level set for particle realxation. */
		ball->addLevelsetMesh();

		/** topology */
		SPHBodyInnerRelation* ball_inner_relation = new SPHBodyInnerRelation(ball);

		/** Write backgroung level set. */
		WriteBodyMeshToPlt write_ball_background_mesh(in_output, ball);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToVtu 		write_ball_to_vtu(in_output, { ball });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { ball });

		/** Random reset the relax solid particle position. */
		RandomizePartilePosition  			random_particles(ball);
		/** bounding particles to insert body surface. */
		relax_dynamics::BodySurfaceBounding
			body_surface_bounding(ball, new NearBodySurface(ball));
		/** Compute the time step for particle relaxation. */
		relax_dynamics::GetTimeStepSize get_relax_timestep(ball);
		/** Physics relaxation algorithm. */
		relax_dynamics::PhysicsRelaxationInner 	relax_process(ball_inner_relation);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		write_ball_background_mesh.WriteToFile(0.0);
		random_particles.parallel_exec(0.25);
		body_surface_bounding.parallel_exec();
		write_ball_to_vtu.WriteToFile(0.0);

		/**
		 * From here the time stepping begines.
		 * Set the starting time.
		 */
		int ite = 0;
		int relax_step = 1000;
		int diffusion_step = 100;
		Real dt = 0.0;
		while (ite < relax_step)
		{

			ball->UpdateCellLinkedList();
			ball_inner_relation->updateConfiguration();

			relax_process.parallel_exec(dt);
			body_surface_bounding.parallel_exec();
			dt = get_relax_timestep.parallel_exec();
			ite++;

			if (ite % 100 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_ball_to_vtu.WriteToFile(Real(ite) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.WriteToFile(0.0);
		return 0;
	}

	return 0;
}
