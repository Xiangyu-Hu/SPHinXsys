/**
 * @file 	shell_shell_collision.cpp
 * @brief 	An elastic shell ball bouncing within a confined shell boundary
 * @details This is a case to test shell->shell collision without impact.
 * @details Both the ball and box are thin shell structures.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h"	//SPHinXsys Library.
using namespace SPH;	//Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
// Real DL = 4.0; 					/**< box length. */
// Real DH = 4.0; 					/**< box height. */
// Real resolution_ref = 0.025; 	/**< reference resolution. */
// Real BW = resolution_ref * 1.; 	/**< wall width for BCs. */
// BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

Real ball_radius = 1.25;			
Real gravity_g = 0.05;

Real PT = 1.0;									  /** Thickness of the square plate. */
Vec2d n_0 = Vec2d(0.0, 1.0);					  /** Pseudo-normal. */
int particle_number = 40;						  /** Particle number in the direction of the length */
int BWD = 1;
Real PL = 10.0;
Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
Real BW = resolution_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
BoundingBox system_domain_bounds(Vec2d(-PL, -PL),
								 Vec2d(PL + BW, PL));
Vec2d ball_center(5. + resolution_ref / 2., 1.6);

//----------------------------------------------------------------------
//	Global paramters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0;					/** Normalized density. */
Real rho0_ball = 1.0;				/** Normalized density. */
Real Youngs_modulus = 5e4;			/** Normalized Youngs Modulus. */
Real poisson = 0.45;				/** Poisson ratio. */
Real physical_viscosity = 200.0;	/** physical damping, here we choose the same value as numerical viscosity. */
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
// class WallBoundary : public ThinStructure
// {
// public:
// 	WallBoundary(SPHSystem &sph_system, std::string body_name) 
// 	: ThinStructure(sph_system, body_name)
// 	{
// 		std::vector<Vecd> outer_wall_shape;
// 		outer_wall_shape.push_back(Vecd(-BW, -BW));
// 		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
// 		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
// 		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
// 		outer_wall_shape.push_back(Vecd(-BW, -BW));

// 		std::vector<Vecd> inner_wall_shape;
// 		inner_wall_shape.push_back(Vecd(0.0, 0.0));
// 		inner_wall_shape.push_back(Vecd(0.0, DH));
// 		inner_wall_shape.push_back(Vecd(DL, DH));
// 		inner_wall_shape.push_back(Vecd(DL, 0.0));
// 		inner_wall_shape.push_back(Vecd(0.0, 0.0));

// 		MultiPolygon multi_polygon;
// 		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
// 		multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
// 		body_shape_.add<MultiPolygonShape>(multi_polygon);
// 	}
// };
/** Define application dependent particle generator for thin structure. */
class BallParticleGenerator : public ParticleGeneratorDirect
{
public:
	BallParticleGenerator() : ParticleGeneratorDirect()
	{
		for (Real t = 0; t <= 1.65*Pi*ball_radius; t += 0.75 * resolution_ref) 
		{
            Real x = ball_radius * cos(t) + ball_center[0];
            Real y = ball_radius * sin(t) + ball_center[1];
			positions_volumes_.push_back(std::make_pair(Vecd(x, y), resolution_ref));
		}
	}
};

/** Define application dependent particle generator for thin structure. */
class PlateParticleGeneratorWall : public ParticleGeneratorDirect
{
public:
	PlateParticleGeneratorWall() : ParticleGeneratorDirect()
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			Real x = resolution_ref * i - BW + resolution_ref * 0.5;
			positions_volumes_.push_back(std::make_pair(Vecd(x, 0.0), resolution_ref));
		}
	}
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_->pos_n_[index_i][0] < 0.0 || base_particles_->pos_n_[index_i][0] > PL)
		{
			body_part_particles_.push_back(index_i);
		}
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
	ThinStructure ball(sph_system, "Ball", makeShared<SPHAdaptation>(1.15, 1.0));
	ShellParticles ball_particles(ball,
										makeShared<LinearElasticSolid>(rho0_ball, Youngs_modulus, poisson),
										makeShared<BallParticleGenerator>(), PT);

		/** Creat a plate body. */
	ThinStructure wall_boundary(sph_system, "Wall", makeShared<SPHAdaptation>(1.15, 1.0));
	/** Creat particles for the elastic body. */
	ShellParticles wall_particles(wall_boundary,
										makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
										makeShared<PlateParticleGeneratorWall>(), PT);
	wall_particles.addAVariableToWrite<Vecd>("PriorAcceleration");


	// WallBoundary wall_boundary(sph_system, "Wall");
	// SharedPtr<ParticleGenerator> wall_particle_generator = makeShared<ParticleGeneratorLattice>();
	// if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
	// 	wall_particle_generator = makeShared<ParticleGeneratorReload>(in_output, wall_boundary.getBodyName());
	// SharedPtr<LinearElasticSolid> wall_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	// ShellParticles wall_particles(wall_boundary, wall_material, wall_particle_generator, BW);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner ball_inner(ball);
	BodyRelationInner wall_inner(wall_boundary);
	SolidBodyRelationContact ball_contact(ball, {&wall_boundary});
	SolidBodyRelationContact wall_ball_contact(wall_boundary, {&ball});
	// wall_particles.addAVariableToWrite<Vec3d>("PriorAcceleration");
	ball_particles.addAVariableToWrite<Vec3d>("PriorAcceleration");
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	TimeStepInitialization 	ball_initialize_timestep(ball, gravity);
	TimeStepInitialization 	wall_initialize_timestep(wall_boundary);
	thin_structure_dynamics::ShellCorrectConfiguration ball_corrected_configuration(ball_inner);
	thin_structure_dynamics::ShellCorrectConfiguration wall_corrected_configuration(wall_inner);
	thin_structure_dynamics::ShellAcousticTimeStepSize ball_get_time_step_size(ball);
	/** stress relaxation for the balls. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf ball_stress_relaxation_first_half(ball_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf ball_stress_relaxation_second_half(ball_inner);
	thin_structure_dynamics::ShellStressRelaxationFirstHalf wall_stress_relaxation_first_half(wall_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf wall_stress_relaxation_second_half(wall_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity ball_update_contact_density(ball_contact);
	solid_dynamics::ShellContactDensity wall_ball_update_contact_density(wall_ball_contact);
	solid_dynamics::ShellShellContactForce ball_compute_solid_contact_forces(ball_contact);
	solid_dynamics::ShellShellContactForce wall_compute_solid_contact_forces(wall_ball_contact);
	/** Damping */
	DampingWithRandomChoice<DampingPairwiseInner<Vec2d>>
		ball_position_damping(ball_inner, 0.2, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<Vec2d>>
		ball_rotation_damping(ball_inner, 0.2, "AngularVelocity", physical_viscosity);

	DampingWithRandomChoice<DampingPairwiseInner<Vec2d>>
		wall_position_damping(wall_inner, 0.2, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<Vec2d>>
		wall_rotation_damping(wall_inner, 0.2, "AngularVelocity", physical_viscosity);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(wall_boundary, "BoundaryGeometry");
	thin_structure_dynamics::ConstrainShellBodyRegion constrain_holder(wall_boundary, boundary_geometry);
	// thin_structure_dynamics::ConstrainShellBodyRegion	constrain_holder(wall_boundary, holder);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp	body_states_recording(in_output, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	ball_corrected_configuration.parallel_exec();
	wall_corrected_configuration.parallel_exec();
	
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
				wall_initialize_timestep.parallel_exec();
				if (ite % 100 == 0) 
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				ball_update_contact_density.parallel_exec();
				ball_compute_solid_contact_forces.parallel_exec();

				wall_ball_update_contact_density.parallel_exec();
				wall_compute_solid_contact_forces.parallel_exec();

				ball_stress_relaxation_first_half.parallel_exec(dt);
				ball_position_damping.parallel_exec(dt);
				ball_rotation_damping.parallel_exec(dt);
				ball_stress_relaxation_second_half.parallel_exec(dt);
				
				wall_stress_relaxation_first_half.parallel_exec(dt);
				constrain_holder.parallel_exec(dt);
				wall_position_damping.parallel_exec(dt);
				wall_rotation_damping.parallel_exec(dt);
				constrain_holder.parallel_exec(dt);
				wall_stress_relaxation_second_half.parallel_exec(dt);

				ball.updateCellLinkedList();
				ball_contact.updateConfiguration();
				wall_boundary.updateCellLinkedList();
				wall_ball_contact.updateConfiguration();

				ite++;
				Real dt_free = ball_get_time_step_size.parallel_exec();
				dt = dt_free;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
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
	return 0;
}
