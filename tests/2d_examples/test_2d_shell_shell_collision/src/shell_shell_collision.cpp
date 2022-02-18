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
Real DL = 4.0; 					/**< box length. */
Real DH = 4.0; 					/**< box height. */
Real resolution_ref = 0.025; 	/**< reference resolution. */
Real BW = resolution_ref * 1.; 	/**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center(2.0, 0.55);
Real ball_radius = 0.5;			
Real gravity_g = 0.05;
//----------------------------------------------------------------------
//	Global paramters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0;					/** Normalized density. */
Real rho0_ball = 0.1;				/** Normalized density. */
Real Youngs_modulus = 5e4;			/** Normalized Youngs Modulus. */
Real poisson = 0.45;				/** Poisson ratio. */
Real physical_viscosity = 200.0;	/** physical damping, here we choose the same value as numerical viscosity. */
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
		for (Real t = 0; t <= 4.*Pi*ball_radius; t += 1.75 * resolution_ref) 
		{
            Real x = ball_radius * cos(t) + ball_center[0];
            Real y = ball_radius * sin(t) + ball_center[1];
			positions_volumes_.push_back(std::make_pair(Vecd(x, y), resolution_ref));
		}
	}
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
										makeShared<BallParticleGenerator>(), BW);

	WallBoundary wall_boundary(sph_system, "Wall");
	SharedPtr<ParticleGenerator> wall_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		wall_particle_generator = makeShared<ParticleGeneratorReload>(in_output, wall_boundary.getBodyName());
	SharedPtr<LinearElasticSolid> wall_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	ShellParticles solid_particles(wall_boundary, wall_material, wall_particle_generator, BW);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner ball_inner(ball);
	BodyRelationInner wall_inner(wall_boundary);
	SolidBodyRelationContact ball_contact(ball, {&wall_boundary});
	SolidBodyRelationContact wall_ball_contact(wall_boundary, {&ball});
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
	/** Geometry and generation of the holder. */
	std::vector<Vecd> holder_shape;
		holder_shape.push_back(Vecd(-BW, -BW));
		holder_shape.push_back(Vecd(-BW, BW));
		holder_shape.push_back(Vecd(DL/4., BW));
		holder_shape.push_back(Vecd(DL/4., -BW));
		holder_shape.push_back(Vecd(-BW, -BW));
	std::vector<Vecd> holder_shape2;
		holder_shape2.push_back(Vecd(3.*DL/4., -BW));
		holder_shape2.push_back(Vecd(3.*DL/4., BW));
		holder_shape2.push_back(Vecd(DL + BW, BW));
		holder_shape2.push_back(Vecd(DL + BW, -BW));
		holder_shape2.push_back(Vecd(3.*DL/4., -BW));
	MultiPolygon multi_polygon_holder;
	multi_polygon_holder.addAPolygon(holder_shape, ShapeBooleanOps::add);
	multi_polygon_holder.addAPolygon(holder_shape2, ShapeBooleanOps::add);
	MultiPolygonShape holder_multibody_shape(multi_polygon_holder);
	BodyRegionByParticle holder(wall_boundary, "Holder", holder_multibody_shape);
	thin_structure_dynamics::ConstrainShellBodyRegion	constrain_holder(wall_boundary, holder);
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
