/**
 * @file 	collision.cpp
 * @brief 	two soft balls with and without internal damping bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"	//SPHinXsys Library.
using namespace SPH;	//Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0; 					/**< box length. */
Real DH = 4.0; 					/**< box height. */
Real resolution_ref = 0.025; 	/**< reference resolution. */
Real BW = resolution_ref * 4; 	/**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center_1(2.0, 2.0);
Vec2d ball_center_2(6.0, 2.0);
Real ball_radius = 0.5;			
Real gravity_g = 1.0;
//----------------------------------------------------------------------
//	Global paramters on material properties
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
	WallBoundary(SPHSystem &sph_system, std::string body_name) : SolidBody(sph_system, body_name)
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

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
class FreeBall : public SolidBody
{
public:
	FreeBall(SPHSystem& system, std::string body_name) : SolidBody(system, body_name)
	{
		ComplexShape original_body_shape;
		original_body_shape.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
class DampingBall : public SolidBody
{
public:
	DampingBall(SPHSystem& system, std::string body_name) : SolidBody(system, body_name)
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addACircle(ball_center_2, ball_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
//----------------------------------------------------------------------
//	Case dependent material properties
//----------------------------------------------------------------------
class WallMaterial : public LinearElasticSolid
{
public:
	WallMaterial() : LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
class BallMaterial : public NeoHookeanSolid
{
public:
	BallMaterial() : NeoHookeanSolid()
	{
		rho0_ 	= rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
//----------------------------------------------------------------------
//	Observer bodies with cases-dependent observation points.
//	The observation particles have zero volume.
//----------------------------------------------------------------------
class FreeBallObserver : public FictitiousBody
{
public:
	FreeBallObserver(SPHSystem& system, std::string body_name) : FictitiousBody(system, body_name)
	{
		body_input_points_volumes_.push_back(std::make_pair(ball_center_1, 0.0));
	}
};
class DampingBallObserver : public FictitiousBody
{
public:
	DampingBallObserver(SPHSystem& system, std::string body_name) : FictitiousBody(system, body_name)
	{
		body_input_points_volumes_.push_back(std::make_pair(ball_center_2, 0.0));
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
	sph_system.reload_particles_ = true;
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
	WallBoundary *wall_boundary = new WallBoundary(sph_system, "Wall");
	WallMaterial* wall_material = new WallMaterial();
	SolidParticles 	solid_particles(wall_boundary, wall_material);

	FreeBall* free_ball = new FreeBall(sph_system, "FreeBall");
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_) free_ball->useParticleGeneratorReload();
	BallMaterial* free_ball_material = new BallMaterial();
	ElasticSolidParticles 	free_ball_particles(free_ball, free_ball_material);

	DampingBall* damping_ball = new DampingBall(sph_system, "DampingBall");
	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_) damping_ball->useParticleGeneratorReload();
	BallMaterial* damping_ball_material = new BallMaterial();
	ElasticSolidParticles 	damping_ball_particles(damping_ball, damping_ball_material);

	FreeBallObserver* free_ball_observer = new FreeBallObserver(sph_system, "FreeBallObserver");
	BaseParticles 	free_ball_observer_particles(free_ball_observer);
	DampingBallObserver* damping_ball_observer = new DampingBallObserver(sph_system, "DampingBallObserver");
	BaseParticles 	damping_ball_observer_particles(damping_ball_observer);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner* free_ball_inner = new BodyRelationInner(free_ball);
		BodyRelationInner* damping_ball_inner = new BodyRelationInner(damping_ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		RandomizePartilePosition  			free_ball_random_particles(free_ball);
		RandomizePartilePosition  			damping_ball_random_particles(damping_ball);
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		relax_dynamics::RelaxationStepInner damping_ball_relaxation_step_inner(damping_ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtu 		write_ball_state(in_output, sph_system.real_bodies_);
		ReloadParticleIO 		write_particle_reload_files(in_output, { free_ball, damping_ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		free_ball_random_particles.parallel_exec(0.25);
		damping_ball_random_particles.parallel_exec(0.25);
		write_ball_state.writeToFile(0.0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begines.
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
	BodyRelationInner*   free_ball_inner = new BodyRelationInner(free_ball);
	SolidBodyRelationContact* free_ball_contact = new SolidBodyRelationContact(free_ball, {wall_boundary});
	BodyRelationInner*   damping_ball_inner = new BodyRelationInner(damping_ball);
	SolidBodyRelationContact* damping_ball_contact = new SolidBodyRelationContact(damping_ball, {wall_boundary});
	BodyRelationContact* free_ball_observer_contact = new BodyRelationContact(free_ball_observer, { free_ball });
	BodyRelationContact* damping_all_observer_contact = new BodyRelationContact(damping_ball_observer, { damping_ball });
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	TimeStepInitialization 	free_ball_initialize_timestep(free_ball, &gravity);
	TimeStepInitialization 	damping_ball_initialize_timestep(damping_ball, &gravity);
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
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		damping(damping_ball_inner, 0.5, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtu	body_states_recording(in_output, sph_system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		free_ball_displacement_recording("Position", in_output, free_ball_observer_contact);
	ObservedQuantityRecording<indexVector, Vecd>
		damping_ball_displacement_recording("Position", in_output, damping_all_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	free_ball_corrected_configuration.parallel_exec();
	damping_ball_corrected_configuration.parallel_exec();
	/** Initial states output. */
	body_states_recording.writeToFile(0);
	free_ball_displacement_recording.writeToFile(0);
	damping_ball_displacement_recording.writeToFile(0);
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
	return 0;
}
