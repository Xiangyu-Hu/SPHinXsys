/**
 * @file 	collision.cpp
 * @brief 	a soft ball with initial velocity bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 20.0; /**< box length. */
Real DH = 13.0; /**< box height. */
Real L = 1.0;
Real slop_h = 11.55;
Real resolution_ref = L / 10.0; /**< reference particle spacing. */
Real BW = resolution_ref * 4;	/**< wall width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(25, 15));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e5;
Real poisson = 0.45;
Real gravity_g = 9.8;
Real physical_viscosity = 1000000.0;
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
// Wall boundary body definition.
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> wall_shape{Vecd(0, 0), Vecd(0, slop_h), Vecd(DL, slop_h), Vecd(0, 0)};
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(wall_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/** Definition of the ball as a elastic structure. */
class Cubic : public SolidBody
{
public:
	Cubic(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> cubic_shape;
		cubic_shape.push_back(Vecd(BW, slop_h + resolution_ref));
		cubic_shape.push_back(Vecd(BW, slop_h + L + resolution_ref));
		cubic_shape.push_back(Vecd(BW + L, slop_h + L + resolution_ref));
		cubic_shape.push_back(Vecd(BW + L, slop_h + resolution_ref));
		cubic_shape.push_back(Vecd(BW, slop_h + resolution_ref));
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(cubic_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Observer particle generator.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(7.2, 9.8), 0.0));
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
	/** tag for run particle relaxation for the initial body fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** tag for computation start with relaxed body fitted particles distribution */
	sph_system.reload_particles_ = true;
	/** tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	/** output environment. */
	In_Output in_output(sph_system);
/** handle command line arguments. */
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av);
#endif
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	Cubic free_cubic(sph_system, "FreeBall");
	SharedPtr<LinearElasticSolid> cubic_material = makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	ElasticSolidParticles free_cubic_particles(free_cubic, cubic_material);

	WallBoundary wall_boundary(sph_system, "Wall");
	SolidParticles solid_particles(wall_boundary, makeShared<Solid>(rho0_s, cubic_material->ContactStiffness()));

	ObserverBody free_cubic_observer(sph_system, "FreeBallObserver");
	ObserverParticles free_cubic_observer_particles(free_cubic_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner free_cubic_inner(free_cubic);
	SolidBodyRelationContact free_cubic_contact(free_cubic, {&wall_boundary});
	BodyRelationContact free_cubic_observer_contact(free_cubic_observer, {&free_cubic});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	TimeStepInitialization free_cubic_initialize_timestep(free_cubic, gravity);
	/** Kernel correction. */
	solid_dynamics::CorrectConfiguration free_cubic_corrected_configuration(free_cubic_inner);
	/** Time step size. */
	solid_dynamics::AcousticTimeStepSize free_cubic_get_time_step_size(free_cubic);
	/** stress relaxation for the solid body. */
	solid_dynamics::StressRelaxationFirstHalf free_cubic_stress_relaxation_first_half(free_cubic_inner);
	solid_dynamics::StressRelaxationSecondHalf free_cubic_stress_relaxation_second_half(free_cubic_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation free_cubic_update_contact_density(free_cubic_contact);
	solid_dynamics::ContactForce free_cubic_compute_solid_contact_forces(free_cubic_contact);
	/** Damping*/
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		damping(free_cubic_inner, 0.5, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	/** Observer and output. */
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexVector, Vecd>>
		write_free_cubic_displacement("Position", in_output, free_cubic_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	GlobalStaticVariables::physical_time_ = 0.0;
	Transformd transform(-0.5235, Vecd(0));
	solid_particles.ParticleTranslationAndRotation(transform);
	free_cubic_particles.ParticleTranslationAndRotation(transform);
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	free_cubic_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Initial states output.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile(0);
	write_free_cubic_displacement.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 2.5;
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
				free_cubic_initialize_timestep.parallel_exec();
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: " << dt
							  << "\n";
				}
				free_cubic_update_contact_density.parallel_exec();
				free_cubic_compute_solid_contact_forces.parallel_exec();
				free_cubic_stress_relaxation_first_half.parallel_exec(dt);
				free_cubic_stress_relaxation_second_half.parallel_exec(dt);

				free_cubic.updateCellLinkedList();
				free_cubic_contact.updateConfiguration();

				ite++;
				dt = free_cubic_get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_free_cubic_displacement.writeToFile(ite);
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

	write_free_cubic_displacement.newResultTest();

	return 0;
}
