/**
 * @file 	3d_self_contact.cpp
 * @brief 	This is the test to check self contact for solid dynamics
 * @author 	Xiangyu Hu
 */

#include "sphinxsys.h"
// case file to setup the test case
#include "3d_self_contact.h"

using namespace SPH;

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	system.run_particle_relaxation_ = false;
	// Tag for reload initially repaxed particles.
	system.reload_particles_ = true;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	// output environment
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	Coil coil(system, "Coil");
	SharedPtr<ParticleGenerator> coil_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		coil_particle_generator = makeShared<ParticleGeneratorReload>(in_output, coil.getBodyName());
	ElasticSolidParticles coil_particles(coil, makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson), coil_particle_generator);

	StationaryPlate stationary_plate(system, "StationaryPlate");
	SolidParticles moving_plate_particles(stationary_plate, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner coil_inner(coil);
	SolidBodyRelationSelfContact coil_self_contact(coil);
	SolidBodyRelationContact coil_contact(coil_self_contact, {&stationary_plate});
	//----------------------------------------------------------------------
	//	check whether run particle relaxation for body fitted particle distribution.
	//----------------------------------------------------------------------
	if (system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		// Random reset the insert body particle position.
		RandomizePartilePosition random_inserted_body_particles(coil);
		// Write the particle reload files.
		ReloadParticleIO write_particle_reload_files(in_output, {&coil});
		// A  Physics relaxation step.
		relax_dynamics::RelaxationStepInner relaxation_step_inner(coil_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_states.writeToFile(0);
		//----------------------------------------------------------------------
		//	Particle relaxation loop.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_states.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		// Output particles position for reload.
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	This section define all numerical methods will be used in this case.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -1.0, 0.0));
	// initialize a time step
	TimeStepInitialization initialization_with_gravity(coil, gravity);
	// Corrected configuration for reproducing rigid rotation.
	solid_dynamics::CorrectConfiguration corrected_configuration(coil_inner);
	// Time step size
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(coil);
	//stress relaxation.
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(coil_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(coil_inner);
	// Algorithms for solid-solid contacts.
	solid_dynamics::ContactDensitySummation coil_update_contact_density(coil_contact);
	solid_dynamics::ContactForce coil_compute_solid_contact_forces(coil_contact);
	solid_dynamics::SelfContactDensitySummation coil_self_contact_density(coil_self_contact);
	solid_dynamics::SelfContactForce coil_self_contact_forces(coil_self_contact);
	// Damping the velocity field for quasi-static solution
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>
		coil_damping(coil_inner, 0.2, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	From here the time stepping begins.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	// apply initial condition
	corrected_configuration.parallel_exec();
	write_states.writeToFile(0);
	// Setup time stepping control parameters.
	int ite = 0;
	Real end_time = 10.0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	// Statistics for computing time.
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			initialization_with_gravity.parallel_exec();
			// contact dynamics.
			coil_self_contact_density.parallel_exec();
			coil_self_contact_forces.parallel_exec();
			coil_update_contact_density.parallel_exec();
			coil_compute_solid_contact_forces.parallel_exec();
			// Stress relaxation and damping.
			stress_relaxation_first_half.parallel_exec(dt);
			coil_damping.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			//update particle neighbor relations for contact dynamics
			coil.updateCellLinkedList();
			coil_self_contact.updateConfiguration();
			coil_contact.updateConfiguration();
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
