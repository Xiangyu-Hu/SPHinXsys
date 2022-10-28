/**
 * @file 	freestream_flow_around_cylinder_case.h
 * @brief 	This is the case file for the test of free-stream flow.
 * @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
#include "mr_free_stream_around_cylinder.h"

using namespace SPH;

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	sph_system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	sph_system.restart_step_ = 0;
	/** handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineComponentLevelSetShape("OuterBoundary");
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
	cylinder.defineAdaptationRatios(1.15, 2.0);
	cylinder.defineBodyLevelSetShape();
	cylinder.defineParticlesAndMaterial<SolidParticles, Solid>();
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
		: cylinder.generateParticles<ParticleGeneratorLattice>();

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex(water_block_inner, {&cylinder});
	ContactRelation cylinder_contact(cylinder, {&water_block});
	ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		/** body topology only for particle relaxation */
		InnerRelation cylinder_inner(cylinder);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
		SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
		BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
		ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
		relax_dynamics::RelaxationStepComplex relaxation_step_complex(water_block_complex, "OuterBoundary", true);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.parallel_exec(0.25);
		random_water_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		relaxation_step_complex.surface_bounding_.parallel_exec();
		write_real_body_states.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			relaxation_step_complex.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_real_body_states.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		/** Output results. */
		write_real_body_particle_reload_files.writeToFile(0);
		return 0;
	}

	return 0;
}
