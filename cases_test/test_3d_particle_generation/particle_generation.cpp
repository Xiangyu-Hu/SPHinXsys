 /**
 * @file 	particle_generation.cpp
 * @brief 	This is the test of using levelset to generate particles relax particles (3D).
 * @details We use this case to test the particle generation and relaxation by levelset for a complex geometry. 
 *			Before particle generation, we clean the sharp corners of the model. 

 * @author 	Yongchuan Yu and Xiangyu Hu
 * @version 0.1
 */

#include "sphinxsys.h"

/** case file to setup the test case */
#include "case.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(domain_lower_bound, domain_upper_bound, dp_0);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);

	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	ImportedModel* imported_model = new ImportedModel(system, "ImportedModel", 0);
	SolidParticles imported_model_particles(imported_model);

	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_real_body_states_to_vtu(in_output, system.real_bodies_);

	/**
	 * @brief 	Body relation map.
	 * @details The contact map gives the topological connections between the bodies.
	 * 			Basically the the range of bodies to build neighbor particle lists.
	 */
	SPHBodyInnerRelation* imported_model_inner = new SPHBodyInnerRelation(imported_model);

	/**
	 * @brief 	Methods used for particle relaxation.
	 */
	 /** Random reset the insert body particle position. */
	RandomizePartilePosition  random_inserted_body_particles(imported_model);
	/** Write the body state to Vtu file. */
	WriteBodyStatesToVtu 		write_inserted_body_to_vtu(in_output, { imported_model });
	/** Write the particle reload files. */
	WriteReloadParticle 		write_particle_reload_files(in_output, { imported_model });

	/** bounding particles to insert body surface. */
	relax_dynamics::BodySurfaceBounding
		body_surface_bounding(imported_model, new NearBodySurface(imported_model));

	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner);
	/** finalizing  particle number density and inital position after relaxatoin. */
	relax_dynamics::FinalizingParticleRelaxation finalizing_imported_model_particles(imported_model);
	/**
	  * @brief 	Particle relaxation starts here.
	  */
	body_surface_bounding.parallel_exec();
	random_inserted_body_particles.parallel_exec(0.25);
	write_real_body_states_to_vtu.WriteToFile(0.0);

	/** relax particles of the insert body. */
	int ite_p = 0;
	while (ite_p < 1000)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			cout << fixed << setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			write_inserted_body_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
		}
	}
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
	finalizing_imported_model_particles.parallel_exec();

	return 0;
}


