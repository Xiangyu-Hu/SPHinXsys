/**
* @file 	airfoil_2d.cpp
* @brief 	This is the test of using levelset to generate particles relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the re-initialization

* @author 	Yongchuan Yu and Xiangyu Hu
* @version 0.1
*/

#include "sphinxsys.h"

/** case file to setup the test case */
#include "airfoil_2d.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(Vec2d(-(DL1 + BW), -(DH + BW)), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
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
	Airfoil* airfoil = new Airfoil(system, "Airfoil", 0);
	SolidParticles airfoil_particles(airfoil);

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
	SPHBodyInnerRelation* airfoil_inner = new SPHBodyInnerRelation(airfoil);

	/**
	 * @brief 	Methods used for particle relaxation.
	 */
	 /** Random reset the insert body particle position. */
	RandomizePartilePosition  random_inserted_body_particles(airfoil);
	/** Write the body state to Vtu file. */
	WriteBodyStatesToVtu 		write_inserted_body_to_vtu(in_output, { airfoil });
	/** Write the particle reload files. */
	WriteReloadParticle 		write_particle_reload_files(in_output, { airfoil });

	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(airfoil_inner);
	/** finalizing  particle number density and inital position after relaxatoin. */
	relax_dynamics::FinalizingParticleRelaxation finalizing_airfoil_particles(airfoil);
	/**
	  * @brief 	Particle relaxation starts here.
	  */
	random_inserted_body_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	write_real_body_states_to_vtu.WriteToFile(0.0);
	/** relax particles of the insert body. */
	int ite_p = 0;
	while (ite_p < 1000)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			cout << fixed << setprecision(9) << "Relaxation steps for the airfoil N = " << ite_p << "\n";
			write_inserted_body_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
		}
	}
	std::cout << "The physics relaxation process of airfoil finish !" << std::endl;
	finalizing_airfoil_particles.parallel_exec();

	return 0;
}
