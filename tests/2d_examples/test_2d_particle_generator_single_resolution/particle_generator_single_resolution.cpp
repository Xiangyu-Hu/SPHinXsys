/**
* @file 	particle_generator_single_resolution.cpp
* @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the re-initialization

* @author 	Yongchuan Yu and Xiangyu Hu
* @version 0.1
*/

#include "sphinxsys.h"

/** case file to setup the test case */
#include "particle_generator_single_resolution.h"

using namespace SPH;

int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	//handle command line arguments
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	/** output environment. */
	In_Output 	in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	InputBody* inputbody = new InputBody(system, "Airfoil");
	SolidParticles inputbody_particles(inputbody);
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtu inputbody_recording_to_vtu(in_output, { inputbody });
	MeshRecordingToPlt 	cell_linked_list_recording(in_output, inputbody, inputbody->cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BaseBodyRelationInner* inputbody_inner = new BodyRelationInner(inputbody);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_inputbody_particles(inputbody);
	relax_dynamics::RelaxationStepInner relaxation_step_inner(inputbody_inner, true);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	random_inputbody_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	inputbody->updateCellLinkedList();
	//----------------------------------------------------------------------
	//	First output before the simulation.
	//----------------------------------------------------------------------
	inputbody_recording_to_vtu.writeToFile(0.0);
	cell_linked_list_recording.writeToFile(0.0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the airfoil N = " << ite_p << "\n";
			inputbody_recording_to_vtu.writeToFile(ite_p);
		}
	}
	std::cout << "The physics relaxation process of airfoil finish !" << std::endl;

	return 0;
}
