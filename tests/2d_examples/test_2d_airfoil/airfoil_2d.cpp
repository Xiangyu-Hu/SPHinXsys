/**
* @file 	airfoil_2d.cpp
* @brief 	This is the test of using levelset to generate body fitted SPH particles.
* @details	We use this case to test the particle generation and relaxation with a complex geometry (2D).
*			Before the particles are generated, we clean the sharp corners and other unresolvable surfaces.
* @author 	Yongchuan Yu and Xiangyu Hu
*/

#include "sphinxsys.h"

#include "airfoil_2d.h" //Header to setup the test case

using namespace SPH;

int main(int ac, char *av[])
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
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	Airfoil airfoil(system, "Airfoil");
	SolidParticles airfoil_particles(airfoil, makeShared<ParticleGeneratorMultiResolution>());
	airfoil_particles.addAVariableToWrite<indexScalar, Real>("SmoothingLengthRatio");
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp airfoil_recording_to_vtp(in_output, {&airfoil});
	MeshRecordingToPlt cell_linked_list_recording(in_output, airfoil, airfoil.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInnerVariableSmoothingLength airfoil_inner(airfoil);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_airfoil_particles(airfoil);
	relax_dynamics::RelaxationStepInner relaxation_step_inner(airfoil_inner, true);
	relax_dynamics::UpdateSmoothingLengthRatioByBodyShape update_smoothing_length_ratio(airfoil);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	random_airfoil_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	update_smoothing_length_ratio.parallel_exec();
	airfoil.updateCellLinkedList();
	//----------------------------------------------------------------------
	//	First output before the simulation.
	//----------------------------------------------------------------------
	airfoil_recording_to_vtp.writeToFile(0);
	cell_linked_list_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 2000)
	{
		update_smoothing_length_ratio.parallel_exec();
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the airfoil N = " << ite_p << "\n";
			airfoil_recording_to_vtp.writeToFile(ite_p);
		}
	}
	std::cout << "The physics relaxation process of airfoil finish !" << std::endl;

	return 0;
}
