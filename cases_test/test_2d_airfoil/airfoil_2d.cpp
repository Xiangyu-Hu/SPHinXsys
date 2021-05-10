/**
* @file 	airfoil_2d.cpp
* @brief 	This is the test of using levelset to generate body fitted SPH particles.
* @details	We use this case to test the particle generation and relaxation with a complex geometry (2D).
*			Before the particles are generated, we clean the sharp corners and other unresolvable surfaces.
* @author 	Yongchuan Yu and Xiangyu Hu
*/

#include "sphinxsys.h"

//Case file to setup the test case
#include "airfoil_2d.h"

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
	Airfoil* airfoil = new Airfoil(system, "Airfoil");
	SolidParticles airfoil_particles(airfoil);
	airfoil_particles.addAVariableToWrite<indexScalar, Real>("SmoothingLengthRatio");
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	
	WriteBodyStatesToVtu write_airfoil_to_vtu(in_output, { airfoil });
	WriteMeshToPlt 	write_mesh_cell_linked_list(in_output, airfoil, airfoil->mesh_cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BaseInnerBodyRelation* airfoil_inner = new InnerBodyRelationVariableSmoothingLength(airfoil);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_airfoil_particles(airfoil);
	relax_dynamics::RelaxationStepInner relaxation_step_inner(airfoil_inner, true);
	relax_dynamics::UpdateSmoothingLengthRatioByBodyShape update_smoothing_length_ratio(airfoil);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_airfoil_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	update_smoothing_length_ratio.parallel_exec();
	write_airfoil_to_vtu.WriteToFile(0.0);
	airfoil->updateCellLinkedList();
	write_mesh_cell_linked_list.WriteToFile(0.0);
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
			write_airfoil_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
		}
	}
	std::cout << "The physics relaxation process of airfoil finish !" << std::endl;

	return 0;
}
