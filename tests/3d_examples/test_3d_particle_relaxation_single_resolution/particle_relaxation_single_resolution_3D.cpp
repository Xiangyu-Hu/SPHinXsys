/**
 * @file 	particle_relaxation_single_resolution.cpp
 * @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
 * @details We use this case to test the particle generation and relaxation for a complex geometry. 
 *			Before particle generation, we clean the sharp corners of the model. 
 * @author 	Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/SPHinXsys.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-2.3, -0.1, -0.3);
Vec3d domain_upper_bound(2.3, 4.5, 0.3);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 80.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	define a body from the imported model.
//----------------------------------------------------------------------
class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_shape_stl(full_path_to_file, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_shape_stl, true);
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
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
	ImportedModel imported_model(system, "ImportedModel");
	SolidParticles imported_model_particles(imported_model);
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(in_output, {imported_model});
	MeshRecordingToPlt write_cell_linked_list(in_output, imported_model, imported_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner imported_model_inner(imported_model);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	write_imported_model_to_vtp.writeToFile(0.0);
	imported_model.updateCellLinkedList();
	write_cell_linked_list.writeToFile(0.0);
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
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			write_imported_model_to_vtp.writeToFile(ite_p);
		}
	}
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;

	return 0;
}
