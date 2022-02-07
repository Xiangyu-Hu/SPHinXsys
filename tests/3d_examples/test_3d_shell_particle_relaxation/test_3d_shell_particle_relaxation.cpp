/**
* @file 	test_3d_shell_particle_relaxation.cpp
* @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex thin structures geometry (3D).
* @author 	Dong Wu and Xiangyu Hu
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_geometry = "./input/SPHinXsys.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(12, 14, 446);
Vec3d domain_upper_bound(1315, 1317, 1302);
Real dp_0 = 25.0;
Real thickness = 50.0;
Real level_set_refinement_ratio = dp_0 / (0.1 * thickness);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	For material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 1.3024653e6;	                     /** Normalized Youngs Modulus. */
Real poisson = 0.3; 			                         /** Poisson ratio. */
//----------------------------------------------------------------------
//	Define the body.
//----------------------------------------------------------------------
class ImportedShellModel : public ThinStructure
{
public:
	ImportedShellModel(SPHSystem &system, const std::string body_name)
		: ThinStructure(system, body_name, makeShared<SPHAdaptation>(1.15, 1.0, 0.75, level_set_refinement_ratio))
	{
		/** Geometry definition. */
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_geometry_shape(full_path_to_geometry, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_geometry_shape, true, false);
	}
};
//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	//----------------------------------------------------------------------
	//	Tag for run particle relaxation for the initial body fitted distribution.
	//----------------------------------------------------------------------
	system.run_particle_relaxation_ = true;
	//----------------------------------------------------------------------
	//	handle command line arguments.
	//----------------------------------------------------------------------
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	/** output environment. */
	In_Output 	in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	ImportedShellModel imported_model(system, "ImportedShellModel");
	ShellParticles imported_model_particles(imported_model,
											makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
											makeShared<ShellParticleGeneratorLattice>(thickness), thickness);
	imported_model_particles.addAVariableToWrite<Vecd>("InitialNormalDirection");
	imported_model_particles.addAVariableToWrite<Vecd>("NormalDirection");
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(in_output, { imported_model });
	MeshRecordingToPlt 	write_mesh_cell_linked_list(in_output, imported_model, imported_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner imported_model_inner(imported_model);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(imported_model_inner, thickness, level_set_refinement_ratio);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	write_imported_model_to_vtp.writeToFile(0.0);
	imported_model.updateCellLinkedList();
	write_mesh_cell_linked_list.writeToFile(0.0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
			write_imported_model_to_vtp.writeToFile(ite_p);
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	relaxation_step_inner.mid_surface_bounding_.calculateNormalDirection(); 
	write_imported_model_to_vtp.writeToFile(ite_p);
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;

	return 0;
}


