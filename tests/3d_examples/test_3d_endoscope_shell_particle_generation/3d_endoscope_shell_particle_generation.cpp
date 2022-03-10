/**
* @file 	test_3d_endoscope_shell_particle_generation.cpp
* @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation for an encoscope in stomach and esophagus case.
* @author 	Massoud Rezavand, Dong Wu and Xiangyu Hu
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_geometry = "./input/stomach_esophagus.stl";
std::string full_path_to_file_endoscope = "./input/endoscope.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-40.0, -100.0, -10.0);
Vec3d domain_upper_bound(40.0, 100.0, 830.0);
Real dp_0 = 1.5;
Real thickness = 1.23;
Real level_set_refinement_ratio = dp_0 / (0.2 * thickness);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	For material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
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

class endoscope : public SolidBody
{
public:
	endoscope(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_shape_stl(full_path_to_file_endoscope, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_shape_stl, true);
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
	system.reload_particles_ = true;
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
	ImportedShellModel imported_model(system, "EsophagusShellModel");
	ShellParticles imported_model_particles(imported_model,
											makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
											// makeShared<ParticleGeneratorReload>(in_output, imported_model.getBodyName()), //for reloading
											// or 
											makeShared<ShellParticleGeneratorLattice>(thickness),// for generation
											thickness);

	endoscope endoscope_model(system, "EndoscopeModel");
	ElasticSolidParticles endoscope_particles(endoscope_model, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(in_output, { imported_model});
	BodyStatesRecordingToVtp write_endoscope_to_vtp(in_output, { endoscope_model});
	ReloadParticleIO write_particle_reload_files(in_output, { &imported_model });
	ReloadParticleIO write_endoscope_reload_files(in_output, { &endoscope_model });
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
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the esophagus_stomach body N = " << ite_p << "\n";
			write_imported_model_to_vtp.writeToFile(ite_p);
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	write_imported_model_to_vtp.writeToFile(ite_p);
	std::cout << "The physics relaxation process of esophagus_stomach model finished !" << std::endl;
	/**
	 * @brief Particle generation for endoscope side
	 */
	BodyRelationInner endoscope_inner(endoscope_model);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_endoscope_particles(endoscope_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner endoscope_relaxation_step_inner(endoscope_inner, true);
	MeshRecordingToPlt write_endoscope_cell_linked_list(in_output, endoscope_model, endoscope_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_endoscope_particles.parallel_exec(0.25);
	endoscope_relaxation_step_inner.surface_bounding_.parallel_exec();
	write_endoscope_to_vtp.writeToFile(0.0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_2 = 0;
	while (ite_2 < 1000)
	{
		endoscope_relaxation_step_inner.parallel_exec();
		ite_2 += 1;
		if (ite_2 % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the endoscope model N = " << ite_2 << "\n";
			write_endoscope_to_vtp.writeToFile(ite_2);
		}
	}
	write_endoscope_to_vtp.writeToFile(ite_2);
	std::cout << "The physics relaxation process of endoscope model finished !" << std::endl;
	/** Output results. */
	write_endoscope_reload_files.writeToFile(0);
	write_particle_reload_files.writeToFile(0);

	return 0;
}


