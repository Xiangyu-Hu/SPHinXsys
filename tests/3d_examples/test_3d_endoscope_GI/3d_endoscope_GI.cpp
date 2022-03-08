/**
* @file 	test_3d_endoscope_GI.cpp
* @brief 	This is a test to simulate an endoscope stirred in the GI system.
* @details	The encoscope will be leaded through esophagus and stomach.
* @author 	Massoud Rezavand, Virtonomy GmbH
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
Real dp_0 = 0.5;
Real thickness = 1.0;	
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
	endoscope endoscope_model(system, "Endoscope");
	ElasticSolidParticles endoscope_particles(endoscope_model, 
											  makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
											  makeShared<ParticleGeneratorReload>(in_output, endoscope_model.getBodyName()));
	std::cout <<"Endoscope reloaded !" << std::endl;

	ImportedShellModel imported_model(system, "EsophagusStomach");
	ShellParticles imported_model_particles(imported_model,
											makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
											makeShared<ParticleGeneratorReload>(in_output, imported_model.getBodyName()),
											thickness);
	std::cout <<"EsophagusStomach reloaded !" << std::endl;
	
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_system_state_to_vtp(in_output, system.real_bodies_);
	write_system_state_to_vtp.writeToFile(0.0);

	return 0;
}


