/**
 * @file 	3d_roof_analytical.cpp
 * @brief 	Shell verificaiton  incl. refinement study
 * @details Roof shell verification case with relaxed shell particles
 * @author 	Bence Rochlitz
 * @ref 	ANSYS Workbench Verification Manual, Release 15.0, November 2013, VMMECH069: Barrel Vault Roof Under Self Weight
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

Real to_rad(Real angle){return angle*Pi/180;}
static const int simtk_res = 1;

void relax_shell(RealBody& plate_body, Real thickness, Real level_set_refinement_ratio)
{
	// BUG: apparently only works if dp < thickness, otherwise ShellNormalDirectionPrediction::correctNormalDirection() throws error

	InnerRelation imported_model_inner(plate_body);
	SimpleDynamics<RandomizeParticlePosition> random_imported_model_particles(plate_body);
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(imported_model_inner, thickness, level_set_refinement_ratio);
	relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(imported_model_inner, thickness);	
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	plate_body.updateCellLinkedList();
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	shell_normal_prediction.exec();
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
}

int main(int ac, char *av[])
{
	// main geometric parameters
	Real radius = 25;
	Real length = 50;
	Real thickness = 0.25;
	Real teta = 40;
	Real arc = radius*to_rad(teta);
	// resolution
	Real dp = 2;

	// material
	Real rho = 36.7347;
	Real E = 4.32e8;
	Real mu = 0.3;
	auto material = makeShared<LinearElasticSolid>(rho, E, mu);

	{// generate particle positions
		// 1. to create the roof geometry we create the initial shell particles based on a plate
		// 2. then we warp the plate according to the given radius and angle

		// 1. Plate
		// fake thickness for fast levelset generation - we just need particle positions
		Real fake_th = 4;
		Real level_set_refinement_ratio = dp / (0.1 * fake_th);
		// plate dimensions
		Vec3d plate_dim(2*arc+dp, fake_th, length+dp);
		// shape
		auto plate_shape = makeShared<ComplexShape>("plate_shape");
		plate_shape->add<TriangleMeshShapeBrick>(plate_dim/2, simtk_res, Vec3d(0));
		SPHSystem system(plate_shape->getBounds(), dp);
		// body and particles
		RealBody plate_body(system, plate_shape);
		plate_body.defineBodyLevelSetShape(level_set_refinement_ratio)->correctLevelSetSign();
		plate_body.defineParticlesWithMaterial<ShellParticles>(material.get());
		plate_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(fake_th);
		relax_shell(plate_body, fake_th, level_set_refinement_ratio);
		// output
		IOEnvironment io_env(system);
		BodyStatesRecordingToVtp vtp_output(io_env, {plate_body});
		vtp_output.writeToFile();

		// 2. Warping

	}
	
	return 0;
}
