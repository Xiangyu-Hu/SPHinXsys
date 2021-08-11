#include <gtest/gtest.h>
#include "test_structural_simulation_class.h"

TEST(SoftCubeContact, NeoHookean)
{
	Real scale_stl = 0.001;
	Real end_time = 0.4;

	Real rho_0 = 1265.0; // Gheorghe 2019 
	Real poisson = 0.45; // nearly incompressible
	Real Youngs_modulus = 5e4; // Sommer 2015
	Real physical_viscosity = 200.0; //physical damping
	Real pressure = 6250; // 10N force on the surface
	int particles_in_thickness = 8;
	Real res = 40.0 / particles_in_thickness;

	/** STL IMPORT PARAMETERS */
	std::string relative_input_path = "./input/"; //path definition for linux
	std::vector<std::string> imported_stl_list = { "cube.stl" };
	std::vector<Vec3d> translation_list = { Vec3d(0) };
	std::vector<Real> resolution_list = { res };
	NeoHookeanSolid material = NeoHookeanSolid(rho_0, Youngs_modulus, poisson);
	std::vector<LinearElasticSolid> material_model_list = { material };

	TriangleMeshShape cube("./input/cube.stl", Vec3d(0), scale_stl);
	BoundingBox bbox = cube.findBounds();

	BoundingBox fixed_part = bbox;
	Real x_fix = bbox.first[0] + 0.01; // left part is fixed
	fixed_part.second[0] = x_fix;

	BoundingBox force_part = bbox;
	Real x_min = bbox.second[0] - 2 * res; // right side: force application, 2 particle layers
	force_part.first[0] = x_min;
	
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		{}
	};
	input.body_indeces_fixed_constraint_region_ = vector<ConstrainedRegionPair>{ ConstrainedRegionPair(0, fixed_part) };
	input.surface_pressure_tuple_ = StdVec<PressureTuple>{ PressureTuple(0, pressure, Vec3d(0.1, 0.02, 0.02), end_time * 0.1) };
	input.particle_relaxation_list_ = { false };

	//=================================================================================================//
	TestStructuralSimulation sim(input);
	sim.TestRunSimulation(end_time);

	// check which particles have pressure
	// only one layer of particles on one side should have
	int particles_ref = particles_in_thickness * particles_in_thickness;
	StdLargeVec<bool>& apply_pressure_to_particle_ = sim.get_surface_pressure_()[0].get()->GetApplyPressureToParticle();
	int particles_with_pressure = 0;
	for (size_t i = 0; i < apply_pressure_to_particle_.size(); i++)
	{
		if (apply_pressure_to_particle_[i]) particles_with_pressure++;
	}
	EXPECT_EQ(particles_ref, particles_with_pressure);

	/* this part */
	Real displ_max = sim.getMaxDisplacement(0);
	Real displ_max_ref = 0.00621; // in mm, absolute max displacement, this reference is the first simulation results, not validated solution
	EXPECT_NEAR(displ_max, displ_max_ref, displ_max_ref * 0.05); // 5% tolerance
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
