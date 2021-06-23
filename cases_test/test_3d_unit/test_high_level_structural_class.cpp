#include <gtest/gtest.h>
#include "solid_structural_simulation_class.h"

Real tolerance = 1e-6;

TEST(StructuralSimulation, ExpandBoundingBox)
{
	BoundingBox bb(Vec3d(0), Vec3d(0));
	BoundingBox bb_2(Vec3d(-0.5, -1.0, -2.0), Vec3d(1.05, 1.1, 2.2));
	BoundingBox bb_3(Vec3d(1.0, -10.0, -20.0), Vec3d(10.05, 0.0, 0.0));

	ExpandBoundingBox(&bb, &bb_2);
	ExpandBoundingBox(&bb, &bb_3);

	BoundingBox bb_ref(Vec3d(-0.5, -10.0, -20.0), Vec3d(10.05, 1.1, 2.2));
	EXPECT_EQ(bb, bb_ref);
}

class TestStructuralSimulation : StructuralSimulation
{
public:
	TestStructuralSimulation(StructuralSimulationInput& input) : StructuralSimulation(input){};
	BoundingBox Get_system_domain_bounds_(){ return system_.system_domain_bounds_; };
};
TEST(StructuralSimulation, StructuralSimulation)
{
	/** INPUT PARAMETERS */
	Real scale_stl = 0.001 / 4; // diameter of 0.025 m
	Real resolution_mass = 8.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e4;
	Real physical_viscosity = 200;
	Real rho_0 = 1000;
	Real end_time = 0.1;
	Vecd final_position_center = Vecd(0, 0, 0.1);
	/** STL IMPORT PARAMETERS */
	std::string relative_input_path = "./input/"; //path definition for linux
	std::vector<std::string> imported_stl_list = { "ball_mass.stl" };
	std::vector<Vec3d> translation_list = { Vec3d(0) };
	std::vector<Real> resolution_list = { resolution_mass};
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	std::vector<LinearElasticSolid> material_model_list = { material };
	/** INPUT DECLERATION */
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
	Real scale_system_bounds = 10;
	input.scale_system_boundaries_ = scale_system_bounds;
	input.position_solid_body_tuple_ = { PositionSolidBodyTuple(0, end_time * 0.9, final_position_center ) };

	/** SIMULATION MODEL */
	TestStructuralSimulation sim(input);

	Real ball_radius = 100 * scale_stl * 0.5;
	BoundingBox test_bounds(Vec3d(-ball_radius * scale_system_bounds), Vec3d(ball_radius * scale_system_bounds));

	for (int i = 0; i < 3; i++)
	{
		EXPECT_NEAR(sim.Get_system_domain_bounds_().first[i], test_bounds.first[i], abs(test_bounds.first[i] * tolerance));
		EXPECT_NEAR(sim.Get_system_domain_bounds_().second[i], test_bounds.second[i], abs(test_bounds.first[i] * tolerance));
	}
	
	/** START SIMULATION */
	//sim.RunSimulation(end_time);
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
