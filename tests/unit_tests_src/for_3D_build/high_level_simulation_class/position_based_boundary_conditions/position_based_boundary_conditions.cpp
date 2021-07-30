#include <gtest/gtest.h>
#include "test_structural_simulation_class.h"

TEST(StructuralSimulation, expandBoundingBox)
{
	BoundingBox bb(Vec3d(0), Vec3d(0));
	BoundingBox bb_2(Vec3d(-0.5, -1.0, -2.0), Vec3d(1.05, 1.1, 2.2));
	BoundingBox bb_3(Vec3d(1.0, -10.0, -20.0), Vec3d(10.05, 0.0, 0.0));

	expandBoundingBox(&bb, &bb_2);
	expandBoundingBox(&bb, &bb_3);
	
	BoundingBox bb_ref(Vec3d(-0.5, -10.0, -20.0), Vec3d(10.05, 1.1, 2.2));

	for (size_t i = 0; i < 3; i++)
	{
		EXPECT_EQ(bb.first[i], bb_ref.first[i]);
		EXPECT_EQ(bb.second[i], bb_ref.second[i]);
	}
}

TEST(StructuralSimulation, PositionSolidBodyTuple)
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
	int number_of_bodies = 1;
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
	input.position_solid_body_tuple_ = { PositionSolidBodyTuple(0, 0.0, end_time * 0.75, final_position_center ),
											PositionSolidBodyTuple(0, end_time * 0.75, end_time, final_position_center ) };
	//=================================================================================================//

	//=================================================================================================//
	/** SIMULATION MODEL */
	TestStructuralSimulation sim(input);
	//=================================================================================================//

	//=================================================================================================//
	// test scaleTranslationAndResolution();
	EXPECT_EQ(sim.get_translation_list_().size(), sim.get_resolution_list_().size());
	for (size_t i = 0; i < translation_list.size(); i++)
	{	
		EXPECT_EQ(sim.get_translation_list_()[i], translation_list[i] * scale_stl);
		EXPECT_EQ(sim.get_resolution_list_()[i], resolution_list[i] * scale_stl);
	}
	EXPECT_EQ(sim.get_system_resolution_(), resolution_mass * scale_stl);
	//=================================================================================================//
	// test createBodyMeshList();
	EXPECT_EQ(sim.get_body_mesh_list_().size(), number_of_bodies);
	//=================================================================================================//
	// test calculateSystemBoundaries();
	Real ball_radius = 100 * scale_stl * 0.5;
	BoundingBox test_bounds(Vec3d(-ball_radius * scale_system_bounds), Vec3d(ball_radius * scale_system_bounds));
	for (size_t i = 0; i < 3; i++)
	{
		EXPECT_NEAR(sim.get_system_().system_domain_bounds_.first[i], test_bounds.first[i], abs(test_bounds.first[i] * tolerance));
		EXPECT_NEAR(sim.get_system_().system_domain_bounds_.second[i], test_bounds.second[i], abs(test_bounds.first[i] * tolerance));
	}
	//=================================================================================================//
	// test InitializeElasticSolidBodies();
	EXPECT_EQ(sim.get_solid_body_list_().size(), number_of_bodies);
	//=================================================================================================//
	// test InitializeAllContacts();
	EXPECT_EQ(sim.get_contacting_body_pairs_list_().size(), 0);
	EXPECT_EQ(sim.get_contact_list_().size(), 0);
	EXPECT_EQ(sim.get_contact_density_list_().size(), 0);
	EXPECT_EQ(sim.get_contact_force_list_().size(), 0);
	//=================================================================================================//
	// test Boundary Conditions
	EXPECT_EQ(sim.get_position_solid_body_().size(), 2);
	EXPECT_EQ(sim.get_position_solid_body_tuple_().size(), 2);
	//=================================================================================================//

	//=================================================================================================//
	/** START SIMULATION */
	sim.TestRunSimulation(end_time);
	StdLargeVec<Vecd>& pos_0 = sim.get_position_solid_body_()[0]->GetParticlePos0();
	StdLargeVec<Vecd>& pos_n = sim.get_position_solid_body_()[0]->GetParticlePosN();

	for (size_t index = 0; index < pos_0.size(); index++)
	{
		for (size_t i = 0; i < 3; i++)
		{
			Vecd displ = pos_n[index] - pos_0[index];
			EXPECT_NEAR(displ[i], final_position_center[i], final_position_center.norm() * tolerance);
		}
	}
	//=================================================================================================//
}

TEST(StructuralSimulation, PositionScaleSolidBodyTuple)
{
	Real scale_stl = 0.001;
	Real resolution_cylinder = 0.4;
	Real end_time_simulation = 0.05;
	Real end_time_position = end_time_simulation * 0.9;
	Real rho_0 = 1000.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e4;
	Real physical_viscosity = 200;
	string relative_input_path = "./input/"; //path definition for linux
	string cylinder_stl = "cylinder.stl";
	vector<string> imported_stl_list = { cylinder_stl };
	vector<Vec3d> translation_list = { Vec3d(0) };
	vector<Real> resolution_list = { resolution_cylinder };
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	vector<LinearElasticSolid> material_model_list = { material };
	vector<array<int, 2>> contacting_bodies_list = {};
	StructuralSimulationInput input
	{
		relative_input_path,
		imported_stl_list,
		scale_stl,
		translation_list,
		resolution_list,
		material_model_list,
		physical_viscosity,
		contacting_bodies_list,
	};
	Real scale = 0.9;
	input.position_scale_solid_body_tuple_ = { PositionScaleSolidBodyTuple(0, 0.0, end_time_position, scale) };

	//=================================================================================================//
	TestStructuralSimulation sim (input);
	sim.TestRunSimulation(end_time_simulation);
	//=================================================================================================//

	StdLargeVec<Vecd>& pos_0 = sim.get_position_scale_solid_body_()[0]->GetParticlePos0();
	StdLargeVec<Vecd>& pos_n = sim.get_position_scale_solid_body_()[0]->GetParticlePosN();

	string name = "./input/cylinder.stl";
	TriangleMeshShape cylinder_mesh(name, translation_list[0] * scale_stl, scale_stl);
	BoundingBox bounding_box = cylinder_mesh.findBounds();
	Vec3d center = (bounding_box.first + bounding_box.second) * 0.5;

	for (size_t index = 0; index < pos_0.size(); index++)
	{
		for (size_t i = 0; i < 3; i++)
		{
			Vec3d displ = center + (pos_0[index] - center) * scale;
			EXPECT_NEAR(pos_n[index][i], displ[i], displ.norm() * tolerance);
		}
	}
}

TEST(StructuralSimulation, TranslateSolidBodyTuple)
{
	Real scale_stl = 0.001 / 4; // diameter of 0.025 m
	Real resolution_mass = 8.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e4;
	Real physical_viscosity = 200;
	Real rho_0 = 1000;
	Real end_time = 0.1;

		/** STL IMPORT PARAMETERS */
	std::string relative_input_path = "./input/"; //path definition for linux
	std::vector<std::string> imported_stl_list = { "ball_mass.stl" };
	std::vector<Vec3d> translation_list = { Vec3d(0) };
	std::vector<Real> resolution_list = { resolution_mass};
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	std::vector<LinearElasticSolid> material_model_list = { material };
	
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
	Vecd translation_vector = Vec3d(0.0, 0.0, 0.1);
	input.translation_solid_body_tuple_ = { TranslateSolidBodyTuple(0, end_time * 0.32, end_time, translation_vector) };

	//=================================================================================================//
	TestStructuralSimulation sim (input);
	sim.TestRunSimulation(end_time);
	//=================================================================================================//

StdLargeVec<Vecd>& pos_0 = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_0_;
StdLargeVec<Vecd>& pos_n = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_n_;

	for (size_t index = 0; index < pos_0.size(); index++)
	{
		Vec3d end_pos = pos_0[index] + translation_vector;
		EXPECT_NEAR(pos_n[index][2], end_pos[2], end_pos.norm() * 1e-2);
	}
}

TEST(StructuralSimulation, TranslateSolidBodyPartTuple)
{
	Real scale_stl = 0.001 / 4; // diameter of 0.025 m
	Real resolution_mass = 8.0;
	Real poisson = 0.35;
	Real Youngs_modulus = 1e4;
	Real physical_viscosity = 200;
	Real rho_0 = 1000;
	Real end_time = 0.1;

	/** STL IMPORT PARAMETERS */
	std::string relative_input_path = "./input/"; //path definition for linux
	std::vector<std::string> imported_stl_list = { "ball_mass.stl" };
	std::vector<Vec3d> translation_list = { Vec3d(0) };
	std::vector<Real> resolution_list = { resolution_mass};
	LinearElasticSolid material = LinearElasticSolid(rho_0, Youngs_modulus, poisson);
	std::vector<LinearElasticSolid> material_model_list = { material };

	TriangleMeshShape ball_mesh("./input/ball_mass.stl", translation_list[0] * scale_stl, scale_stl);
	BoundingBox bbox = ball_mesh.findBounds();
	Real z_limit = 0.75 * bbox.first[2] + 0.25 * bbox.second[2]; // only apply to the bottom 25% of the ball
	bbox.second[2] = z_limit;
	
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
	Vecd translation_vector = Vec3d(0.0, 0.0, 0.02);
	input.translation_solid_body_part_tuple_ = { TranslateSolidBodyPartTuple(0, end_time * 0.124, end_time, translation_vector, bbox) };

	//=================================================================================================//
	TestStructuralSimulation sim (input);
	sim.TestRunSimulation(end_time);
	//=================================================================================================//

	StdLargeVec<Vecd>& pos_0 = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_0_;
	StdLargeVec<Vecd>& pos_n = sim.get_solid_body_list_()[0].get()->getElasticSolidParticles()->pos_n_;

	for (size_t index = 0; index < pos_0.size(); index++)
	{
		if (pos_0[index][2] < z_limit)
		{
			Vec3d end_pos = pos_0[index] + translation_vector;
			EXPECT_NEAR(pos_n[index][2], end_pos[2], end_pos.norm() * 0.05);
		}
		else
		{
			Real z_limit_end  = z_limit + translation_vector[2] * 0.95; // z_limit + translation_vector[2] = 0.01375 is the actual limit, but some particle go below this level
			EXPECT_GT(pos_n[index][2], z_limit_end);
		}
	}
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
