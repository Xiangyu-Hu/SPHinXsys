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
	ASSERT_EQ(bb, bb_ref);
}

class TestStructuralSimulation : StructuralSimulation
{
public:
	TestStructuralSimulation(StructuralSimulationInput& input) : StructuralSimulation(input){};
	void TestRunSimulation(Real end_time){ RunSimulation(end_time); };
	// input members
	string Get_relative_input_path_(){ return relative_input_path_; };
	vector<string> Get_imported_stl_list_(){ return imported_stl_list_; };
	Real Get_scale_stl_(){ return scale_stl_; };
	vector<Vec3d> Get_translation_list_(){ return translation_list_; };
	Real Get_system_resolution_(){ return system_resolution_; };
	vector<Real> Get_resolution_list_(){ return resolution_list_; };
	vector<LinearElasticSolid> Get_material_model_list_(){ return material_model_list_; };
	Real Get_physical_viscosity_(){ return physical_viscosity_; };
	// internal members
	SPHSystem Get_system_(){ return system_; };
	Real Get_scale_system_boundaries_(){ return scale_system_boundaries_; };
	In_Output Get_in_output_(){ return in_output_; };
	// other
	vector<TriangleMeshShape> Get_body_mesh_list_(){ return body_mesh_list_; };
	vector<shared_ptr<SolidBodyForSimulation>> Get_solid_body_list_(){ return solid_body_list_; };
	vector<array<int, 2>> Get_contacting_bodies_list_(){ return contacting_bodies_list_; };
	vector<shared_ptr<SolidContactBodyRelation>> Get_contact_list_(){ return contact_list_; };
	vector<shared_ptr<solid_dynamics::ContactDensitySummation>> Get_contact_density_list_(){ return contact_density_list_; };
	vector<shared_ptr<solid_dynamics::ContactForce>> Get_contact_force_list_(){ return contact_force_list_; };
	// for InitializeATimeStep
	vector<shared_ptr<TimeStepInitialization>> Get_initialize_gravity_(){ return initialize_gravity_; };
	vector<GravityPair> Get_non_zero_gravity_(){ return non_zero_gravity_; };
	// for AccelerationForBodyPartInBoundingBox
	vector<shared_ptr<solid_dynamics::AccelerationForBodyPartInBoundingBox>> Get_acceleration_bounding_box_(){ return acceleration_bounding_box_; };
	vector<AccelTuple> Get_acceleration_bounding_box_tuple_(){ return acceleration_bounding_box_tuple_; };
	// for SpringDamperConstraintParticleWise
	vector<shared_ptr<solid_dynamics::SpringDamperConstraintParticleWise>> Get_spring_damper_constraint_(){ return spring_damper_constraint_; };
	vector<SpringDamperTuple> Get_spring_damper_tuple_(){ return spring_damper_tuple_; };
	// for ConstrainSolidBodyRegion
	vector<shared_ptr<solid_dynamics::ConstrainSolidBodyRegion>> Get_fixed_constraint_(){ return fixed_constraint_; };
	vector<int> Get_body_indeces_fixed_constraint_(){ return body_indeces_fixed_constraint_; };
	// for PositionSolidBody
	vector<shared_ptr<solid_dynamics::PositionSolidBody>> Get_position_solid_body_(){ return position_solid_body_; };
	vector<PositionSolidBodyTuple> Get_position_solid_body_tuple_(){ return position_solid_body_tuple_; };
	// for PositionScaleSolidBody
	vector<shared_ptr<solid_dynamics::PositionScaleSolidBody>> Get_position_scale_solid_body_(){ return position_scale_solid_body_; };
	vector<PositionScaleSolidBodyTuple> Get_position_scale_solid_body_tuple_(){ return position_scale_solid_body_tuple_; };
};

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
	// test ScaleTranslationAndResolution();
	ASSERT_EQ(sim.Get_translation_list_().size(), sim.Get_resolution_list_().size());
	for (unsigned int i = 0; i < translation_list.size(); i++)
	{	
		ASSERT_EQ(sim.Get_translation_list_()[i], translation_list[i] * scale_stl);
		ASSERT_EQ(sim.Get_resolution_list_()[i], resolution_list[i] * scale_stl);
	}
	ASSERT_EQ(sim.Get_system_resolution_(), resolution_mass * scale_stl);
	//=================================================================================================//
	// test CreateBodyMeshList();
	ASSERT_EQ(sim.Get_body_mesh_list_().size(), number_of_bodies);
	//=================================================================================================//
	// test CalculateSystemBoundaries();
	Real ball_radius = 100 * scale_stl * 0.5;
	BoundingBox test_bounds(Vec3d(-ball_radius * scale_system_bounds), Vec3d(ball_radius * scale_system_bounds));
	for (int i = 0; i < 3; i++)
	{
		ASSERT_NEAR(sim.Get_system_().system_domain_bounds_.first[i], test_bounds.first[i], abs(test_bounds.first[i] * tolerance));
		ASSERT_NEAR(sim.Get_system_().system_domain_bounds_.second[i], test_bounds.second[i], abs(test_bounds.first[i] * tolerance));
	}
	//=================================================================================================//
	// test InitializeElasticSolidBodies();
	ASSERT_EQ(sim.Get_solid_body_list_().size(), number_of_bodies);
	//=================================================================================================//
	// test InitializeAllContacts();
	ASSERT_EQ(sim.Get_contacting_bodies_list_().size(), 0);
	ASSERT_EQ(sim.Get_contact_list_().size(), 0);
	ASSERT_EQ(sim.Get_contact_density_list_().size(), 0);
	ASSERT_EQ(sim.Get_contact_force_list_().size(), 0);
	//=================================================================================================//
	// test Boundary Conditions
	ASSERT_EQ(sim.Get_position_solid_body_().size(), 2);
	ASSERT_EQ(sim.Get_position_solid_body_tuple_().size(), 2);
	//=================================================================================================//

	//=================================================================================================//
	/** START SIMULATION */
	sim.TestRunSimulation(end_time);
	StdLargeVec<Vecd>& pos_0 = sim.Get_position_solid_body_()[0]->GetParticlePos0();
	StdLargeVec<Vecd>& pos_n = sim.Get_position_solid_body_()[0]->GetParticlePosN();

	for (unsigned int index = 0; index < pos_0.size(); index++)
	{
		for (int i = 0; i < 3; i++)
		{
			Vecd displ = pos_n[index] - pos_0[index];
			ASSERT_NEAR(displ[i], final_position_center[i], final_position_center.norm() * tolerance);
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

	StdLargeVec<Vecd>& pos_0 = sim.Get_position_scale_solid_body_()[0]->GetParticlePos0();
	StdLargeVec<Vecd>& pos_n = sim.Get_position_scale_solid_body_()[0]->GetParticlePosN();

	string name = "./input/cylinder.stl";
	TriangleMeshShape cylinder_mesh(name, translation_list[0] * scale_stl, scale_stl);
	BoundingBox bounding_box = cylinder_mesh.findBounds();
	Vec3d center = (bounding_box.first + bounding_box.second) * 0.5;

	for (unsigned int index = 0; index < pos_0.size(); index++)
	{
		for (int i = 0; i < 3; i++)
		{
			Vec3d displ = center + (pos_0[index] - center) * scale;
			ASSERT_NEAR(pos_n[index][i], displ[i], displ.norm() * tolerance);
		}
	}
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
