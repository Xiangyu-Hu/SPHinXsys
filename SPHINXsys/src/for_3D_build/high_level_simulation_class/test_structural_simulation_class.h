#include "solid_structural_simulation_class.h"

Real tolerance = 1e-6;

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
	vector<array<int, 2>> Get_contacting_body_pairs_list_(){ return contacting_body_pairs_list_; };
	vector<pair<array<int, 2>, array<Real, 2>>> Get_time_dep_contacting_body_pairs_list_(){ return time_dep_contacting_body_pairs_list_; };
	vector<shared_ptr<SolidBodyRelationContact>> Get_contact_list_(){ return contact_list_; };
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
    // for TranslateSolidBody
	vector<shared_ptr<solid_dynamics::TranslateSolidBody>> Get_translation_solid_body_(){ return translation_solid_body_; };
	vector<TranslateSolidBodyTuple> Get_translation_solid_body_tuple_(){ return translation_solid_body_tuple_; };
};