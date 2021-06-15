/**
* @file 	solid_structural_simulation_class.h
* @brief 	solid structural simulation class definition
* @details	solid structural simulation class for general structural simulations
* @author 	Bence Z. Rochlitz
*/

#ifndef SOLID_STRUCTURAL_SIMULATION_CLASS_H
#define SOLID_STRUCTURAL_SIMULATION_CLASS_H

#include "sphinxsys.h"
#include <algorithm>

using namespace SPH;
using namespace std;
using IndexPair = pair<int, int>;
using GravityPair = pair<int, Vec3d>;
using AccelTuple = tuple<int, BoundingBox, Vec3d>;
using SpringDamperTuple = tuple<int, Vec3d, Real>;

class BodyPartByParticleTriMesh : public BodyPartByParticle
{
public:
	BodyPartByParticleTriMesh(SPHBody* body, string body_part_name, TriangleMeshShape* triangle_mesh_shape);
	~BodyPartByParticleTriMesh();
};

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, string body_name, TriangleMeshShape* triangle_mesh_shape, Real resolution);
	~ImportedModel();
};

class SolidBodyForSimulation
{
private:
	ImportedModel imported_model_;
	//LinearElasticSolid material_model_;
	ElasticSolidParticles elastic_solid_particles_;
	InnerBodyRelation inner_body_relation_;

	solid_dynamics::CorrectConfiguration correct_configuration_;
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half_;
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half_;
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>> damping_random_;

public:
	SolidBodyForSimulation(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, Real resolution, Real physical_viscosity, LinearElasticSolid& material_model);
	~SolidBodyForSimulation(){};

	ImportedModel* GetImportedModel() { return &imported_model_; };
	//LinearElasticSolid* GetMaterialModel() { return &material_model_; };
	ElasticSolidParticles* GetElasticSolidParticles() { return &elastic_solid_particles_; };
	InnerBodyRelation* GetInnerBodyRelation() { return &inner_body_relation_; };

	solid_dynamics::CorrectConfiguration* GetCorrectConfiguration() { return &correct_configuration_; };
	solid_dynamics::StressRelaxationFirstHalf* GetStressRelaxationFirstHalf() { return &stress_relaxation_first_half_; };
	solid_dynamics::StressRelaxationSecondHalf* GetStressRelaxationSecondHalf() { return &stress_relaxation_second_half_; };
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>* GetDampingWithRandomChoice() { return &damping_random_; };
};

void ExpandBoundingBox(BoundingBox* original, BoundingBox* additional);

void RelaxParticlesSingleResolution(In_Output* in_output,
									bool write_particles_to_file,
									ImportedModel* imported_model,
									ElasticSolidParticles* imported_model_particles,
									InnerBodyRelation* imported_model_inner);


class StructuralSimulationInput
{
public:
	string relative_input_path_;
	vector<string> imported_stl_list_;
	Real scale_stl_;
	vector<Vec3d> translation_list_;
	vector<Real> resolution_list_;
	vector<LinearElasticSolid> material_model_list_;
	Real physical_viscosity_;
	vector<IndexPair> contacting_bodies_list_;
	// boundary conditions
	vector<GravityPair> non_zero_gravity_;
	vector<AccelTuple> acceleration_bounding_box_tuple_;
	vector<SpringDamperTuple> spring_damper_tuple_;
	vector<int> body_indeces_fixed_constraint_;

	StructuralSimulationInput(
		string relative_input_path,
		vector<string> imported_stl_list,
		Real scale_stl,
		vector<Vec3d> translation_list,
		vector<Real> resolution_list,
		vector<LinearElasticSolid> material_model_list,
		Real physical_viscosity,
		vector<IndexPair> contacting_bodies_list
		);
};

class StructuralSimulation
	{
	private:
		// input members
		string relative_input_path_;
		vector<string> imported_stl_list_;
		Real scale_stl_;
		vector<Vec3d> translation_list_;
		Real default_resolution_;
		vector<Real> resolution_list_;
		vector<LinearElasticSolid> material_model_list_;
		Real physical_viscosity_;

		// internal members
		SPHSystem system_;
		In_Output in_output_;

		vector<TriangleMeshShape> body_mesh_list_;
		vector<TriangleMeshShape> primitive_shape_list_;

		vector<SolidBodyForSimulation*> solid_body_list_;

		vector<IndexPair> contacting_bodies_list_;
		vector<SolidContactBodyRelation*> contact_list_;
		vector<solid_dynamics::ContactDensitySummation*> contact_density_list_;
		vector<solid_dynamics::ContactForce*> contact_force_list_;

		// for InitializeGravity
		vector<InitializeATimeStep*> initialize_gravity_;
		vector<GravityPair> non_zero_gravity_;
		// for AddAccelerationForBodyPartInBoundingBox
		vector<solid_dynamics::AccelerationForBodyPartInBoundingBox*> acceleration_bounding_box_;
		vector<AccelTuple> acceleration_bounding_box_tuple_;
		// for AddSpringDamperConstraintParticleWise
		vector<solid_dynamics::SpringDamperConstraintParticleWise*> spring_damper_constraint_;
		vector<SpringDamperTuple> spring_damper_tuple_;
		// for AddSpringDamperConstraintParticleWise
		vector<solid_dynamics::ConstrainSolidBodyRegion*> fixed_constraint_;
		vector<int> body_indeces_fixed_constraint_;
		
		// for constructor, the order is important
		void ScaleTranslationAndResolution();
		void CreateBodyMeshList();
		void CalculateSystemBoundaries();
		void InitializeElasticSolidBodies();
		void InitializeContactBetweenTwoBodies(int first, int second);
		void InitializeAllContacts();

		// for InitializeBoundaryConditions
		void InitializeGravity();
		void InitializeAccelerationForBodyPartInBoundingBox();
		void InitializeSpringDamperConstraintParticleWise();
		void InitializeConstrainSolidBodyRegion();

		// for RunSimulation, the order is important
		void ExecuteCorrectConfiguration();
		void ExecuteInitializeATimeStep();
		void ExecuteAccelerationForBodyPartInBoundingBox();
		void ExecuteSpringDamperConstraintParticleWise();
		void ExecuteContactDensitySummation();
		void ExecuteContactForce();
		void ExecuteStressRelaxationFirstHalf(Real dt);
		void ExecuteConstrainSolidBodyRegion();
		void ExecuteDamping(Real dt);
		void ExecuteStressRelaxationSecondHalf(Real dt);
		void ExecuteUpdateCellLinkedList();
		void ExecuteContactUpdateConfiguration();
		void RunSimulationStep(int &ite, Real &dt, Real &integration_time);

	public:
		StructuralSimulation(StructuralSimulationInput* input);
 		~StructuralSimulation();

		void RunSimulation(Real end_time);
	};

#endif //SOLID_STRUCTURAL_SIMULATION_CLASS_H