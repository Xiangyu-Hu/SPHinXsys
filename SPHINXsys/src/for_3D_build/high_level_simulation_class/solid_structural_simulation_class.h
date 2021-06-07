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

class BodyPartByParticleTriMesh : public BodyPartByParticle
{
public:
	BodyPartByParticleTriMesh(SPHBody* body, std::string body_part_name, TriangleMeshShape* triangle_mesh_shape);
	virtual ~BodyPartByParticleTriMesh() {};
};

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, std::string body_name, TriangleMeshShape* triangle_mesh_shape, Real resolution);
	~ImportedModel(){};
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
	SolidBodyForSimulation(SPHSystem &system, std::string body_name, TriangleMeshShape& triangle_mesh_shape, Real resolution, Real physical_viscosity, LinearElasticSolid& material_model);
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


struct SolidStructuralSimulationInput
{
	std::string relative_input_path;
	std::vector<std::string> imported_stl_list;
	Real scale_stl;
	std::vector<Vec3d> translation_list;
	Real default_resolution;
	std::vector<Real> resolution_list;
	std::vector<LinearElasticSolid> material_model_list;
	Real physical_viscosity;
};

class SolidStructuralSimulation
	{
	private:
		// input members
		std::string relative_input_path_;
		std::vector<std::string> imported_stl_list_;
		Real scale_stl_;
		std::vector<Vec3d> translation_list_;
		Real default_resolution_;
		std::vector<Real> resolution_list_;
		Real physical_viscosity_;

		// internal members
		SPHSystem system_;
		In_Output in_output_;

		std::vector<TriangleMeshShape> body_mesh_list_;
		std::vector<TriangleMeshShape> primitive_shape_list_;

		std::vector<SolidBodyForSimulation*> solid_body_list_;

		std::vector<std::pair<int, int>> contacting_bodies_list_;
		std::vector<SolidContactBodyRelation*> contact_list_;
		std::vector<solid_dynamics::ContactDensitySummation*> contact_density_list_;
		std::vector<solid_dynamics::ContactForce*> contact_force_list_;

		// for InitializeGravity
		std::vector<InitializeATimeStep*> initialize_gravity_;
		std::vector<int> body_indeces_gravity_;
		std::vector<Vec3d*> gravity_;
		// for AddAccelerationForBodyPartInBoundingBox
		std::vector<solid_dynamics::AccelerationForBodyPartInBoundingBox*> acceleration_for_body_part_;
		std::vector<int> body_indeces_accelerations_;
		std::vector<BoundingBox*> bounding_boxes_;
		std::vector<Vec3d> accelerations_;
		// for AddSpringDamperConstraintParticleWise
		std::vector<solid_dynamics::SpringDamperConstraintParticleWise*> spring_damper_contraint_;
		std::vector<int> body_indeces_spring_damper_;
		std::vector<Vec3d> stiffnesses_;
		std::vector<Real> damping_ratios_;
		// for AddSpringDamperConstraintParticleWise
		std::vector<solid_dynamics::ConstrainSolidBodyRegion*> fixed_contraint_;
		std::vector<int> body_indeces_fixed_contraint_;
		
		// for PreprocessSimulation, the order is important
		void ImportSTLModelsAndAddPrimitives();
		void CalculateSystemBoundaries(); //for SetupSystem
		void SetupSystem();
		void InitializeElasticBodies(bool write_particle_relaxation);
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
		SolidStructuralSimulation(SolidStructuralSimulationInput* input);
 		~SolidStructuralSimulation();

		//add primitive shapes
		void AddPrimitiveCuboid(Vec3d halfsize_cuboid, Vec3d translation, Real resolution, LinearElasticSolid& material);
		
		// get data from private members
		TriangleMeshShape* GetBodyMesh(int body_index) { return &body_mesh_list_[body_index]; };

		// add contacting bodies
		void AddContactPair(int first_id, int second_id);

		// boundary conditions for user
		void AddGravity(int body_index, Vec3d* gravity);
		void AddAccelerationForBodyPartInBoundingBox(int body_index, BoundingBox* bounding_box, Vec3d acceleration);
		void AddSpringDamperConstraintParticleWise(int body_index, Vec3d stiffness, Real damping_ratio);
		void AddConstrainSolidBodyRegion(int body_index);

		void InitializeBoundaryConditions()
		{
			InitializeAllContacts();
			InitializeGravity();
			InitializeAccelerationForBodyPartInBoundingBox();
			InitializeSpringDamperConstraintParticleWise();
			InitializeConstrainSolidBodyRegion();
		};
		void RunSimulation(Real end_time);
	};

#endif //SOLID_STRUCTURAL_SIMULATION_CLASS_H