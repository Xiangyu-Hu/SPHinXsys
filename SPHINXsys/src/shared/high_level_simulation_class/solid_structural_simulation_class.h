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

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, std::string body_name, TriangleMeshShape* triangle_mesh_shape, Real resolution);
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
	std::vector<std::string>* imported_stl_list;
	Real scale_stl;
	std::vector<Vec3d>* translation_list;
	Real default_resolution;
	std::vector<Real>* resolution_list;
	std::vector<LinearElasticSolid*>* material_model_list;
	Real physical_viscosity;
};

class SolidStructuralSimulation
	{
	private:
		// input members
		std::string relative_input_path_;
		std::vector<std::string>* imported_stl_list_;
		Real scale_stl_;
		std::vector<Vec3d>* translation_list_;
		Real default_resolution_;
		std::vector<Real>* resolution_list_;
		std::vector<LinearElasticSolid*>* material_model_list_;
		Real physical_viscosity_;

		// internal members
		SPHSystem* system_;
		In_Output* in_output_;

		std::vector<TriangleMeshShape*> body_mesh_list_;

		std::vector<ImportedModel*> imported_model_list_;
		std::vector<ElasticSolidParticles*> imported_model_particles_list_;
		std::vector<InnerBodyRelation*> imported_model_inner_list_;

		std::vector<solid_dynamics::CorrectConfiguration*> correct_configuration_list_;
		std::vector<solid_dynamics::StressRelaxationFirstHalf*> stress_relaxation_first_half_list_;
		std::vector<solid_dynamics::StressRelaxationSecondHalf*> stress_relaxation_second_half_list_;
		std::vector<DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>*> damping_list_;

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
		void ImportSTLModels();
		BoundingBox* CalculateSystemBoundaries(); //for SetupSystem
		void SetupSystem();
		void InitializeElasticBodies();
		void InitializeContactBetweenTwoBodies(int first, int second);

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
 		SolidStructuralSimulation(SolidStructuralSimulationInput* input)
		{
			relative_input_path_ = input->relative_input_path;
			imported_stl_list_ = input->imported_stl_list;
			scale_stl_ = input->scale_stl;
			translation_list_ = input->translation_list;
			default_resolution_ = input->default_resolution;
			resolution_list_ = input->resolution_list;
			material_model_list_ = input->material_model_list;
			physical_viscosity_ = input->physical_viscosity;
		};
		virtual ~SolidStructuralSimulation() {};
		
		// get data from private members
		TriangleMeshShape* GetBodyMesh(int body_index) { return body_mesh_list_[body_index]; };

		// boundary conditions for user
		void AddGravity(int body_index, Vec3d* gravity);
		void AddAccelerationForBodyPartInBoundingBox(int body_index, BoundingBox* bounding_box, Vec3d acceleration);
		void AddSpringDamperConstraintParticleWise(int body_index, Vec3d stiffness, Real damping_ratio);
		void AddConstrainSolidBodyRegion(int body_index);

		// high level functions for user
		void PreprocessSimulation()
		{
			ImportSTLModels();
			SetupSystem();
			InitializeElasticBodies();
			InitializeContactBetweenTwoBodies(0, 1);
		};
		void InitializeBoundaryConditions()
		{
			InitializeGravity();
			InitializeAccelerationForBodyPartInBoundingBox();
			InitializeSpringDamperConstraintParticleWise();
			InitializeConstrainSolidBodyRegion();
		};
		void RunSimulation(Real end_time);
	};

#endif //SOLID_STRUCTURAL_SIMULATION_CLASS_H