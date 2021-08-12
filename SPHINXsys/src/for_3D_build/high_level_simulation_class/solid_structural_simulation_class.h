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
#include <memory>
#include <vector>

using namespace SPH;
using namespace std;
using GravityPair = pair<int, Vec3d>;
using AccelTuple = tuple<int, BoundingBox, Vec3d>;
using ForceTuple = tuple<int, BoundingBox, Vec3d, Real>;
using PressureTuple = tuple<int, Vec3d, StdVec<array<Real, 2>>>;
using SpringDamperTuple = tuple<int, Vec3d, Real>;
using ConstrainedRegionPair = pair<int, BoundingBox>;
using PositionSolidBodyTuple = tuple<int, Real, Real, Vec3d>;
using PositionScaleSolidBodyTuple = tuple<int, Real, Real, Real>;
using TranslateSolidBodyTuple = tuple<int, Real, Real, Vec3d>;
using TranslateSolidBodyPartTuple = tuple<int, Real, Real, Vec3d, BoundingBox>;

class BodyPartByParticleTriMesh : public BodyPartByParticle
{
public:
	BodyPartByParticleTriMesh(SPHBody* body, string body_part_name, TriangleMeshShape* triangle_mesh_shape);
	~BodyPartByParticleTriMesh();
};

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, string body_name, TriangleMeshShape* triangle_mesh_shape, ParticleAdaptation* particle_adaptation);
	~ImportedModel();
};

class SolidBodyForSimulation
{
private:
	ImportedModel imported_model_;
	//LinearElasticSolid material_model_;
	ElasticSolidParticles elastic_solid_particles_;
	BodyRelationInner inner_body_relation_;

	solid_dynamics::CorrectConfiguration correct_configuration_;
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half_;
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half_;
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>> damping_random_;

public:
	SolidBodyForSimulation(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, ParticleAdaptation& particle_adaptation, Real physical_viscosity, LinearElasticSolid& material_model);
	~SolidBodyForSimulation(){};

	ImportedModel* getImportedModel() { return &imported_model_; };
	//LinearElasticSolid* GetMaterialModel() { return &material_model_; };
	ElasticSolidParticles* getElasticSolidParticles() { return &elastic_solid_particles_; };
	BodyRelationInner* getInnerBodyRelation() { return &inner_body_relation_; };

	solid_dynamics::CorrectConfiguration* getCorrectConfiguration() { return &correct_configuration_; };
	solid_dynamics::StressRelaxationFirstHalf* getStressRelaxationFirstHalf() { return &stress_relaxation_first_half_; };
	solid_dynamics::StressRelaxationSecondHalf* getStressRelaxationSecondHalf() { return &stress_relaxation_second_half_; };
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>* getDampingWithRandomChoice() { return &damping_random_; };
};

void expandBoundingBox(BoundingBox* original, BoundingBox* additional);

void relaxParticlesSingleResolution(In_Output* in_output,
									bool write_particles_to_file,
									ImportedModel* imported_model,
									ElasticSolidParticles* imported_model_particles,
									BodyRelationInner* imported_model_inner);


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
	vector<array<int, 2>> contacting_body_pairs_list_;
	vector<pair<array<int, 2>, array<Real, 2>>> time_dep_contacting_body_pairs_list_;
	// scale system boundaries
	Real scale_system_boundaries_;
	// particle relaxation
	vector<bool> particle_relaxation_list_;
	bool write_particle_relaxation_data_;
	// boundary conditions
	vector<GravityPair> non_zero_gravity_;
	vector<AccelTuple> acceleration_bounding_box_tuple_;
	vector<ForceTuple> force_in_body_region_tuple_;
	vector<PressureTuple> surface_pressure_tuple_;
	vector<SpringDamperTuple> spring_damper_tuple_;
	vector<int> body_indeces_fixed_constraint_;
	vector<ConstrainedRegionPair> body_indeces_fixed_constraint_region_;
	vector<PositionSolidBodyTuple> position_solid_body_tuple_;
	vector<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
	vector<TranslateSolidBodyTuple> translation_solid_body_tuple_;
	vector<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

	StructuralSimulationInput(
		string relative_input_path,
		vector<string> imported_stl_list,
		Real scale_stl,
		vector<Vec3d> translation_list,
		vector<Real> resolution_list,
		vector<LinearElasticSolid> material_model_list,
		Real physical_viscosity,
		vector<array<int, 2>> contacting_bodies_list
	);
};

class StructuralSimulation
	{
	protected:
		// mandatory input
		string relative_input_path_;
		vector<string> imported_stl_list_;
		Real scale_stl_;
		vector<Vec3d> translation_list_;
		vector<Real> resolution_list_;
		vector<LinearElasticSolid> material_model_list_;
		Real physical_viscosity_;
		vector<array<int, 2>> contacting_body_pairs_list_;
		vector<pair<array<int, 2>, array<Real, 2>>> time_dep_contacting_body_pairs_list_; //optional: time dependent contact
		vector<bool> particle_relaxation_list_; // optional: particle relaxation
		bool write_particle_relaxation_data_;

		// internal members
		Real system_resolution_;
		SPHSystem system_;
		Real scale_system_boundaries_;
		In_Output in_output_;

		vector<TriangleMeshShape> body_mesh_list_;
		vector<ParticleAdaptation> particle_adaptation_list_;
		vector<shared_ptr<SolidBodyForSimulation>> solid_body_list_;
		vector<shared_ptr<solid_dynamics::UpdateElasticNormalDirection>> particle_normal_update_;

		vector<shared_ptr<SolidBodyRelationContact>> contact_list_;
		vector<shared_ptr<solid_dynamics::ContactDensitySummation>> contact_density_list_;
		vector<shared_ptr<solid_dynamics::ContactForce>> contact_force_list_;

		// for initializeATimeStep
		vector<shared_ptr<TimeStepInitialization>> initialize_gravity_;
		vector<GravityPair> non_zero_gravity_;
		// for AccelerationForBodyPartInBoundingBox
		vector<shared_ptr<solid_dynamics::AccelerationForBodyPartInBoundingBox>> acceleration_bounding_box_;
		vector<AccelTuple> acceleration_bounding_box_tuple_;
		// for ForceInBodyRegion
		vector<shared_ptr<solid_dynamics::ForceInBodyRegion>> force_in_body_region_;
		vector<ForceTuple> force_in_body_region_tuple_;
		// for SurfacePressureFromSource
		vector<shared_ptr<solid_dynamics::SurfacePressureFromSource>> surface_pressure_;
		vector<PressureTuple> surface_pressure_tuple_;
		// for SpringDamperConstraintParticleWise
		vector<shared_ptr<solid_dynamics::SpringDamperConstraintParticleWise>> spring_damper_constraint_;
		vector<SpringDamperTuple> spring_damper_tuple_;
		// for ConstrainSolidBody
		vector<shared_ptr<solid_dynamics::ConstrainSolidBodyRegion>> fixed_constraint_body_;
		vector<int> body_indeces_fixed_constraint_;
		// for ConstrainSolidBodyRegion
		vector<shared_ptr<solid_dynamics::ConstrainSolidBodyRegion>> fixed_constraint_region_;
		vector<ConstrainedRegionPair> body_indeces_fixed_constraint_region_;
		// for PositionSolidBody
		vector<shared_ptr<solid_dynamics::PositionSolidBody>> position_solid_body_;
		vector<PositionSolidBodyTuple> position_solid_body_tuple_;
		// for PositionScaleSolidBody
		vector<shared_ptr<solid_dynamics::PositionScaleSolidBody>> position_scale_solid_body_;
		vector<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
		// for TranslateSolidBody
		vector<shared_ptr<solid_dynamics::TranslateSolidBody>> translation_solid_body_;
		vector<TranslateSolidBodyTuple> translation_solid_body_tuple_;
		// for TranslateSolidBodyPart
		vector<shared_ptr<solid_dynamics::TranslateSolidBodyPart>> translation_solid_body_part_;
		vector<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

		// iterators
		int iteration_;

		// data storage
		vector<Real> von_mises_stress_max_;
		
		// for constructor, the order is important
		void scaleTranslationAndResolution();
		void setSystemResolutionMax();
		void createBodyMeshList();
		void createParticleAdaptationList();
		void calculateSystemBoundaries();
		void initializeElasticSolidBodies();
		void initializeContactBetweenTwoBodies(int first, int second);
		void initializeAllContacts();

		// for initializeBoundaryConditions
		void initializeGravity();
		void initializeAccelerationForBodyPartInBoundingBox();
		void initializeForceInBodyRegion();
		void initializeSurfacePressure();
		void initializeSpringDamperConstraintParticleWise();
		void initializeConstrainSolidBody();
		void initializeConstrainSolidBodyRegion();
		void initializePositionSolidBody();
		void initializePositionScaleSolidBody();
		void initializeTranslateSolidBody();
		void initializeTranslateSolidBodyPart();

		// for runSimulation, the order is important
		void executeCorrectConfiguration();
		void executeUpdateElasticNormalDirection();
		void executeinitializeATimeStep();
		void executeAccelerationForBodyPartInBoundingBox();
		void executeForceInBodyRegion();
		void executeSurfacePressure();
		void executeSpringDamperConstraintParticleWise();
		void executeContactDensitySummation();
		void executeContactForce();
		void executeStressRelaxationFirstHalf(Real dt);
		void executeConstrainSolidBody();
		void executeConstrainSolidBodyRegion();
		void executePositionSolidBody(Real dt);
		void executePositionScaleSolidBody(Real dt);
		void executeTranslateSolidBody(Real dt);
		void executeTranslateSolidBodyPart(Real dt);
		void executeDamping(Real dt);
		void executeStressRelaxationSecondHalf(Real dt);
		void executeUpdateCellLinkedList();
		void executeContactUpdateConfiguration();

		void initializeSimulation();

		void runSimulationStep(Real &dt, Real &integration_time);

	public:
		StructuralSimulation(StructuralSimulationInput& input);
		~StructuralSimulation();

		//For c++
		void runSimulation(Real end_time);
	
		//For JS
		double runSimulationFixedDurationJS(int number_of_steps);
	};

#endif //SOLID_STRUCTURAL_SIMULATION_CLASS_H