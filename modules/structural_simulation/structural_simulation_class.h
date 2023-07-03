/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	structural_simulation_class.h
 * @brief 	The structural simulation module is licensed under the Aladdin Free Public License
 * (https://spdx.org/licenses/Aladdin.html) regarding usage for medical device development.
 * Commercial use for medical device development is not permitted. This does not apply to applications in other fields.
 * @details	solid structural simulation class for general structural simulations
 * @author 	Bence Z. Rochlitz - Virtonomy GmbH, Xiangyu Hu
 */

#ifndef SOLID_STRUCTURAL_SIMULATION_CLASS_H
#define SOLID_STRUCTURAL_SIMULATION_CLASS_H

#include "sphinxsys.h"

using namespace SPH;

using GravityPair = std::pair<int, Vec3d>;
using AccelTuple = std::tuple<int, BoundingBox, Vec3d>;
using ForceTuple = std::tuple<int, BoundingBox, Vec3d, Real>;
using PressureTuple = std::tuple<int, SharedPtr<TriangleMeshShape>, Vec3d, StdVec<std::array<Real, 2>>>;
using SpringDamperTuple = std::tuple<int, Vec3d, Real>;
/**
 * @brief SurfaceSpringTuple
 * int: body index
 * TriangleMeshShape*: the body part, the normal spring is applied to
 * bool: if true, the "outer" surface is considered (particle normals > 90° from the particle-source point vector), if false, the "inner" surface
 * Vec3d: source point to relate inner and outer surface
 * Real: normal spring stiffness
 * Real: damping coefficient
 */
using SurfaceSpringTuple = std::tuple<int, SharedPtr<TriangleMeshShape>, bool, Vec3d, Real, Real>;
using ConstrainedRegionPair = std::pair<int, BoundingBox>;
using PositionSolidBodyTuple = std::tuple<int, Real, Real, Vec3d>;
using PositionScaleSolidBodyTuple = std::tuple<int, Real, Real, Real>;
using TranslateSolidBodyTuple = std::tuple<int, Real, Real, Vec3d>;
using TranslateSolidBodyPartTuple = std::tuple<int, Real, Real, Vec3d, BoundingBox>;

#ifdef __EMSCRIPTEN__
struct StlData
{
    string name;
    uintptr_t ptr;
};

using StlList = StdVec<StlData>;
#else
using StlList = StdVec<std::string>;
#endif

class BodyPartFromMesh : public BodyRegionByParticle
{
  public:
    BodyPartFromMesh(SPHBody &body, SharedPtr<TriangleMeshShape> triangle_mesh_shape_ptr);
    ~BodyPartFromMesh(){};
};

class SolidBodyFromMesh : public SolidBody
{
  public:
    SolidBodyFromMesh(SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
                      SharedPtr<SaintVenantKirchhoffSolid> material_model, StdLargeVec<Vec3d> &pos_0, StdLargeVec<Real> &volume);
    ~SolidBodyFromMesh(){};
};

class SolidBodyForSimulation
{
  private:
    SolidBodyFromMesh solid_body_from_mesh_;
    InnerRelation inner_body_relation_;

    SimpleDynamics<NormalDirectionFromBodyShape> initial_normal_direction_;
    InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration_;
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half_;
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>> damping_random_;

  public:
    // no particle reload --> direct generator
    SolidBodyForSimulation(
        SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
        Real physical_viscosity, SharedPtr<SaintVenantKirchhoffSolid> material_model, StdLargeVec<Vec3d> &pos_0, StdLargeVec<Real> &volume);
    ~SolidBodyForSimulation(){};

    SolidBodyFromMesh *getSolidBodyFromMesh() { return &solid_body_from_mesh_; };
    ElasticSolidParticles *getElasticSolidParticles() { return DynamicCast<ElasticSolidParticles>(this, &solid_body_from_mesh_.getBaseParticles()); };
    InnerRelation *getInnerBodyRelation() { return &inner_body_relation_; };

    SimpleDynamics<NormalDirectionFromBodyShape> *getInitialNormalDirection() { return &initial_normal_direction_; };
    InteractionWithUpdate<CorrectedConfigurationInner> *getCorrectConfiguration() { return &correct_configuration_; };
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> *getStressRelaxationFirstHalf() { return &stress_relaxation_first_half_; };
    Dynamics1Level<solid_dynamics::Integration2ndHalf> *getStressRelaxationSecondHalf() { return &stress_relaxation_second_half_; };
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>> *getDampingWithRandomChoice() { return &damping_random_; };
};

void expandBoundingBox(BoundingBox *original, BoundingBox *additional);

void relaxParticlesSingleResolution(IOEnvironment &io_environment,
                                    bool write_particles_to_file,
                                    SolidBodyFromMesh &solid_body_from_mesh,
                                    ElasticSolidParticles &solid_body_from_mesh_particles,
                                    InnerRelation &solid_body_from_mesh_inner);

static inline Real getPhysicalViscosityGeneral(Real rho, Real youngs_modulus, Real length_scale, Real shape_constant = 1.0)
{
    // the physical viscosity is defined in the paper pf prof. Hu
    // https://arxiv.org/pdf/2103.08932.pdf
    // physical_viscosity = (beta / 4) * sqrt(rho * Young's modulus) * L
    // beta: shape constant --> how to define it? - it's 1 for now.. TODO
    // L: length scale of the problem --> 10 mm roughly
    return shape_constant / 4.0 * std::sqrt(rho * youngs_modulus) * length_scale;
}

class StructuralSimulationInput
{
  public:
    std::string relative_input_path_;
    StlList imported_stl_list_;
    Real scale_stl_;
    StdVec<Vec3d> translation_list_;
    StdVec<Real> resolution_list_;
    StdVec<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list_;
    StdVec<Real> physical_viscosity_;
    StdVec<IndexVector> contacting_body_pairs_list_;
    StdVec<std::pair<std::array<int, 2>, std::array<Real, 2>>> time_dep_contacting_body_pairs_list_;
    // scale system boundaries
    Real scale_system_boundaries_;
    // particle relaxation
    StdVec<bool> particle_relaxation_list_;
    bool write_particle_relaxation_data_;
    // boundary conditions
    StdVec<GravityPair> non_zero_gravity_;
    StdVec<AccelTuple> acceleration_bounding_box_tuple_;
    StdVec<ForceTuple> force_in_body_region_tuple_;
    StdVec<PressureTuple> surface_pressure_tuple_;
    StdVec<SpringDamperTuple> spring_damper_tuple_;
    StdVec<SurfaceSpringTuple> surface_spring_tuple_;
    StdVec<int> body_indices_fixed_constraint_;
    StdVec<ConstrainedRegionPair> body_indices_fixed_constraint_region_;
    StdVec<PositionSolidBodyTuple> position_solid_body_tuple_;
    StdVec<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
    StdVec<TranslateSolidBodyTuple> translation_solid_body_tuple_;
    StdVec<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

    StructuralSimulationInput(
        std::string relative_input_path,
        StlList imported_stl_list,
        Real scale_stl,
        StdVec<Vec3d> translation_list,
        StdVec<Real> resolution_list,
        StdVec<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list,
        StdVec<Real> physical_viscosity,
        StdVec<IndexVector> contacting_bodies_list);
};

class StructuralSimulation
{
  private:
    UniquePtrsKeeper<SurfaceContactRelation> contact_relation_ptr_keeper_;
    UniquePtrsKeeper<Gravity> gravity_ptr_keeper_;
    UniquePtrsKeeper<BodyPartFromMesh> body_part_tri_mesh_ptr_keeper_;

  protected:
    // mandatory input
    std::string relative_input_path_;
    StlList imported_stl_list_;
    Real scale_stl_;
    StdVec<Vec3d> translation_list_;
    StdVec<Real> resolution_list_;
    StdVec<SharedPtr<TriangleMeshShape>> body_mesh_list_;
    StdVec<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list_;
    StdVec<Real> physical_viscosity_;
    StdVec<IndexVector> contacting_body_pairs_list_;
    StdVec<std::pair<std::array<int, 2>, std::array<Real, 2>>> time_dep_contacting_body_pairs_list_; // optional: time dependent contact
    StdVec<bool> particle_relaxation_list_;                                                          // optional: particle relaxation
    bool write_particle_relaxation_data_;

    // internal members
    Real system_resolution_;
    SPHSystem system_;
    Real scale_system_boundaries_;
    IOEnvironment io_environment_;

    StdVec<SharedPtr<SolidBodyForSimulation>> solid_body_list_;
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection>>> particle_normal_update_;

    StdVec<SharedPtr<SurfaceContactRelation>> contact_list_;
    StdVec<SharedPtr<InteractionDynamics<solid_dynamics::ContactDensitySummation>>> contact_density_list_;
    StdVec<SharedPtr<InteractionDynamics<solid_dynamics::ContactForce>>> contact_force_list_;

    // for initializeATimeStep
    StdVec<SharedPtr<SimpleDynamics<TimeStepInitialization>>> initialize_time_step_;
    StdVec<GravityPair> non_zero_gravity_;
    // for AccelerationForBodyPartInBoundingBox
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::AccelerationForBodyPartInBoundingBox>>> acceleration_bounding_box_;
    StdVec<AccelTuple> acceleration_bounding_box_tuple_;
    // for ForceInBodyRegion
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::ForceInBodyRegion>>> force_in_body_region_;
    StdVec<ForceTuple> force_in_body_region_tuple_;
    // for SurfacePressureFromSource
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::SurfacePressureFromSource>>> surface_pressure_;
    StdVec<PressureTuple> surface_pressure_tuple_;
    // for SpringDamperConstraintParticleWise
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::SpringDamperConstraintParticleWise>>> spring_damper_constraint_;
    StdVec<SpringDamperTuple> spring_damper_tuple_;
    // for SpringNormalOnSurfaceParticles
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::SpringNormalOnSurfaceParticles>>> surface_spring_;
    StdVec<SurfaceSpringTuple> surface_spring_tuple_;
    // for ConstrainSolidBody
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::FixBodyConstraint>>> fixed_constraint_body_;
    StdVec<int> body_indices_fixed_constraint_;
    // for ConstrainSolidBodyRegion
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::FixBodyPartConstraint>>> fixed_constraint_region_;
    StdVec<ConstrainedRegionPair> body_indices_fixed_constraint_region_;
    // for PositionSolidBody
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::PositionSolidBody>>> position_solid_body_;
    StdVec<PositionSolidBodyTuple> position_solid_body_tuple_;
    // for PositionScaleSolidBody
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::PositionScaleSolidBody>>> position_scale_solid_body_;
    StdVec<PositionScaleSolidBodyTuple> position_scale_solid_body_tuple_;
    // for TranslateSolidBody
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::TranslateSolidBody>>> translation_solid_body_;
    StdVec<TranslateSolidBodyTuple> translation_solid_body_tuple_;
    // for TranslateSolidBodyPart
    StdVec<SharedPtr<SimpleDynamics<solid_dynamics::TranslateSolidBodyPart>>> translation_solid_body_part_;
    StdVec<TranslateSolidBodyPartTuple> translation_solid_body_part_tuple_;

    // iterators
    int iteration_;

    // data storage
    StdVec<Real> von_mises_stress_max_;
    StdLargeVec<StdLargeVec<Real>> von_mises_stress_particles_;

    StdVec<Real> von_mises_strain_max_;
    StdLargeVec<StdLargeVec<Real>> von_mises_strain_particles_;

    // for constructor, the order is important
    void scaleTranslationAndResolution();
    void setSystemResolutionMax();
    void createBodyMeshList();
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
    void initializeSpringNormalOnSurfaceParticles();
    void initializeConstrainSolidBody();
    void initializeConstrainSolidBodyRegion();
    void initializePositionSolidBody();
    void initializePositionScaleSolidBody();
    void initializeTranslateSolidBody();
    void initializeTranslateSolidBodyPart();

    // for runSimulation, the order is important
    void executeInitialNormalDirection();
    void executeCorrectConfiguration();
    void executeUpdateElasticNormalDirection();
    void executeInitializeATimeStep();
    void executeAccelerationForBodyPartInBoundingBox();
    void executeForceInBodyRegion();
    void executeSurfacePressure();
    void executeSpringDamperConstraintParticleWise();
    void executeSpringNormalOnSurfaceParticles();
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
    explicit StructuralSimulation(const StructuralSimulationInput &input);
    ~StructuralSimulation();

    StdVec<SharedPtr<SolidBodyForSimulation>> get_solid_body_list_() { return solid_body_list_; };
    Real getMaxDisplacement(int body_index);

    // For c++
    void runSimulation(Real end_time);

    // For JS
    double runSimulationFixedDurationJS(int number_of_steps);
};

class StructuralSimulationJS : public StructuralSimulation
{
  public:
    StructuralSimulationJS(const StructuralSimulationInput &input);
    ~StructuralSimulationJS() = default;

    void runSimulationFixedDuration(int number_of_steps);

    VtuStringData getVtuData();

  private:
    BodyStatesRecordingToVtpString write_states_;
    Real dt;
};

#endif // SOLID_STRUCTURAL_SIMULATION_CLASS_H
