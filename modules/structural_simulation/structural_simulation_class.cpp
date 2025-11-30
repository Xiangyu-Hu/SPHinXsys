/**
 * @file 	structural_simulation_class.cpp
 * @brief 	The structural simulation module is licensed under the Aladdin Free Public License (https://spdx.org/licenses/Aladdin.html) regarding usage for medical device development.
 * Commercial use for medical device development is not permitted. This does not apply to applications in other fields.
 * @details	solid structural simulation class for general structural simulations
 * @author 	Bence Z. Rochlitz - Virtonomy GmbH, Xiangyu Hu
 */

#include "structural_simulation_class.h"

////////////////////////////////////////////////////
/* global functions in StructuralSimulation  */
////////////////////////////////////////////////////

BodyPartFromMesh::BodyPartFromMesh(SPHBody &body, SharedPtr<TriangleMeshShape> triangle_mesh_shape_ptr)
    : BodyRegionByParticle(body, triangle_mesh_shape_ptr)
{
    // set the body domain bounds because it is not set by default
    BoundingBoxd bounds = triangle_mesh_shape_ptr->getBounds();
    setBodyPartBounds(bounds);
}

SolidBodyFromMesh::SolidBodyFromMesh(
    SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
    SharedPtr<SaintVenantKirchhoffSolid> material_model, Vecd *pos_0, Real *volume)
    : SolidBody(system, triangle_mesh_shape)
{
    defineAdaptationRatios(1.15, system.ReferenceResolution() / resolution);
    defineBodyLevelSetShape()->cleanLevelSet();
    defineMaterial<SaintVenantKirchhoffSolid>(*material_model.get());
    generateParticles<BaseParticles, Lattice>();
}

SolidBodyForSimulation::SolidBodyForSimulation(
    SPHSystem &system, SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution,
    Real physical_viscosity, SharedPtr<SaintVenantKirchhoffSolid> material_model, Vecd *pos_0, Real *volume)
    : solid_body_from_mesh_(system, triangle_mesh_shape, resolution, material_model, pos_0, volume),
      inner_body_relation_(solid_body_from_mesh_),
      initial_normal_direction_(SimpleDynamics<NormalDirectionFromBodyShape>(solid_body_from_mesh_)),
      correct_configuration_(inner_body_relation_),
      stress_relaxation_first_half_(inner_body_relation_),
      stress_relaxation_second_half_(inner_body_relation_),
      damping_random_(0.2, inner_body_relation_, "Velocity", physical_viscosity)
{
    initial_normal_direction_.exec();
    std::cout << "  normal initialization done" << std::endl;
}

BoundingBoxd expandBoundingBox(const BoundingBoxd &original, const BoundingBoxd &additional)
{
    BoundingBoxd expanded = original;
    for (int i = 0; i < expanded.lower_.size(); i++)
    {
        if (additional.lower_[i] < expanded.lower_[i])
        {
            expanded.lower_[i] = additional.lower_[i];
        }
        if (additional.upper_[i] > expanded.upper_[i])
        {
            expanded.upper_[i] = additional.upper_[i];
        }
    }
    return expanded;
}

void relaxParticlesSingleResolution(bool write_particle_relaxation_data,
                                    SolidBody &solid_body_from_mesh,
                                    InnerRelation &solid_body_from_mesh_inner)
{
    BodyStatesRecordingToVtp write_solid_body_from_mesh_to_vtp(solid_body_from_mesh);

    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SimpleDynamics<relax_dynamics::RandomizeParticlePosition> random_solid_body_from_mesh_particles(solid_body_from_mesh);
    /** A  Physics relaxation step. */
    relax_dynamics::RelaxationStepLevelSetCorrectionInner relaxation_step_inner(solid_body_from_mesh_inner);
    //----------------------------------------------------------------------
    //	Particle relaxation starts here.
    //----------------------------------------------------------------------
    random_solid_body_from_mesh_particles.exec(0.25);
    relaxation_step_inner.SurfaceBounding().exec();
    if (write_particle_relaxation_data)
    {
        write_solid_body_from_mesh_to_vtp.writeToFile(0.0);
    }
    solid_body_from_mesh.updateCellLinkedList();
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < 1000)
    {
        relaxation_step_inner.exec();
        ite_p += 1;
        if (ite_p % 100 == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
            if (write_particle_relaxation_data)
            {
                write_solid_body_from_mesh_to_vtp.writeToFile(ite_p);
            }
        }
    }
    std::cout << "The physics relaxation process of the imported model finished !" << std::endl;
}

std::tuple<Vecd *, Real *> generateAndRelaxParticlesFromMesh(
    SharedPtr<TriangleMeshShape> triangle_mesh_shape, Real resolution, bool particle_relaxation, bool write_particle_relaxation_data)
{
    BoundingBoxd bb = triangle_mesh_shape->getBounds();
    SPHSystem system(bb, resolution);
    SolidBody model(system, triangle_mesh_shape);
    model.defineBodyLevelSetShape()->cleanLevelSet();
    model.defineMaterial<Solid>();
    model.generateParticles<BaseParticles, Lattice>();

    if (particle_relaxation)
    {
        InnerRelation inner_relation(model);
        relaxParticlesSingleResolution(write_particle_relaxation_data, model, inner_relation);
    }

    return std::tuple<Vecd *, Real *>(model.getBaseParticles().ParticlePositions(), model.getBaseParticles().VolumetricMeasures());
}

BodyPartByParticle *createBodyPartFromMesh(SPHBody &body, const StlList &stl_list, size_t body_index, SharedPtr<TriangleMeshShape> tmesh)
{
#ifdef __EMSCRIPTEN__
    tmesh->setName(stl_list[body_index].name);
    return new BodyPartFromMesh(body, tmesh);
#else
    tmesh->setName(stl_list[body_index]);
    return new BodyPartFromMesh(body, tmesh);
#endif
}

StructuralSimulationInput::StructuralSimulationInput(
    std::string relative_input_path,
    StlList imported_stl_list,
    Real scale_stl,
    StdVec<Vec3d> translation_list,
    StdVec<Real> resolution_list,
    StdVec<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list,
    StdVec<Real> physical_viscosity,
    StdVec<IndexVector> contacting_bodies_list)
    : relative_input_path_(relative_input_path),
      imported_stl_list_(imported_stl_list),
      scale_stl_(scale_stl),
      translation_list_(translation_list),
      resolution_list_(resolution_list),
      material_model_list_(material_model_list),
      physical_viscosity_(physical_viscosity),
      contacting_body_pairs_list_(contacting_bodies_list),
      time_dep_contacting_body_pairs_list_({}),
      particle_relaxation_list_({})
{
    for (size_t i = 0; i < resolution_list_.size(); i++)
    {
        particle_relaxation_list_.push_back(true);
    }
    write_particle_relaxation_data_ = false;
    // scale system boundaries
    scale_system_boundaries_ = 1;
    // boundary conditions
    non_zero_gravity_ = {};
    force_bounding_box_tuple_ = {};
    force_in_body_region_tuple_ = {};
    surface_pressure_tuple_ = {};
    spring_damper_tuple_ = {};
    surface_spring_tuple_ = {};
    body_indices_fixed_constraint_ = {};
    body_indices_fixed_constraint_region_ = {};
    position_solid_body_tuple_ = {};
    position_scale_solid_body_tuple_ = {};
    translation_solid_body_tuple_ = {};
    translation_solid_body_part_tuple_ = {};
};

///////////////////////////////////////
/* StructuralSimulation members */
///////////////////////////////////////

StructuralSimulation::StructuralSimulation(const StructuralSimulationInput &input)
    : // generic input
      relative_input_path_(input.relative_input_path_),
      imported_stl_list_(input.imported_stl_list_),
      scale_stl_(input.scale_stl_),
      translation_list_(input.translation_list_),
      resolution_list_(input.resolution_list_),
      material_model_list_(input.material_model_list_),
      physical_viscosity_(input.physical_viscosity_),
      contacting_body_pairs_list_(input.contacting_body_pairs_list_),
      time_dep_contacting_body_pairs_list_(input.time_dep_contacting_body_pairs_list_),

      // default system, optional: particle relaxation, scale_system_boundaries
      particle_relaxation_list_(input.particle_relaxation_list_),
      write_particle_relaxation_data_(input.write_particle_relaxation_data_),
      system_resolution_(0.0),
      system_(SPHSystem(BoundingBoxd(Vec3d::Zero(), Vec3d::Zero()), system_resolution_)),
      scale_system_boundaries_(input.scale_system_boundaries_),
      physical_time_(*system_.getSystemVariableDataByName<Real>("PhysicalTime")),

      // optional: boundary conditions
      non_zero_gravity_(input.non_zero_gravity_),
      force_bounding_box_tuple_(input.force_bounding_box_tuple_),
      force_in_body_region_tuple_(input.force_in_body_region_tuple_),
      surface_pressure_tuple_(input.surface_pressure_tuple_),
      spring_damper_tuple_(input.spring_damper_tuple_),
      surface_spring_tuple_(input.surface_spring_tuple_),
      body_indices_fixed_constraint_(input.body_indices_fixed_constraint_),
      body_indices_fixed_constraint_region_(input.body_indices_fixed_constraint_region_),
      position_solid_body_tuple_(input.position_solid_body_tuple_),
      position_scale_solid_body_tuple_(input.position_scale_solid_body_tuple_),
      translation_solid_body_tuple_(input.translation_solid_body_tuple_),
      translation_solid_body_part_tuple_(input.translation_solid_body_part_tuple_),

      // iterators
      iteration_(0)
{
    // scaling of translation and resolution
    scaleTranslationAndResolution();
    // set the default resolution to the max in the resolution list
    setSystemResolutionMax();
    // create the body mesh list for triangular mesh shapes storage
    createBodyMeshList();
    // set up the system
    calculateSystemBoundaries();
    system_.setRunParticleRelaxation(true);
    // initialize solid bodies with their properties
    initializeElasticSolidBodies();
    // contacts
    initializeAllContacts();
    // boundary conditions
    initializeGravity();
    initializeExternalForceInBoundingBox();
    initializeForceInBodyRegion();
    initializeSurfacePressure();
    initializeSpringDamperConstraintParticleWise();
    initializeSpringNormalOnSurfaceParticles();
    initializeConstrainSolidBody();
    initializeConstrainSolidBodyRegion();
    initializePositionSolidBody();
    initializePositionScaleSolidBody();
    initializeTranslateSolidBody();
    initializeTranslateSolidBodyPart();
    // initialize simulation
    initializeSimulation();
}

StructuralSimulation::~StructuralSimulation()
{
}

void StructuralSimulation::scaleTranslationAndResolution()
{
    // scale the translation_list_, system_resolution_ and resolution_list_
    for (size_t i = 0; i < translation_list_.size(); i++)
    {
        translation_list_[i] *= scale_stl_;
    }
    system_resolution_ *= scale_stl_;
    for (size_t i = 0; i < resolution_list_.size(); i++)
    {
        resolution_list_[i] *= scale_stl_;
    }
}

void StructuralSimulation::setSystemResolutionMax()
{
    system_resolution_ = 0.0;
    for (size_t i = 0; i < resolution_list_.size(); i++)
    {
        if (system_resolution_ < resolution_list_[i])
        {
            system_resolution_ = resolution_list_[i];
        }
    }
    system_.setReferenceResolution(system_resolution_);
}

void StructuralSimulation::calculateSystemBoundaries()
{
    // calculate system bounds from all bodies
    for (size_t i = 0; i < body_mesh_list_.size(); i++)
    {
        BoundingBoxd additional = body_mesh_list_[i]->getBounds();
        system_.setSystemDomainBounds(expandBoundingBox(system_.getSystemDomainBounds(), additional));
    }
    // scale the system bounds around the center point
    Vecd center_point = (system_.getSystemDomainBounds().lower_ + system_.getSystemDomainBounds().upper_) * 0.5;

    Vecd distance_first = system_.getSystemDomainBounds().lower_ - center_point;
    Vecd distance_second = system_.getSystemDomainBounds().upper_ - center_point;

    system_.getSystemDomainBounds().lower_ = center_point + distance_first * scale_system_boundaries_;
    system_.getSystemDomainBounds().upper_ = center_point + distance_second * scale_system_boundaries_;
}

void StructuralSimulation::createBodyMeshList()
{
    body_mesh_list_ = {};
    for (size_t i = 0; i < imported_stl_list_.size(); i++)
    {
        std::string relative_input_path_copy = relative_input_path_;
#ifdef __EMSCRIPTEN__
        body_mesh_list_.push_back(makeShared<TriangleMeshShapeSTL>(reinterpret_cast<const uint8_t *>(imported_stl_list_[i].ptr), translation_list_[i], scale_stl_, imported_stl_list_[i].name));
#else
        body_mesh_list_.push_back(makeShared<TriangleMeshShapeSTL>(relative_input_path_copy.append(imported_stl_list_[i]), translation_list_[i], scale_stl_, imported_stl_list_[i]));
#endif
    }
}

void StructuralSimulation::initializeElasticSolidBodies()
{
    solid_body_list_ = {};
    particle_normal_update_ = {};
    for (size_t i = 0; i < body_mesh_list_.size(); i++)
    {
        std::string temp_name = "";
#ifdef __EMSCRIPTEN__
        temp_name.append(imported_stl_list_[i].name);
#else  // __EMSCRIPTEN__
        temp_name.append(imported_stl_list_[i]);
#endif // __EMSCRIPTEN__
       // we delete the .stl ending
        temp_name.erase(temp_name.size() - 4);
        // create the initial particles from the triangle mesh shape with particle relaxation option
        std::tuple<Vecd *, Real *> particles =
            generateAndRelaxParticlesFromMesh(body_mesh_list_[i], resolution_list_[i], particle_relaxation_list_[i], write_particle_relaxation_data_);

        // get the particles' initial position and their volume
        Vecd *pos_0 = std::get<0>(particles);
        Real *volume = std::get<1>(particles);

        // create the SolidBodyForSimulation
        solid_body_list_.emplace_back(makeShared<SolidBodyForSimulation>(
            system_, body_mesh_list_[i], resolution_list_[i], physical_viscosity_[i], material_model_list_[i], pos_0, volume));

        // update normal direction of particles
        particle_normal_update_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection>>(*solid_body_list_[i]->getSolidBodyFromMesh()));
    }
}

void StructuralSimulation::initializeContactBetweenTwoBodies(int first, int second)
{
    SolidBodyFromMesh *first_body = solid_body_list_[first]->getSolidBodyFromMesh();
    SolidBodyFromMesh *second_body = solid_body_list_[second]->getSolidBodyFromMesh();

    contact_list_.emplace_back(makeShared<SurfaceContactRelation>(*first_body, RealBodyVector({second_body})));
    contact_list_.emplace_back(makeShared<SurfaceContactRelation>(*second_body, RealBodyVector({first_body})));

    int last = contact_list_.size() - 1;
    contact_density_list_.push_back(makeShared<InteractionDynamics<solid_dynamics::ContactFactorSummation>>(*contact_list_[last - 1]));
    contact_density_list_.push_back(makeShared<InteractionDynamics<solid_dynamics::ContactFactorSummation>>(*contact_list_[last]));

    contact_force_list_.push_back(makeShared<InteractionWithUpdate<solid_dynamics::ContactForce>>(*contact_list_[last - 1]));
    contact_force_list_.push_back(makeShared<InteractionWithUpdate<solid_dynamics::ContactForce>>(*contact_list_[last]));
}

void StructuralSimulation::initializeAllContacts()
{
    contact_list_ = {};
    contact_density_list_ = {};
    contact_force_list_ = {};
    // first place all the regular contacts into the lists
    for (size_t i = 0; i < contacting_body_pairs_list_.size(); i++)
    {
        SolidBodyFromMesh *contact_body = solid_body_list_[i]->getSolidBodyFromMesh();
        RealBodyVector target_list = {};

        for (size_t target_i : contacting_body_pairs_list_[i])
        {
            target_list.emplace_back(solid_body_list_[target_i]->getSolidBodyFromMesh());
        }

        contact_list_.emplace_back(makeShared<SurfaceContactRelation>(*contact_body, target_list));
        int last = contact_list_.size() - 1;
        contact_density_list_.emplace_back(makeShared<InteractionDynamics<solid_dynamics::ContactFactorSummation>>(*contact_list_[last]));
        contact_force_list_.emplace_back(makeShared<InteractionWithUpdate<solid_dynamics::ContactForce>>(*contact_list_[last]));
    }
    // continue appending the lists with the time dependent contacts
    for (size_t i = 0; i < time_dep_contacting_body_pairs_list_.size(); i++)
    {
        int body_1 = time_dep_contacting_body_pairs_list_[i].first[0];
        int body_2 = time_dep_contacting_body_pairs_list_[i].first[1];
        initializeContactBetweenTwoBodies(body_1, body_2); // vector with first element being array with indices
    }
}

void StructuralSimulation::initializeGravity()
{
    // collect all the body indices with non-zero gravity
    StdVec<int> gravity_indices = {};
    for (size_t i = 0; i < non_zero_gravity_.size(); i++)
    {
        gravity_indices.push_back(non_zero_gravity_[i].first);
    }
    // initialize gravity
    initialize_gravity_ = {};
    size_t gravity_index_i = 0; // iterating through gravity_indices
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        // check if i is in indices_gravity
        if (count(gravity_indices.begin(), gravity_indices.end(), i))
        {
            SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(non_zero_gravity_[gravity_index_i].second);
            gravity_list_.emplace_back(non_zero_gravity_[gravity_index_i].second);
            initialize_gravity_.emplace_back(makeShared<SimpleDynamics<GravityForce<Gravity>>>(
                *solid_body_list_[i]->getSolidBodyFromMesh(), gravity_list_.back()));
            gravity_index_i++;
        }
    }
}

void StructuralSimulation::initializeExternalForceInBoundingBox()
{
    force_bounding_box_ = {};
    for (size_t i = 0; i < force_bounding_box_tuple_.size(); i++)
    {
        SolidBody *solid_body = solid_body_list_[std::get<0>(force_bounding_box_tuple_[i])]->getSolidBodyFromMesh();
        force_bounding_box_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::ExternalForceInBoundingBox>>(
            *solid_body, std::get<1>(force_bounding_box_tuple_[i]), std::get<2>(force_bounding_box_tuple_[i])));
    }
}

void StructuralSimulation::initializeForceInBodyRegion()
{
    force_in_body_region_ = {};
    for (size_t i = 0; i < force_in_body_region_tuple_.size(); i++)
    {
        int body_index = std::get<0>(force_in_body_region_tuple_[i]);
        BoundingBoxd bbox = std::get<1>(force_in_body_region_tuple_[i]);
        Vec3d force = std::get<2>(force_in_body_region_tuple_[i]);
        Real end_time = std::get<3>(force_in_body_region_tuple_[i]);

        // get the length of each side to create the box
        Real x_side = bbox.upper_[0] - bbox.lower_[0];
        Real y_side = bbox.upper_[1] - bbox.lower_[1];
        Real z_side = bbox.upper_[2] - bbox.lower_[2];
        Vec3d halfsize_bbox(0.5 * x_side, 0.5 * y_side, 0.5 * z_side);
        // get the center point for translation from the origin
        Vec3d center = (bbox.upper_ + bbox.lower_) * 0.5;
        // SimTK geometric modeling resolution
        int resolution(20);
        // create the triangle mesh of the box
        BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), makeShared<TriangleMeshShapeBrick>(halfsize_bbox, resolution, center, imported_stl_list_[body_index]));
        force_in_body_region_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::ForceInBodyRegion>>(*bp, force, end_time));
    }
}

void StructuralSimulation::initializeSurfacePressure()
{
    surface_pressure_ = {};
    for (size_t i = 0; i < surface_pressure_tuple_.size(); i++)
    {
        int body_index = std::get<0>(surface_pressure_tuple_[i]);
        SharedPtr<TriangleMeshShape> tri_mesh = std::get<1>(surface_pressure_tuple_[i]);
        Vec3d point = std::get<2>(surface_pressure_tuple_[i]);
        StdVec<std::array<Real, 2>> pressure_over_time = std::get<3>(surface_pressure_tuple_[i]);

        BodyPartByParticle *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), tri_mesh);
        surface_pressure_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::SurfacePressureFromSource>>(*bp, point, pressure_over_time));
    }
}

void StructuralSimulation::initializeSpringDamperConstraintParticleWise()
{
    spring_damper_constraint_ = {};
    for (size_t i = 0; i < spring_damper_tuple_.size(); i++)
    {
        SolidBody *solid_body = solid_body_list_[std::get<0>(spring_damper_tuple_[i])]->getSolidBodyFromMesh();
        spring_damper_constraint_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::SpringDamperConstraintParticleWise>>(
            *solid_body, std::get<1>(spring_damper_tuple_[i]), std::get<2>(spring_damper_tuple_[i])));
    }
}

void StructuralSimulation::initializeSpringNormalOnSurfaceParticles()
{
    surface_spring_ = {};
    for (size_t i = 0; i < surface_spring_tuple_.size(); i++)
    {
        int body_index = std::get<0>(surface_spring_tuple_[i]);
        SharedPtr<TriangleMeshShape> tri_mesh = std::get<1>(surface_spring_tuple_[i]);
        bool inner_outer = std::get<2>(surface_spring_tuple_[i]);
        Vec3d point = std::get<3>(surface_spring_tuple_[i]);
        Real spring_stiffness = std::get<4>(surface_spring_tuple_[i]);
        Real damping_coefficient = std::get<5>(surface_spring_tuple_[i]);

        surface_spring_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::SpringNormalOnSurfaceParticles>>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(),
            inner_outer, point, spring_stiffness, damping_coefficient));
    }
}

void StructuralSimulation::initializeConstrainSolidBody()
{
    fixed_constraint_body_ = {};
    for (size_t i = 0; i < body_indices_fixed_constraint_.size(); i++)
    {
        int body_index = body_indices_fixed_constraint_[i];
        fixed_constraint_body_.emplace_back(makeShared<SimpleDynamics<FixBodyConstraint>>(*solid_body_list_[body_index]->getSolidBodyFromMesh()));
    }
}

void StructuralSimulation::initializeConstrainSolidBodyRegion()
{
    fixed_constraint_region_ = {};
    for (size_t i = 0; i < body_indices_fixed_constraint_region_.size(); i++)
    {
        int body_index = body_indices_fixed_constraint_region_[i].first;
        BoundingBoxd bbox = body_indices_fixed_constraint_region_[i].second;

        // get the length of each side to create the box
        Real x_side = bbox.upper_[0] - bbox.lower_[0];
        Real y_side = bbox.upper_[1] - bbox.lower_[1];
        Real z_side = bbox.upper_[2] - bbox.lower_[2];
        Vec3d halfsize_bbox(0.5 * x_side, 0.5 * y_side, 0.5 * z_side);
        // get the center point for translation from the origin
        Vec3d center = (bbox.upper_ + bbox.lower_) * 0.5;
        // SimTK geometric modeling resolution
        int resolution(20);
        // create the triangle mesh of the box
        BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), makeShared<TriangleMeshShapeBrick>(halfsize_bbox, resolution, center, imported_stl_list_[body_index]));
        fixed_constraint_region_.emplace_back(makeShared<SimpleDynamics<FixBodyPartConstraint>>(*bp));
    }
}

void StructuralSimulation::initializePositionSolidBody()
{
    position_solid_body_ = {};
    for (size_t i = 0; i < position_solid_body_tuple_.size(); i++)
    {
        int body_index = std::get<0>(position_solid_body_tuple_[i]);
        Real start_time = std::get<1>(position_solid_body_tuple_[i]);
        Real end_time = std::get<2>(position_solid_body_tuple_[i]);
        Vecd pos_end_center = std::get<3>(position_solid_body_tuple_[i]);
        position_solid_body_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::PositionSolidBody>>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), start_time, end_time, pos_end_center));
    }
}

void StructuralSimulation::initializePositionScaleSolidBody()
{
    position_scale_solid_body_ = {};
    for (size_t i = 0; i < position_scale_solid_body_tuple_.size(); i++)
    {
        int body_index = std::get<0>(position_scale_solid_body_tuple_[i]);
        Real start_time = std::get<1>(position_scale_solid_body_tuple_[i]);
        Real end_time = std::get<2>(position_scale_solid_body_tuple_[i]);
        Real scale = std::get<3>(position_scale_solid_body_tuple_[i]);
        position_scale_solid_body_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::PositionScaleSolidBody>>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), start_time, end_time, scale));
    }
}

void StructuralSimulation::initializeTranslateSolidBody()
{
    translation_solid_body_ = {};
    for (size_t i = 0; i < translation_solid_body_tuple_.size(); i++)
    {
        int body_index = std::get<0>(translation_solid_body_tuple_[i]);
        Real start_time = std::get<1>(translation_solid_body_tuple_[i]);
        Real end_time = std::get<2>(translation_solid_body_tuple_[i]);
        Vecd translation = std::get<3>(translation_solid_body_tuple_[i]);
        translation_solid_body_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::TranslateSolidBody>>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), start_time, end_time, translation));
    }
}

void StructuralSimulation::initializeTranslateSolidBodyPart()
{
    translation_solid_body_part_ = {};
    for (size_t i = 0; i < translation_solid_body_part_tuple_.size(); i++)
    {
        int body_index = std::get<0>(translation_solid_body_part_tuple_[i]);
        Real start_time = std::get<1>(translation_solid_body_part_tuple_[i]);
        Real end_time = std::get<2>(translation_solid_body_part_tuple_[i]);
        Vecd translation = std::get<3>(translation_solid_body_part_tuple_[i]);
        BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(
            *solid_body_list_[body_index]->getSolidBodyFromMesh(), body_mesh_list_[body_index]);

        translation_solid_body_part_.emplace_back(makeShared<SimpleDynamics<solid_dynamics::TranslateSolidBodyPart>>(
            *bp, start_time, end_time, translation));
    }
}

void StructuralSimulation::executeInitialNormalDirection()
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getInitialNormalDirection()->exec();
    }
}

void StructuralSimulation::executeCorrectConfiguration()
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getCorrectConfiguration()->exec();
    }
}

void StructuralSimulation::executeUpdateElasticNormalDirection()
{
    for (size_t i = 0; i < particle_normal_update_.size(); i++)
    {
        particle_normal_update_[i]->exec();
    }
}

void StructuralSimulation::executeInitializeGravity()
{
    for (size_t i = 0; i < initialize_gravity_.size(); i++)
    {
        initialize_gravity_[i]->exec();
    }
}

void StructuralSimulation::executeExternalForceInBoundingBox()
{
    for (size_t i = 0; i < force_bounding_box_.size(); i++)
    {
        force_bounding_box_[i]->exec();
    }
}

void StructuralSimulation::executeForceInBodyRegion()
{
    for (size_t i = 0; i < force_in_body_region_.size(); i++)
    {
        force_in_body_region_[i]->exec();
    }
}

void StructuralSimulation::executeSurfacePressure()
{
    for (size_t i = 0; i < surface_pressure_.size(); i++)
    {
        surface_pressure_[i]->exec();
    }
}

void StructuralSimulation::executeSpringDamperConstraintParticleWise()
{
    for (size_t i = 0; i < spring_damper_constraint_.size(); i++)
    {
        spring_damper_constraint_[i]->exec();
    }
}

void StructuralSimulation::executeSpringNormalOnSurfaceParticles()
{
    for (size_t i = 0; i < surface_spring_.size(); i++)
    {
        surface_spring_[i]->exec();
    }
}

void StructuralSimulation::executeContactFactorSummation()
{
    // number of contacts that are not time dependent: contact pairs * 2
    size_t number_of_general_contacts = contacting_body_pairs_list_.size();
    for (size_t i = 0; i < contact_density_list_.size(); i++)
    {
        if (i < number_of_general_contacts)
        {
            contact_density_list_[i]->exec();
        }
        else
        {
            // index of the time dependent contact body pair
            // for i = 0, 1 --> index = 0, i = 2, 3 --> index = 1, and so on..
            int index = (i - number_of_general_contacts) / 2;
            Real start_time = time_dep_contacting_body_pairs_list_[index].second[0];
            Real end_time = time_dep_contacting_body_pairs_list_[index].second[1];
            if (physical_time_ >= start_time && physical_time_ <= end_time)
            {
                contact_density_list_[i]->exec();
            }
        }
    }
}

void StructuralSimulation::executeContactForce()
{
    // number of contacts that are not time dependent: contact pairs * 2
    size_t number_of_general_contacts = contacting_body_pairs_list_.size();
    for (size_t i = 0; i < contact_force_list_.size(); i++)
    {
        if (i < number_of_general_contacts)
        {
            contact_force_list_[i]->exec();
        }
        else
        {
            // index of the time dependent contact body pair
            // for i = 0, 1 --> index = 0, i = 2, 3 --> index = 1, and so on..
            int index = (i - number_of_general_contacts) / 2;
            Real start_time = time_dep_contacting_body_pairs_list_[index].second[0];
            Real end_time = time_dep_contacting_body_pairs_list_[index].second[1];
            if (physical_time_ >= start_time && physical_time_ <= end_time)
            {
                contact_force_list_[i]->exec();
            }
        }
    }
}

void StructuralSimulation::executeStressRelaxationFirstHalf(Real dt)
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getStressRelaxationFirstHalf()->exec(dt);
    }
}

void StructuralSimulation::executeConstrainSolidBody()
{
    for (size_t i = 0; i < fixed_constraint_body_.size(); i++)
    {
        fixed_constraint_body_[i]->exec();
    }
}

void StructuralSimulation::executeConstrainSolidBodyRegion()
{
    for (size_t i = 0; i < fixed_constraint_region_.size(); i++)
    {
        fixed_constraint_region_[i]->exec();
    }
}

void StructuralSimulation::executePositionSolidBody(Real dt)
{
    for (size_t i = 0; i < position_solid_body_.size(); i++)
    {
        position_solid_body_[i]->exec(dt);
    }
}

void StructuralSimulation::executePositionScaleSolidBody(Real dt)
{
    for (size_t i = 0; i < position_scale_solid_body_.size(); i++)
    {
        position_scale_solid_body_[i]->exec(dt);
    }
}

void StructuralSimulation::executeTranslateSolidBody(Real dt)
{
    for (size_t i = 0; i < translation_solid_body_.size(); i++)
    {
        translation_solid_body_[i]->exec(dt);
    }
}

void StructuralSimulation::executeTranslateSolidBodyPart(Real dt)
{
    for (size_t i = 0; i < translation_solid_body_part_.size(); i++)
    {
        translation_solid_body_part_[i]->exec(dt);
    }
}

void StructuralSimulation::executeDamping(Real dt)
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getDampingWithRandomChoice()->exec(dt);
    }
}

void StructuralSimulation::executeStressRelaxationSecondHalf(Real dt)
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getStressRelaxationSecondHalf()->exec(dt);
    }
}

void StructuralSimulation::executeUpdateCellLinkedList()
{
    for (size_t i = 0; i < solid_body_list_.size(); i++)
    {
        solid_body_list_[i]->getSolidBodyFromMesh()->updateCellLinkedList();
    }
}

void StructuralSimulation::executeContactUpdateConfiguration()
{
    // number of contacts that are not time dependent: contact pairs * 2
    size_t number_of_general_contacts = contacting_body_pairs_list_.size();
    for (size_t i = 0; i < contact_list_.size(); i++)
    {
        // general contacts = contacting_bodies * 2
        if (i < number_of_general_contacts)
        {
            contact_list_[i]->updateConfiguration();
        }
        // time dependent contacts = time dep. contacting_bodies * 2
        else
        {
            // index of the time dependent contact body pair
            // for i = 0, 1 --> index = 0, i = 2, 3 --> index = 1, and so on..
            int index = (i - number_of_general_contacts) / 2;
            Real start_time = time_dep_contacting_body_pairs_list_[index].second[0];
            Real end_time = time_dep_contacting_body_pairs_list_[index].second[1];
            if (physical_time_ >= start_time && physical_time_ <= end_time)
            {
                contact_list_[i]->updateConfiguration();
            }
        }
    }
}

void StructuralSimulation::initializeSimulation()
{
    physical_time_ = 0.0;

    /** INITIALIZE SYSTEM */
    system_.initializeSystemCellLinkedLists();
    system_.initializeSystemConfigurations();

    /** INITIAL CONDITION */
    executeCorrectConfiguration();
}

void StructuralSimulation::runSimulationStep(Real &dt, Real &integration_time)
{
    if (iteration_ % 100 == 0)
        std::cout << "N=" << iteration_ << " Time: " << physical_time_ << "	dt: " << dt << "\n";

    /** UPDATE NORMAL DIRECTIONS */
    executeUpdateElasticNormalDirection();

    /** ACTIVE BOUNDARY CONDITIONS */
    // force (acceleration) based
    executeInitializeGravity();
    executeExternalForceInBoundingBox();
    executeForceInBodyRegion();
    executeSurfacePressure();
    executeSpringDamperConstraintParticleWise();
    executeSpringNormalOnSurfaceParticles();

    /** CONTACT */
    executeContactFactorSummation();
    executeContactForce();

    /** STRESS RELAXATION, DAMPING, POSITIONAL CONSTRAINTS */
    executeStressRelaxationFirstHalf(dt);

    executeConstrainSolidBody();
    executeConstrainSolidBodyRegion();
    executePositionSolidBody(dt);
    executePositionScaleSolidBody(dt);
    executeTranslateSolidBody(dt);
    // velocity based
    executeTranslateSolidBodyPart(dt);

    executeDamping(dt);

    executeConstrainSolidBody();
    executeConstrainSolidBodyRegion();
    executePositionSolidBody(dt);
    executePositionScaleSolidBody(dt);
    executeTranslateSolidBody(dt);
    // velocity based
    executeTranslateSolidBodyPart(dt);

    executeStressRelaxationSecondHalf(dt);
    /** UPDATE TIME STEP SIZE, INCREMENT */
    iteration_++;
    dt = system_.getSmallestTimeStepAmongSolidBodies();
    integration_time += dt;
    physical_time_ += dt;

    /** UPDATE BODIES CELL LINKED LISTS */
    executeUpdateCellLinkedList();

    /** UPDATE CONTACT CONFIGURATION */
    executeContactUpdateConfiguration();
}

void StructuralSimulation::runSimulation(Real end_time)
{
    BodyStatesRecordingToVtp write_states(system_);

    /** Statistics for computing time. */
    write_states.writeToFile(0);
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /** Main loop */
    while (physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            runSimulationStep(dt, integration_time);
        }
        TickCount t2 = TickCount::now();

        // write data to file
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
}

double StructuralSimulation::runSimulationFixedDurationJS(int number_of_steps)
{
    BodyStatesRecordingToVtp write_states(system_);
    physical_time_ = 0.0;

    /** Statistics for computing time. */
    write_states.writeToFile(0);
    int output_period = 100;
    Real dt = 0.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /** Main loop */
    while (iteration_ < number_of_steps)
    {
        Real integration_time = 0.0;
        int output_step = 0;
        while (output_step < output_period)
        {
            runSimulationStep(dt, integration_time);
            output_step++;
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    return tt.seconds();
}

Real StructuralSimulation::getMaxDisplacement(int body_index)
{
    BaseParticles *base_particles = solid_body_list_[body_index].get()->getElasticSolidParticles();
    Vecd *pos = base_particles->ParticlePositions();
    Vecd *pos0 = base_particles->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position");
    Real displ_max = 0;
    for (size_t i = 0; i < base_particles->TotalRealParticles(); i++)
    {
        Real displ = (pos[i] - pos0[i]).norm();
        if (displ > displ_max)
            displ_max = displ;
    }
    return displ_max;
}

StructuralSimulationJS::StructuralSimulationJS(const StructuralSimulationInput &input)
    : StructuralSimulation(input),
      write_states_(system_),
      dt(0.0)
{
    write_states_.writeToFile(0);
    physical_time_ = 0.0;
}

void StructuralSimulationJS::runSimulationFixedDuration(int number_of_steps)
{
    /** Statistics for computing time. */
    TickCount t_start = TickCount::now(); // computation time monitoring
    /** Main loop */
    for (int i = 0; i < number_of_steps; ++i)
    {
        Real integration_time = 0.0;
        runSimulationStep(dt, integration_time);
    }
    TickCount t_end = TickCount::now(); // computation time monitoring
    TimeInterval t_interval = t_end - t_start;
    std::cout << "Total time for computation: " << t_interval.seconds() << " seconds." << std::endl;
}

VtuStringData StructuralSimulationJS::getVtuData()
{
    write_states_.writeToFile();
    const auto vtuData = write_states_.GetVtuData();
    write_states_.clear();
    return vtuData;
}
