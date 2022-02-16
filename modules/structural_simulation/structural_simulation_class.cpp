/**
* @file 	structural_simulation_class.cpp
* @brief 	The structural simulation module is licensed under the Aladdin Free Public License (https://spdx.org/licenses/Aladdin.html) regarding usage for medical device development.
* Commercial use for medical device development is not permitted. This does not apply to applications in other fields.
* @details	solid structural simulation class for general structural simulations
* @author 	Bence Z. Rochlitz - Virtonomy GmbH
*/

#include "structural_simulation_class.h"

////////////////////////////////////////////////////
/* global functions in StructuralSimulation  */
////////////////////////////////////////////////////

BodyPartFromMesh::BodyPartFromMesh(SPHBody &body, const string &body_part_name, TriangleMeshShape &triangle_mesh_shape)
	: BodyRegionByParticle(body, body_part_name, triangle_mesh_shape)
{
	// set the body domain bounds because it is not set by default
	BoundingBox bounds = body_part_shape_.findBounds();
	setBodyPartBounds(bounds);
}

SolidBodyFromMesh::SolidBodyFromMesh(
	SPHSystem &system, const string &body_name, TriangleMeshShape &triangle_mesh_shape,
	SharedPtr<SPHAdaptation> particle_adaptation, StdLargeVec<Vecd> &pos_0, StdLargeVec<Real> &volume)
	: SolidBody(system, body_name, particle_adaptation)
{
	body_shape_.add<LevelSetShape>(this, triangle_mesh_shape, true, false);

	// set the body domain bounds because it is not set by default
	BoundingBox bounds = body_shape_.findBounds();
	setBodyDomainBounds(bounds);
}

SolidBodyForSimulation::SolidBodyForSimulation(
	SPHSystem &system, const string &body_name, TriangleMeshShape &triangle_mesh_shape, SharedPtr<SPHAdaptation> particle_adaptation,
	Real physical_viscosity, SharedPtr<LinearElasticSolid> material_model, StdLargeVec<Vecd> &pos_0, StdLargeVec<Real> &volume)
	: solid_body_from_mesh_(system, body_name, triangle_mesh_shape, particle_adaptation, pos_0, volume),
	  //material_model_(material_model),
	  elastic_solid_particles_(solid_body_from_mesh_, material_model),
	  inner_body_relation_(solid_body_from_mesh_),

	  correct_configuration_(solid_dynamics::CorrectConfiguration(inner_body_relation_)),
	  stress_relaxation_first_half_(solid_dynamics::StressRelaxationFirstHalf(inner_body_relation_)),
	  stress_relaxation_second_half_(solid_dynamics::StressRelaxationSecondHalf(inner_body_relation_)),
	  damping_random_(DampingWithRandomChoice<DampingPairwiseInner<Vec3d>>(inner_body_relation_, 0.2, "Velocity", physical_viscosity))
{
}

void expandBoundingBox(BoundingBox *original, BoundingBox *additional)
{
	for (int i = 0; i < original->first.size(); i++)
	{
		if (additional->first[i] < original->first[i])
		{
			original->first[i] = additional->first[i];
		}
		if (additional->second[i] > original->second[i])
		{
			original->second[i] = additional->second[i];
		}
	}
}

void relaxParticlesSingleResolution(In_Output &in_output,
									bool write_particle_relaxation_data,
									SolidBody &solid_body_from_mesh,
									BodyRelationInner &solid_body_from_mesh_inner)
{
	BodyStatesRecordingToVtp write_solid_body_from_mesh_to_vtp(in_output, solid_body_from_mesh);

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_solid_body_from_mesh_particles(solid_body_from_mesh);
	/** A  Physics relaxation step. */
	relax_dynamics::SolidRelaxationStepInner relaxation_step_inner(solid_body_from_mesh_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_solid_body_from_mesh_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
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
		relaxation_step_inner.parallel_exec();
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

std::tuple<StdLargeVec<Vecd>, StdLargeVec<Real>> generateAndRelaxParticlesFromMesh(
	TriangleMeshShape &triangle_mesh_shape, Real resolution, bool particle_relaxation, bool write_particle_relaxation_data)
{
	BoundingBox bb = triangle_mesh_shape.findBounds();
	SPHSystem system(bb, resolution);

	SharedPtr<SPHAdaptation> particle_adaptation = makeShared<SPHAdaptation>(1.15, 1.0);

	std::string name = "model";
	SolidBody model(system, name, particle_adaptation);
	model.body_shape_.add<LevelSetShape>(&model, triangle_mesh_shape, true, false);

	// set the body domain bounds because it is not set by default
	BoundingBox bounds = model.body_shape_.findBounds();
	model.setBodyDomainBounds(bounds);

	SolidParticles particles(model);

	if (particle_relaxation)
	{
		In_Output in_output(system);
		BodyRelationInner inner_relation(model);
		relaxParticlesSingleResolution(in_output, write_particle_relaxation_data, model, inner_relation);
	}

	return std::tuple<StdLargeVec<Vecd>, StdLargeVec<Real>>(particles.pos_0_, particles.Vol_);
}

StructuralSimulationInput::StructuralSimulationInput(
	const string &relative_input_path,
	const vector<string> &imported_stl_list,
	Real scale_stl,
	const vector<Vec3d> &translation_list,
	const vector<Real> &resolution_list,
	const vector<SharedPtr<LinearElasticSolid>> &material_model_list,
	const StdVec<Real> &physical_viscosity,
	const StdVec<IndexVector> &contacting_bodies_list)
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
	acceleration_bounding_box_tuple_ = {};
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

StructuralSimulation::StructuralSimulation(StructuralSimulationInput &input)
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
	  system_(SPHSystem(BoundingBox(Vec3d(0), Vec3d(0)), system_resolution_)),
	  scale_system_boundaries_(input.scale_system_boundaries_),
	  in_output_(In_Output(system_)),

	  // optional: boundary conditions
	  non_zero_gravity_(input.non_zero_gravity_),
	  acceleration_bounding_box_tuple_(input.acceleration_bounding_box_tuple_),
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
	  iteration_(0),

	  // data storage
	  von_mises_stress_max_({}),
	  von_mises_stress_particles_({})

{
	// scaling of translation and resolution
	scaleTranslationAndResolution();
	// set the default resolution to the max in the resolution list
	setSystemResolutionMax();
	// create the body mesh list for triangular mesh shapes storage
	createBodyMeshList();
	// create the particle adaptions for the bodies
	createParticleAdaptationList();
	// set up the system
	calculateSystemBoundaries();
	system_.run_particle_relaxation_ = true;
	// initialize solid bodies with their properties
	initializeElasticSolidBodies();
	// contacts
	initializeAllContacts();

	// boundary conditions
	initializeGravity();
	initializeAccelerationForBodyPartInBoundingBox();
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
	system_.resolution_ref_ = system_resolution_;
}

void StructuralSimulation::calculateSystemBoundaries()
{
	// calculate system bounds from all bodies
	for (size_t i = 0; i < body_mesh_list_.size(); i++)
	{
		BoundingBox additional = body_mesh_list_[i]->findBounds();
		expandBoundingBox(&system_.system_domain_bounds_, &additional);
	}
	// scale the system bounds around the center point
	Vecd center_point = (system_.system_domain_bounds_.first + system_.system_domain_bounds_.second) * 0.5;

	Vecd distance_first = system_.system_domain_bounds_.first - center_point;
	Vecd distance_second = system_.system_domain_bounds_.second - center_point;

	system_.system_domain_bounds_.first = center_point + distance_first * scale_system_boundaries_;
	system_.system_domain_bounds_.second = center_point + distance_second * scale_system_boundaries_;
}

void StructuralSimulation::createBodyMeshList()
{
	body_mesh_list_ = {};
	for (size_t i = 0; i < imported_stl_list_.size(); i++)
	{
		string relative_input_path_copy = relative_input_path_;
		TriangleMeshShape *tri_mesh_shape = 
			tri_mesh_shape_ptr_keeper_.createPtr<TriangleMeshShapeSTL>(relative_input_path_copy.append(imported_stl_list_[i]), translation_list_[i], scale_stl_);
		body_mesh_list_.push_back(tri_mesh_shape);
	}
}

void StructuralSimulation::createParticleAdaptationList()
{
	particle_adaptation_list_ = {};
	for (size_t i = 0; i < resolution_list_.size(); i++)
	{
		Real system_resolution_ratio = system_resolution_ / resolution_list_[i];
		// for solid bodies, slightly small h_spaing_ratio is used
		particle_adaptation_list_.push_back(makeShared<SPHAdaptation>(1.15, system_resolution_ratio));
	}
}

void StructuralSimulation::initializeElasticSolidBodies()
{
	solid_body_list_ = {};
	particle_normal_update_ = {};
	for (size_t i = 0; i < body_mesh_list_.size(); i++)
	{
		// create the initial particles from the triangle mesh shape with particle relaxation option
		std::tuple<StdLargeVec<Vecd>, StdLargeVec<Real>> particles = generateAndRelaxParticlesFromMesh(*body_mesh_list_[i], resolution_list_[i], particle_relaxation_list_[i], write_particle_relaxation_data_);

		// get the particles' initial position and their volume
		StdLargeVec<Vecd> &pos_0 = std::get<0>(particles);
		StdLargeVec<Real> &volume = std::get<1>(particles);

		// create the SolidBodyForSimulation
		solid_body_list_.emplace_back(make_shared<SolidBodyForSimulation>(
			system_, imported_stl_list_[i], *body_mesh_list_[i], particle_adaptation_list_[i], physical_viscosity_[i], material_model_list_[i], pos_0, volume));

		// update initial normal direction of particles
		solid_body_list_[i]->getElasticSolidParticles()->initializeNormalDirectionFromBodyShape();
		std::cout << "  normal initialization done" << std::endl;
		particle_normal_update_.emplace_back(make_shared<solid_dynamics::UpdateElasticNormalDirection>(*solid_body_list_[i]->getSolidBodyFromMesh()));
	}
}

void StructuralSimulation::initializeContactBetweenTwoBodies(int first, int second)
{
	SolidBodyFromMesh *first_body = solid_body_list_[first]->getSolidBodyFromMesh();
	SolidBodyFromMesh *second_body = solid_body_list_[second]->getSolidBodyFromMesh();
	SolidBodyRelationContact *first_contact = contact_relation_ptr_keeper_.createPtr<SolidBodyRelationContact>(*first_body, RealBodyVector({second_body}));
	SolidBodyRelationContact *second_contact = contact_relation_ptr_keeper_.createPtr<SolidBodyRelationContact>(*second_body, RealBodyVector({first_body}));

	contact_list_.emplace_back(first_contact);
	contact_list_.emplace_back(second_contact);

	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(*first_contact));
	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(*second_contact));

	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(*first_contact));
	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(*second_contact));
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

		SolidBodyRelationContact *contact = contact_relation_ptr_keeper_.createPtr<SolidBodyRelationContact>(*contact_body, target_list);

		contact_list_.emplace_back(contact);
		contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(*contact));
		contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(*contact));
	}
	// continue appending the lists with the time dependent contacts
	for (size_t i = 0; i < time_dep_contacting_body_pairs_list_.size(); i++)
	{
		int body_1 = time_dep_contacting_body_pairs_list_[i].first[0];
		int body_2 = time_dep_contacting_body_pairs_list_[i].first[1];
		initializeContactBetweenTwoBodies(body_1, body_2); //vector with first element being array with indices
	}
}

void StructuralSimulation::initializeGravity()
{
	// collect all the body indices with non-zero gravity
	vector<int> gravity_indices = {};
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
			Gravity *gravity = gravity_ptr_keeper_.createPtr<Gravity>(non_zero_gravity_[gravity_index_i].second);
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(*solid_body_list_[i]->getSolidBodyFromMesh(), *gravity));
			gravity_index_i++;
		}
		else
		{
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(*solid_body_list_[i]->getSolidBodyFromMesh()));
		}
	}
}

void StructuralSimulation::initializeAccelerationForBodyPartInBoundingBox()
{
	acceleration_bounding_box_ = {};
	for (size_t i = 0; i < acceleration_bounding_box_tuple_.size(); i++)
	{
		SolidBody *solid_body = solid_body_list_[get<0>(acceleration_bounding_box_tuple_[i])]->getSolidBodyFromMesh();
		acceleration_bounding_box_.emplace_back(make_shared<solid_dynamics::AccelerationForBodyPartInBoundingBox>(
			*solid_body, get<1>(acceleration_bounding_box_tuple_[i]), get<2>(acceleration_bounding_box_tuple_[i])));
	}
}

void StructuralSimulation::initializeForceInBodyRegion()
{
	force_in_body_region_ = {};
	for (size_t i = 0; i < force_in_body_region_tuple_.size(); i++)
	{
		int body_index = get<0>(force_in_body_region_tuple_[i]);
		BoundingBox bbox = get<1>(force_in_body_region_tuple_[i]);
		Vec3d force = get<2>(force_in_body_region_tuple_[i]);
		Real end_time = get<3>(force_in_body_region_tuple_[i]);

		// get the length of each side to create the box
		Real x_side = bbox.second[0] - bbox.first[0];
		Real y_side = bbox.second[1] - bbox.first[1];
		Real z_side = bbox.second[2] - bbox.first[2];
		Vec3d halfsize_bbox(0.5 * x_side, 0.5 * y_side, 0.5 * z_side);
		// get the center point for translation from the origin
		Vec3d center = (bbox.second + bbox.first) * 0.5;
		// SimTK geometric modeling resolution
		int resolution(20);
		// create the triangle mesh of the box
		TriangleMeshShape *tri_mesh = tri_mesh_shape_ptr_keeper_.createPtr<TriangleMeshShapeBrick>(halfsize_bbox, resolution, center);
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *tri_mesh);
		force_in_body_region_.emplace_back(make_shared<solid_dynamics::ForceInBodyRegion>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, force, end_time));
	}
}

void StructuralSimulation::initializeSurfacePressure()
{
	surface_pressure_ = {};
	for (size_t i = 0; i < surface_pressure_tuple_.size(); i++)
	{
		int body_index = get<0>(surface_pressure_tuple_[i]);
		TriangleMeshShape *tri_mesh = get<1>(surface_pressure_tuple_[i]);
		Vec3d point = get<2>(surface_pressure_tuple_[i]);
		StdVec<array<Real, 2>> pressure_over_time = get<3>(surface_pressure_tuple_[i]);

		BodyPartByParticle *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *tri_mesh);
		surface_pressure_.emplace_back(make_shared<solid_dynamics::SurfacePressureFromSource>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, point, pressure_over_time));
	}
}

void StructuralSimulation::initializeSpringDamperConstraintParticleWise()
{
	spring_damper_constraint_ = {};
	for (size_t i = 0; i < spring_damper_tuple_.size(); i++)
	{
		SolidBody *solid_body = solid_body_list_[get<0>(spring_damper_tuple_[i])]->getSolidBodyFromMesh();
		spring_damper_constraint_.emplace_back(make_shared<solid_dynamics::SpringDamperConstraintParticleWise>(*solid_body, get<1>(spring_damper_tuple_[i]), get<2>(spring_damper_tuple_[i])));
	}
}

void StructuralSimulation::initializeSpringNormalOnSurfaceParticles()
{
	surface_spring_ = {};
	for (size_t i = 0; i < surface_spring_tuple_.size(); i++)
	{
		int body_index = get<0>(surface_spring_tuple_[i]);
		TriangleMeshShape *tri_mesh = get<1>(surface_spring_tuple_[i]);
		bool inner_outer = get<2>(surface_spring_tuple_[i]);
		Vec3d point = get<3>(surface_spring_tuple_[i]);
		Real spring_stiffness = get<4>(surface_spring_tuple_[i]);
		Real damping_coefficient = get<5>(surface_spring_tuple_[i]);

		BodyPartByParticle *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *tri_mesh);
		surface_spring_.emplace_back(make_shared<solid_dynamics::SpringNormalOnSurfaceParticles>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, inner_outer, point, spring_stiffness, damping_coefficient));
	}
}

void StructuralSimulation::initializeConstrainSolidBody()
{
	fixed_constraint_body_ = {};
	for (size_t i = 0; i < body_indices_fixed_constraint_.size(); i++)
	{
		int body_index = body_indices_fixed_constraint_[i];
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *body_mesh_list_[body_index]);
		fixed_constraint_body_.emplace_back(make_shared<solid_dynamics::ConstrainSolidBodyRegion>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp));
	}
}

void StructuralSimulation::initializeConstrainSolidBodyRegion()
{
	fixed_constraint_region_ = {};
	for (size_t i = 0; i < body_indices_fixed_constraint_region_.size(); i++)
	{
		int body_index = body_indices_fixed_constraint_region_[i].first;
		BoundingBox bbox = body_indices_fixed_constraint_region_[i].second;

		// get the length of each side to create the box
		Real x_side = bbox.second[0] - bbox.first[0];
		Real y_side = bbox.second[1] - bbox.first[1];
		Real z_side = bbox.second[2] - bbox.first[2];
		Vec3d halfsize_bbox(0.5 * x_side, 0.5 * y_side, 0.5 * z_side);
		// get the center point for translation from the origin
		Vec3d center = (bbox.second + bbox.first) * 0.5;
		// SimTK geometric modeling resolution
		int resolution(20);
		// create the triangle mesh of the box
		TriangleMeshShape *tri_mesh = tri_mesh_shape_ptr_keeper_.createPtr<TriangleMeshShapeBrick>(halfsize_bbox, resolution, center);

		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *tri_mesh);
		fixed_constraint_region_.emplace_back(make_shared<solid_dynamics::ConstrainSolidBodyRegion>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp));
	}
}

void StructuralSimulation::initializePositionSolidBody()
{
	position_solid_body_ = {};
	for (size_t i = 0; i < position_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(position_solid_body_tuple_[i]);
		Real start_time = get<1>(position_solid_body_tuple_[i]);
		Real end_time = get<2>(position_solid_body_tuple_[i]);
		Vecd pos_end_center = get<3>(position_solid_body_tuple_[i]);
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *body_mesh_list_[body_index]);

		position_solid_body_.emplace_back(make_shared<solid_dynamics::PositionSolidBody>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, start_time, end_time, pos_end_center));
	}
}

void StructuralSimulation::initializePositionScaleSolidBody()
{
	position_scale_solid_body_ = {};
	for (size_t i = 0; i < position_scale_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(position_scale_solid_body_tuple_[i]);
		Real start_time = get<1>(position_scale_solid_body_tuple_[i]);
		Real end_time = get<2>(position_scale_solid_body_tuple_[i]);
		Real scale = get<3>(position_scale_solid_body_tuple_[i]);
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *body_mesh_list_[body_index]);

		position_scale_solid_body_.emplace_back(make_shared<solid_dynamics::PositionScaleSolidBody>(*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, start_time, end_time, scale));
	}
}

void StructuralSimulation::initializeTranslateSolidBody()
{
	translation_solid_body_ = {};
	for (size_t i = 0; i < translation_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(translation_solid_body_tuple_[i]);
		Real start_time = get<1>(translation_solid_body_tuple_[i]);
		Real end_time = get<2>(translation_solid_body_tuple_[i]);
		Vecd translation = get<3>(translation_solid_body_tuple_[i]);
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(
			*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *body_mesh_list_[body_index]);

		translation_solid_body_.emplace_back(make_shared<solid_dynamics::TranslateSolidBody>(
			*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, start_time, end_time, translation));
	}
}

void StructuralSimulation::initializeTranslateSolidBodyPart()
{
	translation_solid_body_part_ = {};
	for (size_t i = 0; i < translation_solid_body_part_tuple_.size(); i++)
	{
		int body_index = get<0>(translation_solid_body_part_tuple_[i]);
		Real start_time = get<1>(translation_solid_body_part_tuple_[i]);
		Real end_time = get<2>(translation_solid_body_part_tuple_[i]);
		Vecd translation = get<3>(translation_solid_body_part_tuple_[i]);
		BoundingBox bbox = get<4>(translation_solid_body_part_tuple_[i]);
		BodyPartFromMesh *bp = body_part_tri_mesh_ptr_keeper_.createPtr<BodyPartFromMesh>(
			*solid_body_list_[body_index]->getSolidBodyFromMesh(), imported_stl_list_[body_index], *body_mesh_list_[body_index]);

		translation_solid_body_part_.emplace_back(make_shared<solid_dynamics::TranslateSolidBodyPart>(
			*solid_body_list_[body_index]->getSolidBodyFromMesh(), *bp, start_time, end_time, translation, bbox));
	}
}

void StructuralSimulation::executeCorrectConfiguration()
{
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->getCorrectConfiguration()->parallel_exec();
	}
}

void StructuralSimulation::executeUpdateElasticNormalDirection()
{
	for (size_t i = 0; i < particle_normal_update_.size(); i++)
	{
		particle_normal_update_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeinitializeATimeStep()
{
	for (size_t i = 0; i < initialize_gravity_.size(); i++)
	{
		initialize_gravity_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeAccelerationForBodyPartInBoundingBox()
{
	for (size_t i = 0; i < acceleration_bounding_box_.size(); i++)
	{
		acceleration_bounding_box_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeForceInBodyRegion()
{
	for (size_t i = 0; i < force_in_body_region_.size(); i++)
	{
		force_in_body_region_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeSurfacePressure()
{
	for (size_t i = 0; i < surface_pressure_.size(); i++)
	{
		surface_pressure_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeSpringDamperConstraintParticleWise()
{
	for (size_t i = 0; i < spring_damper_constraint_.size(); i++)
	{
		spring_damper_constraint_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeSpringNormalOnSurfaceParticles()
{
	for (size_t i = 0; i < surface_spring_.size(); i++)
	{
		surface_spring_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeContactDensitySummation()
{
	// number of contacts that are not time dependent: contact pairs * 2
	size_t number_of_general_contacts = contacting_body_pairs_list_.size();
	for (size_t i = 0; i < contact_density_list_.size(); i++)
	{
		if (i < number_of_general_contacts)
		{
			contact_density_list_[i]->parallel_exec();
		}
		else
		{
			// index of the time dependent contact body pair
			// for i = 0, 1 --> index = 0, i = 2, 3 --> index = 1, and so on..
			int index = (i - number_of_general_contacts) / 2;
			Real start_time = time_dep_contacting_body_pairs_list_[index].second[0];
			Real end_time = time_dep_contacting_body_pairs_list_[index].second[1];
			if (GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
			{
				contact_density_list_[i]->parallel_exec();
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
			contact_force_list_[i]->parallel_exec();
		}
		else
		{
			// index of the time dependent contact body pair
			// for i = 0, 1 --> index = 0, i = 2, 3 --> index = 1, and so on..
			int index = (i - number_of_general_contacts) / 2;
			Real start_time = time_dep_contacting_body_pairs_list_[index].second[0];
			Real end_time = time_dep_contacting_body_pairs_list_[index].second[1];
			if (GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
			{
				contact_force_list_[i]->parallel_exec();
			}
		}
	}
}

void StructuralSimulation::executeStressRelaxationFirstHalf(Real dt)
{
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->getStressRelaxationFirstHalf()->parallel_exec(dt);
	}
}

void StructuralSimulation::executeConstrainSolidBody()
{
	for (size_t i = 0; i < fixed_constraint_body_.size(); i++)
	{
		fixed_constraint_body_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeConstrainSolidBodyRegion()
{
	for (size_t i = 0; i < fixed_constraint_region_.size(); i++)
	{
		fixed_constraint_region_[i]->parallel_exec();
	}
}

void StructuralSimulation::executePositionSolidBody(Real dt)
{
	for (size_t i = 0; i < position_solid_body_.size(); i++)
	{
		position_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::executePositionScaleSolidBody(Real dt)
{
	for (size_t i = 0; i < position_scale_solid_body_.size(); i++)
	{
		position_scale_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::executeTranslateSolidBody(Real dt)
{
	for (size_t i = 0; i < translation_solid_body_.size(); i++)
	{
		translation_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::executeTranslateSolidBodyPart(Real dt)
{
	for (size_t i = 0; i < translation_solid_body_part_.size(); i++)
	{
		translation_solid_body_part_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::executeDamping(Real dt)
{
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->getDampingWithRandomChoice()->parallel_exec(dt);
	}
}

void StructuralSimulation::executeStressRelaxationSecondHalf(Real dt)
{
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->getStressRelaxationSecondHalf()->parallel_exec(dt);
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
			if (GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
			{
				contact_list_[i]->updateConfiguration();
			}
		}
	}
}

void StructuralSimulation::initializeSimulation()
{
	GlobalStaticVariables::physical_time_ = 0.0;

	/** INITIALALIZE SYSTEM */
	system_.initializeSystemCellLinkedLists();
	system_.initializeSystemConfigurations();

	/** INITIAL CONDITION */
	executeCorrectConfiguration();
}

void StructuralSimulation::runSimulationStep(Real &dt, Real &integration_time)
{
	if (iteration_ % 100 == 0)
		cout << "N=" << iteration_ << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";

	/** UPDATE NORMAL DIRECTIONS */
	executeUpdateElasticNormalDirection();

	/** ACTIVE BOUNDARY CONDITIONS */
	// force (acceleration) based
	executeinitializeATimeStep();
	executeAccelerationForBodyPartInBoundingBox();
	executeForceInBodyRegion();
	executeSurfacePressure();
	executeSpringDamperConstraintParticleWise();
	executeSpringNormalOnSurfaceParticles();

	/** CONTACT */
	executeContactDensitySummation();
	executeContactForce();

	/** STRESS RELAXATOIN, DAMPING, POSITIONAL CONSTRAINTS */
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
	GlobalStaticVariables::physical_time_ += dt;

	/** UPDATE BODIES CELL LINKED LISTS */
	executeUpdateCellLinkedList();

	/** UPDATE CONTACT CONFIGURATION */
	executeContactUpdateConfiguration();
}

void StructuralSimulation::runSimulation(Real end_time)
{
	BodyStatesRecordingToVtp write_states(in_output_, system_.real_bodies_);

	/** Statistics for computing time. */
	write_states.writeToFile(0);
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			runSimulationStep(dt, integration_time);
		}
		tick_count t2 = tick_count::now();
		// record data for test
		von_mises_stress_max_.push_back(solid_body_list_[0].get()->getElasticSolidParticles()->getVonMisesStressMax());
		von_mises_stress_particles_.push_back(solid_body_list_[0].get()->getElasticSolidParticles()->getVonMisesStressVector());

		von_mises_strain_max_.push_back(solid_body_list_[0].get()->getElasticSolidParticles()->getVonMisesStrainMax());
		von_mises_strain_particles_.push_back(solid_body_list_[0].get()->getElasticSolidParticles()->getVonMisesStrainVector());
		// write data to file
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
}

double StructuralSimulation::runSimulationFixedDurationJS(int number_of_steps)
{
	BodyStatesRecordingToVtp write_states(in_output_, system_.real_bodies_);
	GlobalStaticVariables::physical_time_ = 0.0;

	/** Statistics for computing time. */
	write_states.writeToFile(0);
	int output_period = 100;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
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
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	return tt.seconds();
}

Real StructuralSimulation::getMaxDisplacement(int body_index)
{
	StdLargeVec<Vecd> &pos_0 = solid_body_list_[body_index].get()->getElasticSolidParticles()->pos_0_;
	StdLargeVec<Vecd> &pos_n = solid_body_list_[body_index].get()->getElasticSolidParticles()->pos_n_;
	Real displ_max = 0;
	for (size_t i = 0; i < pos_0.size(); i++)
	{
		Real displ = (pos_n[i] - pos_0[i]).norm();
		if (displ > displ_max)
			displ_max = displ;
	}
	return displ_max;
}
