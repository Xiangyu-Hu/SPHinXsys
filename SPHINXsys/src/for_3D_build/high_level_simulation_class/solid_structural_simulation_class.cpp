#include "solid_structural_simulation_class.h"

////////////////////////////////////////////////////
/* global functions in StructuralSimulation  */
////////////////////////////////////////////////////

BodyPartByParticleTriMesh::BodyPartByParticleTriMesh(SPHBody* body, string body_part_name, TriangleMeshShape* triangle_mesh_shape)
: BodyPartByParticle(body, body_part_name)
{	
	body_part_shape_ = new ComplexShape(body_part_name);
	body_part_shape_->addTriangleMeshShape(triangle_mesh_shape, ShapeBooleanOps::add);
	tagBodyPart();
}

BodyPartByParticleTriMesh::~BodyPartByParticleTriMesh()
{
	delete body_part_shape_;
}

ImportedModel::ImportedModel(SPHSystem &system, string body_name, TriangleMeshShape* triangle_mesh_shape, ParticleAdaptation* particle_adaptation)
	: SolidBody(system, body_name, particle_adaptation)
{
	ComplexShape original_body_shape;
	original_body_shape.addTriangleMeshShape(triangle_mesh_shape, ShapeBooleanOps::add);
	body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
}

ImportedModel::~ImportedModel()
{
	delete body_shape_;
}

SolidBodyForSimulation::SolidBodyForSimulation(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, ParticleAdaptation& particle_adaptation, Real physical_viscosity, LinearElasticSolid& material_model):
	imported_model_(ImportedModel(system, body_name, &triangle_mesh_shape, &particle_adaptation)),
	//material_model_(material_model),
	elastic_solid_particles_(ElasticSolidParticles(&imported_model_, &material_model)),
	inner_body_relation_(BodyRelationInner(&imported_model_)),

	correct_configuration_(solid_dynamics::CorrectConfiguration(&inner_body_relation_)),
	stress_relaxation_first_half_(solid_dynamics::StressRelaxationFirstHalf(&inner_body_relation_)),
	stress_relaxation_second_half_(solid_dynamics::StressRelaxationSecondHalf(&inner_body_relation_)),
	damping_random_(DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>(&inner_body_relation_, 0.1, "Velocity", physical_viscosity))
{}

void expandBoundingBox(BoundingBox* original, BoundingBox* additional)
{
	for(int i = 0; i < original->first.size(); i++)
	{
		if ( additional->first[i] < original->first[i] )
		{
			original->first[i] = additional->first[i];
		}
		if ( additional->second[i] > original->second[i] )
		{
			original->second[i] = additional->second[i];
		}
	}
}

void relaxParticlesSingleResolution(In_Output* in_output,
									bool write_particle_relaxation_data,
									ImportedModel* imported_model,
									ElasticSolidParticles* imported_model_particles,
									BodyRelationInner* imported_model_inner)
{	

	BodyStatesRecordingToVtu write_imported_model_to_vtu(*in_output, { imported_model });
	MeshRecordingToPlt cell_linked_list_recording(*in_output, imported_model, imported_model->cell_linked_list_);

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::SolidRelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	if (write_particle_relaxation_data)
	{
		write_imported_model_to_vtu.writeToFile(0.0);
	} 
	imported_model->updateCellLinkedList();
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
				write_imported_model_to_vtu.writeToFile(ite_p);
			}
		}
	}
	std::cout << "The physics relaxation process of the imported model finished !" << std::endl;
}

StructuralSimulationInput::StructuralSimulationInput(
	string relative_input_path,
	vector<string> imported_stl_list,
	Real scale_stl,
	vector<Vec3d> translation_list,
	vector<Real> resolution_list,
	vector<LinearElasticSolid> material_model_list,
	Real physical_viscosity,
	vector<array<int, 2>> contacting_bodies_list
	):
	relative_input_path_(relative_input_path),
	imported_stl_list_(imported_stl_list),
	scale_stl_(scale_stl),
	translation_list_(translation_list),
	resolution_list_(resolution_list),
	material_model_list_(material_model_list),
	physical_viscosity_(physical_viscosity),
	contacting_body_pairs_list_(contacting_bodies_list)
{
	//time dependent contact
	time_dep_contacting_body_pairs_list_ = {};
	// particle_relaxation option
	particle_relaxation_list_ = {};
	for (size_t i = 0; i < resolution_list_.size(); i++){ particle_relaxation_list_.push_back(true); }
	write_particle_relaxation_data_ = false;
	// scale system boundaries
	scale_system_boundaries_ = 1;
	// boundary conditions
	non_zero_gravity_ = {};
	acceleration_bounding_box_tuple_ = {};
	spring_damper_tuple_ = {};
	body_indeces_fixed_constraint_ = {};
	position_solid_body_tuple_ = {};
	position_scale_solid_body_tuple_ = {};
	translation_solid_body_tuple_ = {};
	translation_solid_body_part_tuple_ = {};
};

///////////////////////////////////////
/* StructuralSimulation members */
///////////////////////////////////////

StructuralSimulation::StructuralSimulation(StructuralSimulationInput& input):
	// generic input
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
	spring_damper_tuple_(input.spring_damper_tuple_),
	body_indeces_fixed_constraint_(input.body_indeces_fixed_constraint_),
	position_solid_body_tuple_(input.position_solid_body_tuple_),
	position_scale_solid_body_tuple_(input.position_scale_solid_body_tuple_),
	translation_solid_body_tuple_(input.translation_solid_body_tuple_),
	translation_solid_body_part_tuple_(input.translation_solid_body_part_tuple_),

	// iterators
	iteration_(0),

	// data storage
	von_mises_stress_max_({})

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
	initializeSpringDamperConstraintParticleWise();
	initializeConstrainSolidBodyRegion();
	initializePositionSolidBody();
	initializePositionScaleSolidBody();
	initializeTranslateSolidBody();
	initializeTranslateSolidBodyPart();

	// initialize simulation
	initializeSimulation();
}

StructuralSimulation::~StructuralSimulation()
{}

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
		BoundingBox additional = body_mesh_list_[i].findBounds();
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
		body_mesh_list_.push_back(TriangleMeshShape(relative_input_path_copy.append(imported_stl_list_[i]), translation_list_[i], scale_stl_));
	}
}

void StructuralSimulation::createParticleAdaptationList()
{
	particle_adaptation_list_ = {};
	for (size_t i = 0; i < resolution_list_.size(); i++)
	{
		Real system_resolution_ratio = system_resolution_ / resolution_list_[i];
		// for solid bodies, slightly small h_spaing_ratio is used
		particle_adaptation_list_.push_back(ParticleAdaptation(1.15, system_resolution_ratio));
	}
}

void StructuralSimulation::initializeElasticSolidBodies()
{
	solid_body_list_ = {};
	for (size_t i = 0; i < body_mesh_list_.size(); i++)
	{
		solid_body_list_.emplace_back(make_shared<SolidBodyForSimulation>(system_, imported_stl_list_[i], body_mesh_list_[i], particle_adaptation_list_[i], physical_viscosity_, material_model_list_[i]));
		if (particle_relaxation_list_[i])
		{
			relaxParticlesSingleResolution(&in_output_, write_particle_relaxation_data_, solid_body_list_[i]->getImportedModel(), solid_body_list_[i]->getElasticSolidParticles(), solid_body_list_[i]->getInnerBodyRelation());
		}
	}
}

void StructuralSimulation::initializeContactBetweenTwoBodies(int first, int second)
{	
	ImportedModel* first_body = solid_body_list_[first]->getImportedModel();
	ImportedModel* second_body = solid_body_list_[second]->getImportedModel();

	SolidBodyRelationContact* first_contact = new SolidBodyRelationContact(first_body, {second_body});
	SolidBodyRelationContact* second_contact = new SolidBodyRelationContact(second_body, {first_body});

	contact_list_.emplace_back(first_contact);
	contact_list_.emplace_back(second_contact);

	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(first_contact));
	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(second_contact));

	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(first_contact));
	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(second_contact));
}

void StructuralSimulation::initializeAllContacts()
{
	contact_list_ = {};
	contact_density_list_ = {};
	contact_force_list_ = {};
	for (size_t i = 0; i < contacting_body_pairs_list_.size(); i++)
	{
		initializeContactBetweenTwoBodies(contacting_body_pairs_list_[i][0], contacting_body_pairs_list_[i][1]);
	}
	for (size_t i = 0; i < time_dep_contacting_body_pairs_list_.size(); i++)
	{	
		int body_1 = time_dep_contacting_body_pairs_list_[i].first[0];
		int body_2 = time_dep_contacting_body_pairs_list_[i].first[1];
		initializeContactBetweenTwoBodies(body_1, body_2); //vector with first element being array with indices
	}
}

void StructuralSimulation::initializeGravity()
{
	// collect all the body indeces with non-zero gravity
	vector<int> gravity_indeces = {};
	for (size_t i = 0; i < non_zero_gravity_.size(); i++)
	{
		gravity_indeces.push_back(non_zero_gravity_[i].first);
	}
	// initialize gravity
	initialize_gravity_ = {};
	size_t gravity_index_i = 0; // iterating through gravity_indeces
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{	
		// check if i is in indeces_gravity
		if ( count(gravity_indeces.begin(), gravity_indeces.end(), i) )
		{	
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(solid_body_list_[i]->getImportedModel(), new Gravity(non_zero_gravity_[gravity_index_i].second)));
			gravity_index_i++;
		}
		else
		{
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(solid_body_list_[i]->getImportedModel()));
		}
	}
}

void StructuralSimulation::initializeAccelerationForBodyPartInBoundingBox()
{	
	acceleration_bounding_box_ = {};
	for (size_t i = 0; i < acceleration_bounding_box_tuple_.size(); i++)
	{
		SolidBody* solid_body = solid_body_list_[get<0>(acceleration_bounding_box_tuple_[i])]->getImportedModel();
		acceleration_bounding_box_.emplace_back(make_shared<solid_dynamics::AccelerationForBodyPartInBoundingBox>
			(solid_body, get<1>(acceleration_bounding_box_tuple_[i]), get<2>(acceleration_bounding_box_tuple_[i])));
    }
}

void StructuralSimulation::initializeSpringDamperConstraintParticleWise()
{	
	spring_damper_constraint_ = {};
	for (size_t i = 0; i < spring_damper_tuple_.size(); i++)
	{
		SolidBody* solid_body = solid_body_list_[get<0>(spring_damper_tuple_[i])]->getImportedModel();
		spring_damper_constraint_.emplace_back(make_shared<solid_dynamics::SpringDamperConstraintParticleWise>(solid_body, get<1>(spring_damper_tuple_[i]), get<2>(spring_damper_tuple_[i])));
    }
}

void StructuralSimulation::initializeConstrainSolidBodyRegion()
{	
	fixed_constraint_ = {};
	for (size_t i = 0; i < body_indeces_fixed_constraint_.size(); i++)
	{
		int body_index = body_indeces_fixed_constraint_[i];
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->getImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
		fixed_constraint_.emplace_back(make_shared<solid_dynamics::ConstrainSolidBodyRegion>(solid_body_list_[body_index]->getImportedModel(), bp));
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
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->getImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
			
		position_solid_body_.emplace_back(make_shared<solid_dynamics::PositionSolidBody>(solid_body_list_[body_index]->getImportedModel(), bp, start_time, end_time, pos_end_center));
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
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->getImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
				
		position_scale_solid_body_.emplace_back(make_shared<solid_dynamics::PositionScaleSolidBody>(solid_body_list_[body_index]->getImportedModel(), bp, start_time, end_time, scale));
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
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(
			solid_body_list_[body_index]->getImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
			
		translation_solid_body_.emplace_back(make_shared<solid_dynamics::TranslateSolidBody>(
			solid_body_list_[body_index]->getImportedModel(), bp, start_time, end_time, translation));
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
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(
			solid_body_list_[body_index]->getImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
			
		translation_solid_body_part_.emplace_back(make_shared<solid_dynamics::TranslateSolidBodyPart>(
			solid_body_list_[body_index]->getImportedModel(), bp, start_time, end_time, translation, bbox));
	}
}

void StructuralSimulation::executeCorrectConfiguration()
{
	for (size_t i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->getCorrectConfiguration()->parallel_exec();
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

void StructuralSimulation::executeSpringDamperConstraintParticleWise()
{
	for (size_t i = 0; i < spring_damper_constraint_.size(); i++)
	{
		spring_damper_constraint_[i]->parallel_exec();
	}
}

void StructuralSimulation::executeContactDensitySummation()
{
	// number of contacts that are not time dependent: contact pairs * 2
	size_t number_of_general_contacts = contacting_body_pairs_list_.size() * 2;
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
			if(GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
			{
				contact_density_list_[i]->parallel_exec();
			}
		}
	}
}

void StructuralSimulation::executeContactForce()
{
	// number of contacts that are not time dependent: contact pairs * 2
	size_t number_of_general_contacts = contacting_body_pairs_list_.size() * 2;
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
			if(GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
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

void StructuralSimulation::executeConstrainSolidBodyRegion()
{
	for (size_t i = 0; i < fixed_constraint_.size(); i++)
	{
		fixed_constraint_[i]->parallel_exec();
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
		solid_body_list_[i]->getImportedModel()->updateCellLinkedList();
	}
}

void StructuralSimulation::executeContactUpdateConfiguration()
{	
	// number of contacts that are not time dependent: contact pairs * 2
	size_t number_of_general_contacts = contacting_body_pairs_list_.size() * 2;
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
			if(GlobalStaticVariables::physical_time_ >= start_time && GlobalStaticVariables::physical_time_ <= end_time)
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
	if (iteration_ % 100 == 0) cout << "N=" << iteration_ << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";

	/** ACTIVE BOUNDARY CONDITIONS */
	// force (acceleration) based
	executeinitializeATimeStep();
	executeAccelerationForBodyPartInBoundingBox();
	executeSpringDamperConstraintParticleWise();

	/** CONTACT */
	executeContactDensitySummation();
	executeContactForce();

	/** STRESS RELAXATOIN, DAMPING, POSITIONAL CONSTRAINTS */
	executeStressRelaxationFirstHalf(dt);

	executeConstrainSolidBodyRegion();
	executePositionSolidBody(dt);
	executePositionScaleSolidBody(dt);
	executeTranslateSolidBody(dt);
	// velocity based
	executeTranslateSolidBodyPart(dt);

	executeDamping(dt);

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
	BodyStatesRecordingToVtu write_states(in_output_, system_.real_bodies_);

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
		von_mises_stress_max_.push_back(solid_body_list_[0].get()->getElasticSolidParticles()->getMaxVonMisesStress());
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
	BodyStatesRecordingToVtu write_states(in_output_, system_.real_bodies_);
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