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
	inner_body_relation_(InnerBodyRelation(&imported_model_)),

	correct_configuration_(solid_dynamics::CorrectConfiguration(&inner_body_relation_)),
	stress_relaxation_first_half_(solid_dynamics::StressRelaxationFirstHalf(&inner_body_relation_)),
	stress_relaxation_second_half_(solid_dynamics::StressRelaxationSecondHalf(&inner_body_relation_)),
	damping_random_(DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>(&inner_body_relation_, 0.1, "Velocity", physical_viscosity))
{}

void ExpandBoundingBox(BoundingBox* original, BoundingBox* additional)
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

void RelaxParticlesSingleResolution(In_Output* in_output,
									bool write_particles_to_file,
									ImportedModel* imported_model,
									ElasticSolidParticles* imported_model_particles,
									InnerBodyRelation* imported_model_inner)
{	

	BodyStatesRecordingToVtu write_imported_model_to_vtu(*in_output, { imported_model });
	MeshRecordingToPlt mesh_cell_linked_list_recording(*in_output, imported_model, imported_model->mesh_cell_linked_list_);

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
	if (write_particles_to_file)
	{
		write_imported_model_to_vtu.writeToFile(0.0);
	} 
	imported_model->updateCellLinkedList();
	if (write_particles_to_file)
	{
		mesh_cell_linked_list_recording.writeToFile(0.0);
	}
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 500)
	{
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			cout << fixed << setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			if (write_particles_to_file)
			{
				write_imported_model_to_vtu.writeToFile(Real(ite_p) * 1.0e-4);
			}
		}
	}
	cout << "The physics relaxation process of imported model finish !" << endl;
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
	contacting_bodies_list_(contacting_bodies_list)
{
	// particle_relaxation option
	particle_relaxation_list_ = {};
	for (unsigned i = 0; i < resolution_list_.size(); i++){ particle_relaxation_list_.push_back(true); }
	// scale system boundaries
	scale_system_boundaries_ = 1;
	// boundary conditions
	non_zero_gravity_ = {};
	acceleration_bounding_box_tuple_ = {};
	spring_damper_tuple_ = {};
	body_indeces_fixed_constraint_ = {};
	position_solid_body_tuple_ = {};
	position_scale_solid_body_tuple_ = {};
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
	contacting_bodies_list_(input.contacting_bodies_list_),

	// default system, optional: particle relaxation, scale_system_boundaries
	particle_relaxation_list_(input.particle_relaxation_list_),
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
	position_scale_solid_body_tuple_(input.position_scale_solid_body_tuple_)
{
	// scaling of translation and resolution
	ScaleTranslationAndResolution();
	// set the default resolution to the max in the resolution list
	SetSystemResolutionMax();
	// create the body mesh list for triangular mesh shapes storage
	CreateBodyMeshList();
	// create the particle adaptions for the bodies
	CreateParticleAdaptationList();
	// set up the system
	CalculateSystemBoundaries();
	system_.run_particle_relaxation_ = true;
	// initialize solid bodies with their properties
	InitializeElasticSolidBodies();
	// contacts
	InitializeAllContacts();

	// boundary conditions
	InitializeGravity();
	InitializeAccelerationForBodyPartInBoundingBox();
	InitializeSpringDamperConstraintParticleWise();
	InitializeConstrainSolidBodyRegion();
	InitializePositionSolidBody();
	InitializePositionScaleSolidBody();
	InitializeTranslateSolidBody();

	// initialize simulation
	InitSimulation();
}

StructuralSimulation::~StructuralSimulation()
{}

void StructuralSimulation::ScaleTranslationAndResolution()
{
	// scale the translation_list_, system_resolution_ and resolution_list_
	for (unsigned int i = 0; i < translation_list_.size(); i++)
	{
		translation_list_[i] *= scale_stl_;
	}
	system_resolution_ *= scale_stl_;
	for (unsigned int i = 0; i < resolution_list_.size(); i++)
	{
		resolution_list_[i] *= scale_stl_;
	}
}

void StructuralSimulation::SetSystemResolutionMax()
{
	system_resolution_ = 0.0;
	for (unsigned int i = 0; i < resolution_list_.size(); i++)
	{
		if (system_resolution_ < resolution_list_[i])
		{
			system_resolution_ = resolution_list_[i];
		}
	}
	system_.resolution_ref_ = system_resolution_;
}

void StructuralSimulation::CalculateSystemBoundaries()
{	
	// calculate system bounds from all bodies
	for (unsigned int i = 0; i < body_mesh_list_.size(); i++)
	{
		BoundingBox additional = body_mesh_list_[i].findBounds();
		ExpandBoundingBox(&system_.system_domain_bounds_, &additional);
	}
	// scale the system bounds around the center point
	Vecd center_point = (system_.system_domain_bounds_.first + system_.system_domain_bounds_.second) * 0.5;

	Vecd distance_first = system_.system_domain_bounds_.first - center_point;
	Vecd distance_second = system_.system_domain_bounds_.second - center_point;

	system_.system_domain_bounds_.first = center_point + distance_first * scale_system_boundaries_;
	system_.system_domain_bounds_.second = center_point + distance_second * scale_system_boundaries_;
}

void StructuralSimulation::CreateBodyMeshList()
{
	body_mesh_list_ = {};
	for (unsigned int i = 0; i < imported_stl_list_.size(); i++)
	{
		string relative_input_path_copy = relative_input_path_;
		body_mesh_list_.push_back(TriangleMeshShape(relative_input_path_copy.append(imported_stl_list_[i]), translation_list_[i], scale_stl_));
	}
}

void StructuralSimulation::CreateParticleAdaptationList()
{
	particle_adaptation_list_ = {};
	for (unsigned int i = 0; i < resolution_list_.size(); i++)
	{
		particle_adaptation_list_.push_back(ParticleAdaptation());

		Real system_resolution_ratio = resolution_list_[i] / system_resolution_;
		particle_adaptation_list_[i].SetSystemResolutionRatio(system_resolution_ratio);
	}
}

void StructuralSimulation::InitializeElasticSolidBodies()
{
	solid_body_list_ = {};
	for (unsigned int i = 0; i < body_mesh_list_.size(); i++)
	{
		solid_body_list_.emplace_back(make_shared<SolidBodyForSimulation>(system_, imported_stl_list_[i], body_mesh_list_[i], particle_adaptation_list_[i], physical_viscosity_, material_model_list_[i]));
		if (particle_relaxation_list_[i])
		{
			RelaxParticlesSingleResolution(&in_output_, false, solid_body_list_[i]->GetImportedModel(), solid_body_list_[i]->GetElasticSolidParticles(), solid_body_list_[i]->GetInnerBodyRelation());
		}
	}
}

void StructuralSimulation::InitializeContactBetweenTwoBodies(int first, int second)
{	
	ImportedModel* first_body = solid_body_list_[first]->GetImportedModel();
	ImportedModel* second_body = solid_body_list_[second]->GetImportedModel();

	SolidContactBodyRelation* first_contact = new SolidContactBodyRelation(first_body, {second_body});
	SolidContactBodyRelation* second_contact = new SolidContactBodyRelation(second_body, {first_body});

	contact_list_.emplace_back(first_contact);
	contact_list_.emplace_back(second_contact);

	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(first_contact));
	contact_density_list_.emplace_back(make_shared<solid_dynamics::ContactDensitySummation>(second_contact));

	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(first_contact));
	contact_force_list_.emplace_back(make_shared<solid_dynamics::ContactForce>(second_contact));
}

void StructuralSimulation::InitializeAllContacts()
{
	contact_list_ = {};
	contact_density_list_ = {};
	contact_force_list_ = {};
	for (unsigned int i = 0; i < contacting_bodies_list_.size(); i++)
	{
		InitializeContactBetweenTwoBodies(contacting_bodies_list_[i][0], contacting_bodies_list_[i][1]);
	}
}

void StructuralSimulation::InitializeGravity()
{
	// collect all the body indeces with non-zero gravity
	vector<int> body_indeces_gravity = {};
	for (unsigned int i = 0; i < non_zero_gravity_.size(); i++)
	{
		body_indeces_gravity.push_back(non_zero_gravity_[i].first);
	}
	// initialize gravity
	initialize_gravity_ = {};
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{	
		if ( find(body_indeces_gravity.begin(), body_indeces_gravity.end(), i) != body_indeces_gravity.end() )
		{	
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(solid_body_list_[i]->GetImportedModel(), new Gravity(non_zero_gravity_[i].second)));
		}
		else
		{
			initialize_gravity_.emplace_back(make_shared<TimeStepInitialization>(solid_body_list_[i]->GetImportedModel()));
		}
	}
}

void StructuralSimulation::InitializeAccelerationForBodyPartInBoundingBox()
{	
	acceleration_bounding_box_ = {};
	for (unsigned int i = 0; i < acceleration_bounding_box_tuple_.size(); i++)
	{
		SolidBody* solid_body = solid_body_list_[get<0>(acceleration_bounding_box_tuple_[i])]->GetImportedModel();
		acceleration_bounding_box_.emplace_back(make_shared<solid_dynamics::AccelerationForBodyPartInBoundingBox>
			(solid_body, &get<1>(acceleration_bounding_box_tuple_[i]), get<2>(acceleration_bounding_box_tuple_[i])));
    }
}

void StructuralSimulation::InitializeSpringDamperConstraintParticleWise()
{	
	spring_damper_constraint_ = {};
	for (unsigned int i = 0; i < spring_damper_tuple_.size(); i++)
	{
		SolidBody* solid_body = solid_body_list_[get<0>(spring_damper_tuple_[i])]->GetImportedModel();
		spring_damper_constraint_.emplace_back(make_shared<solid_dynamics::SpringDamperConstraintParticleWise>(solid_body, get<1>(spring_damper_tuple_[i]), get<2>(spring_damper_tuple_[i])));
    }
}

void StructuralSimulation::InitializeConstrainSolidBodyRegion()
{	
	fixed_constraint_ = {};
	for (unsigned int i = 0; i < body_indeces_fixed_constraint_.size(); i++)
	{
		int body_index = body_indeces_fixed_constraint_[i];
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->GetImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
		fixed_constraint_.emplace_back(make_shared<solid_dynamics::ConstrainSolidBodyRegion>(solid_body_list_[body_index]->GetImportedModel(), bp));
	}
}

void StructuralSimulation::InitializePositionSolidBody()
{
	position_solid_body_ = {};
	for (unsigned int i = 0; i < position_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(position_solid_body_tuple_[i]);
		Real start_time = get<1>(position_solid_body_tuple_[i]);
		Real end_time = get<2>(position_solid_body_tuple_[i]);
		Vecd pos_end_center = get<3>(position_solid_body_tuple_[i]);
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->GetImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
			
		position_solid_body_.emplace_back(make_shared<solid_dynamics::PositionSolidBody>(solid_body_list_[body_index]->GetImportedModel(), bp, start_time, end_time, pos_end_center));
	}
}

void StructuralSimulation::InitializePositionScaleSolidBody()
{
	position_scale_solid_body_ = {};
	for (unsigned int i = 0; i < position_scale_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(position_scale_solid_body_tuple_[i]);
		Real start_time = get<1>(position_scale_solid_body_tuple_[i]);
		Real end_time = get<2>(position_scale_solid_body_tuple_[i]);
		Real scale = get<3>(position_scale_solid_body_tuple_[i]);
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->GetImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
				
		position_scale_solid_body_.emplace_back(make_shared<solid_dynamics::PositionScaleSolidBody>(solid_body_list_[body_index]->GetImportedModel(), bp, start_time, end_time, scale));
	}
}

void StructuralSimulation::InitializeTranslateSolidBody()
{
	translation_solid_body_ = {};
	for (unsigned int i = 0; i < translation_solid_body_tuple_.size(); i++)
	{
		int body_index = get<0>(translation_solid_body_tuple_[i]);
		Real start_time = get<1>(translation_solid_body_tuple_[i]);
		Real end_time = get<2>(translation_solid_body_tuple_[i]);
		Vecd translation = get<3>(translation_solid_body_tuple_[i]);
		BodyPartByParticleTriMesh* bp = new BodyPartByParticleTriMesh(solid_body_list_[body_index]->GetImportedModel(), imported_stl_list_[body_index], &body_mesh_list_[body_index]);
			
		translation_solid_body_.emplace_back(make_shared<solid_dynamics::TranslateSolidBody>(solid_body_list_[body_index]->GetImportedModel(), bp, start_time, end_time, translation));
	}
}

void StructuralSimulation::ExecuteCorrectConfiguration()
{
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->GetCorrectConfiguration()->parallel_exec();
	}
}

void StructuralSimulation::ExecuteInitializeATimeStep()
{
	for (unsigned int i = 0; i < initialize_gravity_.size(); i++)
	{
		initialize_gravity_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecuteAccelerationForBodyPartInBoundingBox()
{
	for (unsigned int i = 0; i < acceleration_bounding_box_.size(); i++)
	{
		acceleration_bounding_box_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecuteSpringDamperConstraintParticleWise()
{
	for (unsigned int i = 0; i < spring_damper_constraint_.size(); i++)
	{
		spring_damper_constraint_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecuteContactDensitySummation()
{
	for (unsigned int i = 0; i < contact_density_list_.size(); i++)
	{
		contact_density_list_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecuteContactForce()
{
	for (unsigned int i = 0; i < contact_force_list_.size(); i++)
	{
		contact_force_list_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecuteStressRelaxationFirstHalf(Real dt)
{
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->GetStressRelaxationFirstHalf()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteConstrainSolidBodyRegion()
{
	for (unsigned int i = 0; i < fixed_constraint_.size(); i++)
	{
		fixed_constraint_[i]->parallel_exec();
	}
}

void StructuralSimulation::ExecutePositionSolidBody(Real dt)
{
	for (unsigned int i = 0; i < position_solid_body_.size(); i++)
	{
		position_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecutePositionScaleSolidBody(Real dt)
{
	for (unsigned int i = 0; i < position_scale_solid_body_.size(); i++)
	{
		position_scale_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteTranslateSolidBody(Real dt)
{
	for (unsigned int i = 0; i < translation_solid_body_.size(); i++)
	{
		translation_solid_body_[i]->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteDamping(Real dt)
{
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->GetDampingWithRandomChoice()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteStressRelaxationSecondHalf(Real dt)
{
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->GetStressRelaxationSecondHalf()->parallel_exec(dt);
	}
}

void StructuralSimulation::ExecuteUpdateCellLinkedList()
{
	for (unsigned int i = 0; i < solid_body_list_.size(); i++)
	{
		solid_body_list_[i]->GetImportedModel()->updateCellLinkedList();
	}
}

void StructuralSimulation::ExecuteContactUpdateConfiguration()
{
	for (unsigned int i = 0; i < contact_list_.size(); i++)
	{
		contact_list_[i]->updateConfiguration();
	}
}

void StructuralSimulation::RunSimulationStep(int &ite, Real &dt, Real &integration_time)
{
	if (ite % 100 == 0) cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";

	/** ACTIVE BOUNDARY CONDITIONS */
	ExecuteInitializeATimeStep();
	ExecuteAccelerationForBodyPartInBoundingBox();
	ExecuteSpringDamperConstraintParticleWise();

	/** CONTACT */
	ExecuteContactDensitySummation();
	ExecuteContactForce();

	/** STRESS RELAXATOIN, DAMPING, POSITIONAL CONSTRAINTS */
	ExecuteStressRelaxationFirstHalf(dt);

	ExecuteConstrainSolidBodyRegion();
	ExecutePositionSolidBody(dt);
	ExecutePositionScaleSolidBody(dt);
	ExecuteTranslateSolidBody(dt);

	ExecuteDamping(dt);

	ExecuteConstrainSolidBodyRegion();
	ExecutePositionSolidBody(dt);
	ExecutePositionScaleSolidBody(dt);
	ExecuteTranslateSolidBody(dt);

	ExecuteStressRelaxationSecondHalf(dt);
	
	/** UPDATE TIME STEP SIZE, INCREMENT */
	ite++;
	dt = system_.getSmallestTimeStepAmongSolidBodies();
	integration_time += dt;
	GlobalStaticVariables::physical_time_ += dt;
	
	/** UPDATE BODIES CELL LINKED LISTS */
	ExecuteUpdateCellLinkedList();
	
	/** UPDATE CONTACT CONFIGURATION */
	ExecuteContactUpdateConfiguration();
}

void StructuralSimulation::RunSimulation(Real end_time)
{
	BodyStatesRecordingToVtu write_states(in_output_, system_.real_bodies_);
	GlobalStaticVariables::physical_time_ = 0.0;

	/** INITIALALIZE SYSTEM */
	system_.initializeSystemCellLinkedLists();
	system_.initializeSystemConfigurations();

	/** INITIAL CONDITION */
	ExecuteCorrectConfiguration();	

	/** Statistics for computing time. */
	write_states.writeToFile(GlobalStaticVariables::physical_time_);
	int ite = 0;
	Real output_period = 0.1 / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period) 
		{
			RunSimulationStep(ite, dt, integration_time);
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
}

void StructuralSimulation::InitSimulation()
{	
	/** INITIALALIZE SYSTEM */
	system_.initializeSystemCellLinkedLists();
	system_.initializeSystemConfigurations();

	/** INITIAL CONDITION */
	ExecuteCorrectConfiguration();
}

double StructuralSimulation::RunSimulationFixedDurationJS(int number_of_steps)
{
	BodyStatesRecordingToVtu write_states(in_output_, system_.real_bodies_);
	GlobalStaticVariables::physical_time_ = 0.0;
	
	/** Statistics for computing time. */
	write_states.writeToFile(GlobalStaticVariables::physical_time_);
	int output_period = 100;
	int ite = 0;	
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop */
	while (ite < number_of_steps)
	{
		Real integration_time = 0.0;
		int output_step = 0;
		while (output_step < output_period)
		{
			RunSimulationStep(ite, dt, integration_time);
			output_step++;
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	return tt.seconds();
}