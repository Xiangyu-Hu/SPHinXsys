#include "solid_structural_simulation_class.h"

////////////////////////////////////////////////////
/* global functions in SolidStructuralSimulation  */
////////////////////////////////////////////////////

ImportedModel::ImportedModel(SPHSystem &system, std::string body_name, TriangleMeshShape* triangle_mesh_shape, Real resolution)
    : SolidBody(system, body_name, resolution)
{
    ComplexShape original_body_shape;
    original_body_shape.addTriangleMeshShape(triangle_mesh_shape, ShapeBooleanOps::add);
    body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
}

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

	WriteBodyStatesToVtu write_imported_model_to_vtu(*in_output, { imported_model });
	WriteMeshToPlt write_mesh_cell_linked_list(*in_output, imported_model, imported_model->mesh_cell_linked_list_);

	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	if (write_particles_to_file)
	{
		write_imported_model_to_vtu.WriteToFile(0.0);
	} 
	imported_model->updateCellLinkedList();
	if (write_particles_to_file)
	{
		write_mesh_cell_linked_list.WriteToFile(0.0);
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
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			if (write_particles_to_file)
			{
				write_imported_model_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
			}
		}
	}
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
}

///////////////////////////////////////
/* SolidStructuralSimulation members */
///////////////////////////////////////

void SolidStructuralSimulation::ImportSTLModels()
{
	for (auto imported_stl: *imported_stl_list_)
	{	
		std::string relative_input_path_copy = relative_input_path_;
		body_mesh_list_.push_back(new TriangleMeshShape(relative_input_path_copy.append(imported_stl), Vecd(0, 0, 0), scale_stl_));
	}
}

BoundingBox* SolidStructuralSimulation::CalculateSystemBoundaries()
{
	BoundingBox* system_domain_bounds = new BoundingBox(Vecd(0, 0, 0), Vecd(0, 0, 0));
	
	for (auto body_mesh: body_mesh_list_)
	{
		BoundingBox additional = body_mesh->findBounds();
		ExpandBoundingBox(system_domain_bounds, &additional);
	}
	return system_domain_bounds;
}

void SolidStructuralSimulation::SetupSystem()
{
	system_ = new SPHSystem (*CalculateSystemBoundaries(), default_resolution_);
	system_->run_particle_relaxation_ = true;
	in_output_ = new In_Output (*system_);
};

void SolidStructuralSimulation::InitializeElasticBodies()
{	
	int i = 0;
	for (auto body_mesh: body_mesh_list_)
	{	
		std::string imported_stl = (*imported_stl_list_)[i];
		ImportedModel* imported_model = new ImportedModel(*system_, imported_stl, body_mesh, (*resolution_list_)[i]);
		imported_model_list_.push_back(imported_model);

		ElasticSolidParticles* imported_model_particles = new ElasticSolidParticles(imported_model, (*material_model_list_)[i]);
		imported_model_particles_list_.push_back(imported_model_particles);

		InnerBodyRelation* imported_model_inner = new InnerBodyRelation(imported_model);
		imported_model_inner_list_.push_back(imported_model_inner);
		RelaxParticlesSingleResolution(in_output_, false, imported_model, imported_model_particles, imported_model_inner);

		correct_configuration_list_.push_back(new solid_dynamics::CorrectConfiguration(imported_model_inner));
		stress_relaxation_first_half_list_.push_back(new solid_dynamics::StressRelaxationFirstHalf(imported_model_inner));
		stress_relaxation_second_half_list_.push_back(new solid_dynamics::StressRelaxationSecondHalf(imported_model_inner));

		damping_list_.push_back(new DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>(imported_model_inner, 0.1, "Velocity", physical_viscosity_));

		i++;
	}
}

void SolidStructuralSimulation::InitializeContactBetweenTwoBodies(int first, int second)
{	
	ImportedModel* first_body = imported_model_list_[first];
	ImportedModel* second_body = imported_model_list_[second];
	SolidContactBodyRelation* first_contact = new SolidContactBodyRelation(first_body, {second_body});
	SolidContactBodyRelation* second_contact = new SolidContactBodyRelation(second_body, {first_body});
	contact_list_.push_back(first_contact);
	contact_list_.push_back(second_contact);

	contact_density_list_.push_back(new solid_dynamics::ContactDensitySummation (first_contact));
	contact_density_list_.push_back(new solid_dynamics::ContactDensitySummation (second_contact));

	contact_force_list_.push_back(new solid_dynamics::ContactForce (first_contact));
	contact_force_list_.push_back(new solid_dynamics::ContactForce (second_contact));
}

void SolidStructuralSimulation::InitializeGravity()
{
	int i = 0;
	for (auto imported_model: imported_model_list_)
	{	
		if ( std::find(body_indeces_gravity_.begin(), body_indeces_gravity_.end(), i) != body_indeces_gravity_.end() )
		{	
			initialize_gravity_.push_back(new InitializeATimeStep(imported_model, new Gravity(*gravity_[i])));
		}
		else
		{
			initialize_gravity_.push_back(new InitializeATimeStep(imported_model));
		}
		i++;
	}
}

void SolidStructuralSimulation::AddGravity(int body_index, Vecd* gravity)
{
	body_indeces_gravity_.push_back(body_index);
	gravity_.push_back(gravity);
}

void SolidStructuralSimulation::InitializeAccelerationForBodyPartInBoundingBox()
{	
	int i = 0;
	for (auto body_index: body_indeces_accelerations_)
	{
		acceleration_for_body_part_.push_back(new solid_dynamics::AccelerationForBodyPartInBoundingBox(imported_model_list_[body_index], bounding_boxes_[i], accelerations_[i]));
	}
	i++;
}

void SolidStructuralSimulation::AddAccelerationForBodyPartInBoundingBox(int body_index, BoundingBox* bounding_box, Vecd acceleration)
{
	body_indeces_accelerations_.push_back(body_index);
	bounding_boxes_.push_back(bounding_box);
	accelerations_.push_back(acceleration);
}

void SolidStructuralSimulation::InitializeSpringDamperConstraintParticleWise()
{	
	int i = 0;
	for (auto body_index: body_indeces_spring_damper_)
	{
		spring_damper_contraint_.push_back(new solid_dynamics::SpringDamperConstraintParticleWise(imported_model_list_[body_index], stiffnesses_[i], damping_ratios_[i]));
	}
	i++;
}

void SolidStructuralSimulation::AddSpringDamperConstraintParticleWise(int body_index, Vecd stiffness, Real damping_ratio)
{
	body_indeces_spring_damper_.push_back(body_index);
	stiffnesses_.push_back(stiffness);
	damping_ratios_.push_back(damping_ratio);
}

void SolidStructuralSimulation::InitializeConstrainSolidBodyRegion()
{	
	for (auto body_index: body_indeces_fixed_contraint_)
	{
		BodyPartByParticle* bp = new BodyPartByParticle(imported_model_list_[body_index], (*imported_stl_list_)[body_index], body_mesh_list_[body_index]);
		fixed_contraint_.push_back(new solid_dynamics::ConstrainSolidBodyRegion(imported_model_list_[body_index], bp));
	}
}

void SolidStructuralSimulation::AddConstrainSolidBodyRegion(int body_index)
{
	body_indeces_fixed_contraint_.push_back(body_index);
}

void SolidStructuralSimulation::ExecuteCorrectConfiguration()
{
	for (auto correct_configuration: correct_configuration_list_)
	{
		correct_configuration->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteInitializeATimeStep()
{
	for (auto ig: initialize_gravity_)
	{
		ig->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteAccelerationForBodyPartInBoundingBox()
{
	for (auto acc: acceleration_for_body_part_)
	{
		acc->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteSpringDamperConstraintParticleWise()
{
	for (auto sd: spring_damper_contraint_)
	{
		sd->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteContactDensitySummation()
{
	for (auto cd: contact_density_list_)
	{
		cd->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteContactForce()
{
	for (auto cf: contact_force_list_)
	{
		cf->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteStressRelaxationFirstHalf(Real dt)
{
	for (auto sr: stress_relaxation_first_half_list_)
	{
		sr->parallel_exec(dt);
	}
}

void SolidStructuralSimulation::ExecuteConstrainSolidBodyRegion()
{
	for (auto fc: fixed_contraint_)
	{
		fc->parallel_exec();
	}
}

void SolidStructuralSimulation::ExecuteDamping(Real dt)
{
	for (auto da: damping_list_)
	{
		da->parallel_exec(dt);
	}
}

void SolidStructuralSimulation::ExecuteStressRelaxationSecondHalf(Real dt)
{
	for (auto sr: stress_relaxation_second_half_list_)
	{
		sr->parallel_exec(dt);
	}
}

void SolidStructuralSimulation::ExecuteUpdateCellLinkedList()
{
	for (auto im: imported_model_list_)
	{
		im->updateCellLinkedList();
	}
}

void SolidStructuralSimulation::ExecuteContactUpdateConfiguration()
{
	for (auto cl: contact_list_)
	{
		cl->updateConfiguration();
	}
}

void SolidStructuralSimulation::RunSimulationStep(int &ite, Real &dt, Real &integration_time)
{
	if (ite % 100 == 0) cout << "N=" << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";

	/** ACTIVE BOUNDARY CONDITIONS */
	ExecuteInitializeATimeStep();
	ExecuteAccelerationForBodyPartInBoundingBox();
	ExecuteSpringDamperConstraintParticleWise();

	/** CONTACT */
	ExecuteContactDensitySummation();
	ExecuteContactForce();

	/** STRESS RELAXATOIN, DAMPING, PASSIVE CONSTRAINTS */
	ExecuteStressRelaxationFirstHalf(dt);
	ExecuteConstrainSolidBodyRegion();
	ExecuteDamping(dt);
	ExecuteConstrainSolidBodyRegion();
	ExecuteStressRelaxationSecondHalf(dt);
	
	/** UPDATE TIME STEP SIZE, INCREMENT */
	ite++;
	dt = system_->getSmallestTimeStepAmongSolidBodies();
	integration_time += dt;
	GlobalStaticVariables::physical_time_ += dt;
	
	/** UPDATE BODIES CELL LINKED LISTS */
	ExecuteUpdateCellLinkedList();
	
	/** UPDATE CONTACT CONFIGURATION */
	ExecuteContactUpdateConfiguration();
}

void SolidStructuralSimulation::RunSimulation(Real end_time)
{
	WriteBodyStatesToVtu write_states(*in_output_, system_->real_bodies_);
	GlobalStaticVariables::physical_time_ = 0.0;
	
	/** INITIALALIZE SYSTEM */
	system_->initializeSystemCellLinkedLists();
	system_->initializeSystemConfigurations();

	/** INITIAL CONDITION */
	ExecuteCorrectConfiguration();
	
	/** Statistics for computing time. */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
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
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
}
