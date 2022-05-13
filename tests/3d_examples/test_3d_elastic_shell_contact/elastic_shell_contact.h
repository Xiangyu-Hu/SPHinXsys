#ifndef ELASTIC_SHELL_CONTACT_H
#define ELASTIC_SHELL_CONTACT_H

#include "sphinxsys.h"
using namespace SPH;

class ShellBodyFromMesh : public ThinStructure
{
public:
    ShellBodyFromMesh(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape, SharedPtr<SPHAdaptation> particle_adaptation)
	: ThinStructure(system, body_name, particle_adaptation)
    {
        body_shape_.add<LevelSetShape>(this, triangle_mesh_shape, true, false);
        // set the body domain bounds because it is not set by default
        BoundingBox bounds = body_shape_.findBounds();
        setBodyDomainBounds(bounds);
    }
};

class SolidBodyFromMesh : public SolidBody
{
public:
	SolidBodyFromMesh(SPHSystem &system, string body_name, TriangleMeshShape& triangle_mesh_shape,
	 shared_ptr<SPHAdaptation> particle_adaptation)
	: SolidBody(system, body_name, particle_adaptation)
	{
		body_shape_.add<LevelSetShape>(this, triangle_mesh_shape, true, false);
		// set the body domain bounds because it is not set by default
		BoundingBox bounds = body_shape_.findBounds();
		setBodyDomainBounds(bounds);
	}
};

void ShellParticleRelaxation(In_Output& in_output, ShellBodyFromMesh& imported_model, BodyRelationInner& imported_model_inner, Real thickness, Real level_set_refinement_ratio)
{
	BodyStatesRecordingToVtp write_imported_model_to_vtp(in_output, { imported_model });
	MeshRecordingToPlt 	write_mesh_cell_linked_list(in_output, imported_model, imported_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(imported_model_inner, thickness, level_set_refinement_ratio);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	//write_imported_model_to_vtp.writeToFile(0.0);
	imported_model.updateCellLinkedList();
	write_mesh_cell_linked_list.writeToFile(0.0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
			//write_imported_model_to_vtp.writeToFile(ite_p);
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	relaxation_step_inner.mid_surface_bounding_.calculateNormalDirection();
	//write_imported_model_to_vtp.writeToFile(ite_p);
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;
}

class ShellDynamicsInitialCondition
	: public thin_structure_dynamics::ShellDynamicsInitialCondition
{
public:
	explicit ShellDynamicsInitialCondition(SolidBody &solid_body)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(solid_body){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial pseudo-normal. */
		n_[index_i] = n_0_[index_i];
		pseudo_n_[index_i] = n_0_[index_i];
		transformation_matrix_[index_i] = getTransformationMatrix(n_0_[index_i]);
	};
};

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

void ElasticShellContact(bool elastic_shell, Real shell_resolution, Real shell_thickness, Real solid_resolution)
{
    // Currently running it in these units
    Real unit_mm = 0.001; // mm, MPa, N, kg

	// Global parameters
	Real end_time = 2.0;
	Real scale = 1.0; // in mm currently

	Real resolution_ref = solid_resolution * scale;
	Real shell_res = shell_resolution * scale;
	Real thickness = shell_thickness * scale;
	Real res_factor = resolution_ref / shell_res;
	Real level_set_refinement_ratio = (resolution_ref / res_factor) / (0.1 * thickness);

	Real rho0_s = 1000.0 * std::pow(unit_mm, 3);
	Real poisson = 0.45;
	Real Youngs_modulus = 5e5 * std::pow(unit_mm, 2);
	Real Youngs_modulus_shell = 5e5 * std::pow(unit_mm, 2);
	Real physical_viscosity = 200.0 * std::pow(unit_mm, 2);

	Gravity gravity_solid(Vec3d(10.0, 0.0, 0.0));
	Gravity gravity_shell(Vec3d(0.0, 0.0, 0.0));

	// Import meshes
	TriangleMeshShapeSTL mesh_1("./input/solid_cube.stl", Vecd(0.0)*scale, scale);
	TriangleMeshShapeSTL mesh_2("./input/curved_shell.stl", Vecd(20.0, 0.0, 0.0)*scale, scale);

	BoundingBox system_bb = mesh_1.findBounds();
	BoundingBox bb2 = mesh_2.findBounds();
	expandBoundingBox(&system_bb, &bb2);

	// System
	SPHSystem system(system_bb, resolution_ref);
	In_Output in_output(system);

	// Create solid and shell objects
	SolidBodyFromMesh myocardium_body(system, "Solid_cube", mesh_1, makeShared<SPHAdaptation>());
	ElasticSolidParticles myocardium_particles(myocardium_body, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	BodyRelationInner myocardium_body_inner(myocardium_body);

	shared_ptr<SPHAdaptation> shell_adaptation = makeShared<SPHAdaptation>(1.15, res_factor, 0.75, level_set_refinement_ratio);
	ShellBodyFromMesh shell_body (system, "Shell_plate", mesh_2, shell_adaptation);
	ShellParticles shell_body_particles(
		shell_body,
		makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus_shell, poisson),
		makeShared<ShellParticleGeneratorLattice>(thickness), thickness
	);
	BodyRelationInner shell_body_inner(shell_body);
	ShellParticleRelaxation(in_output, shell_body, shell_body_inner, thickness, level_set_refinement_ratio);

	// Contacts
	SolidBodyRelationContact myocardium_shell_contact(myocardium_body, {&shell_body});
	SolidBodyRelationContact shell_myocardium_contact(shell_body, {&myocardium_body});
	solid_dynamics::ShellContactDensity myocardium_shell_contact_density(myocardium_shell_contact);
	solid_dynamics::ContactDensitySummation shell_myocardium_contact_density(shell_myocardium_contact);
	solid_dynamics::ContactForce myocardium_shell_contact_forces(myocardium_shell_contact);
	solid_dynamics::ContactForce shell_myocardium_contact_forces(shell_myocardium_contact);

	// SOLID BODY
	TimeStepInitialization myocardium_initialize_gravity(myocardium_body, gravity_solid);
	solid_dynamics::CorrectConfiguration corrected_configuration(myocardium_body_inner);
	solid_dynamics::KirchhoffStressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>> muscle_damping(myocardium_body_inner, 0.2, "Velocity", physical_viscosity);
	
    // // Constraint - not used
	// TriangleMeshShapeBrick holder_shape(Vec3d(10.0, 100.0, 100.0)*scale, 20, Vec3d(-10.0, 0.0, 0.0)*scale);
	// BodyRegionByParticle holder(myocardium_body, "Holder", holder_shape);
	// solid_dynamics::ConstrainSolidBodyRegion constrain_holder(myocardium_body, holder);

	// SHELL BODY
	TimeStepInitialization initialize_shell(shell_body, gravity_shell);
	thin_structure_dynamics::ShellCorrectConfiguration corrected_configuration_2(shell_body_inner);
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size_2(shell_body);
	thin_structure_dynamics::ShellStressRelaxationFirstHalf stress_relaxation_first_half_2(shell_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf stress_relaxation_second_half_2(shell_body_inner);
	DampingWithRandomChoice<DampingPairwiseInner<Vec3d>> shell_damping(shell_body_inner, 0.2, "Velocity", physical_viscosity);

	// Output
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);

	// Initialization
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	// Update the normal directions, initial conditions, and configations
	ShellDynamicsInitialCondition initial_normal_update(shell_body);
	initial_normal_update.parallel_exec();
	corrected_configuration.parallel_exec();
	corrected_configuration_2.parallel_exec();
	write_states.writeToFile(0);

	// Step parameters
	int ite = 0;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;

	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			/** Gravity. */
			myocardium_initialize_gravity.parallel_exec(dt);
			initialize_shell.parallel_exec(dt);

			myocardium_shell_contact_density.parallel_exec();
			shell_myocardium_contact_density.parallel_exec();
			myocardium_shell_contact_forces.parallel_exec();
			if (elastic_shell)
				shell_myocardium_contact_forces.parallel_exec();

			/** Stress relaxation and damping. */
			stress_relaxation_first_half.parallel_exec(dt);
			if (elastic_shell)
				stress_relaxation_first_half_2.parallel_exec(dt);

			// constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
			if (elastic_shell)
				shell_damping.parallel_exec(dt);
			// constrain_holder.parallel_exec(dt);

			stress_relaxation_second_half.parallel_exec(dt);
			if (elastic_shell)
				stress_relaxation_second_half_2.parallel_exec(dt);

			ite++;
			dt = 0.1 * computing_time_step_size_2.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			myocardium_body.updateCellLinkedList();
			shell_body.updateCellLinkedList();

			myocardium_shell_contact.updateConfiguration();
			shell_myocardium_contact.updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;


	// the shell particles can blow up sometimes, when the particle relaxation is not proper
	// max strain should be small, less than 1e-3
	EXPECT_LT(shell_body_particles.getVonMisesStrainMax(), 1e-3);
	// check that there is contact, the solid body particles don't go through the shell
	Real shell_x_min = Infinity;
	for (auto pos: shell_body_particles.pos_n_)
		if (pos[0] < shell_x_min)
			shell_x_min = pos[0];
	for (auto pos: myocardium_particles.pos_n_)
		EXPECT_LT(pos[0], shell_x_min);
}

#endif // ELASTIC_SHELL_CONTACT_H