/**
 * @file 	dw_3d_spherical_cap.cpp
 * @brief 	This is the test of using levelset to generate shell particles with single resolution and relax particles.
 * @details	We use this case to test the particle generation and relaxation by levelset for a complex thin structures geometry (3D).
 * @author 	Dong Wu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 1.0e-3;								 
Real R = 577.0 * scale;	
Real d = 10.41 * scale;
Real theta = 26.67 / 180 * Pi;
Vec3d domain_lower_bound(0.0, 508.0 * scale, 0.0);
Vec3d domain_upper_bound(523.0 * scale, 583.0 * scale, 533.0 * scale);
Real resolution_ref = 0.5 * d;
Real center = 261.3234 * scale;
Real constrained_circle = powerN(center - 1.5 * resolution_ref, 2);

// Observer location
StdVec<Vecd> observation_location = { Vecd(center, R, center) };
// level set resolution much higher than that of particles is required
Real level_set_refinement_ratio = resolution_ref / (0.1 * d);
//For material properties of the solid.
Real rho0_s = 2618.0;
Real Youngs_modulus = 7.24e10;
Real poisson = 0.3;
Real physical_viscosity = 1.0e3;

Real q = -4.14e6;
Real time_to_full_external_force = 0.0;

Real gravitational_acceleration = 0.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Define the body shape.
//----------------------------------------------------------------------
class ImportedShellModel: public ComplexShape
{
public:
	explicit ImportedShellModel(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>("./input/cap.stl", Vecd(0.0) * scale, scale);
	}
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry() {};

private:
	void tagManually(size_t index_i)
	{
		if (powerN(base_particles_.pos_[index_i][0] - center, 2) 
			+ powerN(base_particles_.pos_[index_i][2] - center, 2) > constrained_circle)
		//if (base_particles_->pos_n_[index_i][0] < center)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};
/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force
			? current_time * global_acceleration_ / time_to_full_external_force
			: global_acceleration_;
	}
};
/**
 *  The main program
 */
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	system.run_particle_relaxation_ = false;
	/** Tag for starting with relaxed body-fitted particles distribution */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition */
	system.restart_step_ = 0;
	/** Handle command line arguments. */
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody plate_body(system, makeShared<ImportedShellModel>("ImportedShellModel"));
	plate_body.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	if (!system.run_particle_relaxation_ && system.reload_particles_)
	{
		plate_body.generateParticles<ParticleGeneratorReload>(io_environment, plate_body.getName());
	}
	else
	{
		plate_body.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(io_environment);
	    //here dummy linear elastic solid is use because no solid dynamics in particle relaxation
	    plate_body.generateParticles<ThickSurfaceParticleGeneratorLattice>(d);
	}

	if (!system.run_particle_relaxation_ && !system.reload_particles_)
	{
		std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
		return 0;
	}
	plate_body.addBodyStateForRecording<Vec3d>("PriorAcceleration");
	plate_body.addBodyStateForRecording<Matd>("CorrectionMatrix");
	plate_body.addBodyStateForRecording<Vecd>("PseudoNormal");
	plate_body.addBodyStateForRecording<Vecd>("NormalDirection");

	/** Define Observer. */
	ObserverBody plate_observer(system, "PlateObserver");
	plate_observer.defineParticlesAndMaterial();
	plate_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation plate_body_inner(plate_body);
	ContactRelation plate_observer_contact(plate_observer, { &plate_body });
    //----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	if (system.run_particle_relaxation_)
	{
		SimpleDynamics <RandomizeParticlePosition>  random_imported_model_particles(plate_body);
	    // A  Physics relaxation step. 
	    relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(plate_body_inner, d, level_set_refinement_ratio);
	    relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(plate_body_inner, d);
		plate_body.addBodyStateForRecording<int>("UpdatedIndicator");
	    //----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_relaxed_particles(io_environment, system.real_bodies_);
		MeshRecordingToPlt write_mesh_cell_linked_list(io_environment, plate_body.getCellLinkedList());
		ReloadParticleIO write_particle_reload_files(io_environment, {&plate_body});
		//----------------------------------------------------------------------
	    //	Particle relaxation starts here.
	    //----------------------------------------------------------------------
	    random_imported_model_particles.parallel_exec(0.25);
	    relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	    write_relaxed_particles.writeToFile(0.0);
	    plate_body.updateCellLinkedList();
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
	    		write_relaxed_particles.writeToFile(ite_p);
	    	}
	    	relaxation_step_inner.parallel_exec();
	    	ite_p += 1;
	    }
		std::cout << "The physics relaxation process of imported model finish !" << std::endl;
		shell_normal_prediction.exec();
		write_relaxed_particles.writeToFile(ite_p);
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	SimpleDynamics<TimeStepInitialization> initialize_external_force(plate_body,
		makeShared<TimeDependentExternalForce>(Vec3d(0.0, q / (d * rho0_s) - gravitational_acceleration, 0.0)));

	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration(plate_body_inner);
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half(plate_body_inner);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half(plate_body_inner);
	BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion, BoundaryGeometry> 
		constrain_holder(boundary_geometry);
	//solid_dynamics::ConstrainSolidBodyRegionVelocity
	//	constrain_holder(plate_body, boundary_geometry, Vecd(1.0, 0.0, 1.0));
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_position_damping(0.1, plate_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_rotation_damping(0.1, plate_body_inner, "AngularVelocity", physical_viscosity);

	ObservedQuantityRecording<Vecd>
		write_plate_max_displacement("Position", io_environment, plate_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();

	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_plate_max_displacement.writeToFile(0);

	int ite = 0;
	Real end_time = 0.001;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			//plate_position_damping.parallel_exec(dt);
			//plate_rotation_damping.parallel_exec(dt);
			//constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_plate_max_displacement.writeToFile(ite);
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
