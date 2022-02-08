/**
 * @file 	pkj_lv_electrocontraction.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			Chi Zhang
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 *@version 0.3
 *			Here, the coupling with Purkinje network will be condcuted.
 */
/**  SPHinXsys Library. */
#include "sphinxsys.h"
#include "pkj_lv_electrocontraction.h"
/** Namespace cite here. */
using namespace SPH;
/** 
 * The main program. 
 */
int main(int ac, char *av[])
{
	/** 
	 * Build up context -- a SPHSystem. 
	 */
	SPHSystem system(system_domain_bounds, dp_0);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** Tag for reload initially repaxed particles. */
	system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** in- and output environment. */
	In_Output in_output(system);

	/** Creat a heart body for physiology, material and particles */
	HeartBody physiology_body(system, "ExcitationHeart");
	SharedPtr<ParticleGenerator> physiology_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		physiology_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, physiology_body.getBodyName());
	AlievPanfilowModel muscle_reaction_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	SharedPtr<LocalMonoFieldElectroPhysiology> myocardium_excitation =
		makeShared<LocalMonoFieldElectroPhysiology>(muscle_reaction_model, diffusion_coff, bias_coff, fiber_direction);
	ElectroPhysiologyParticles physiology_articles(physiology_body, myocardium_excitation, physiology_body_particle_generator);

	/** Creat a heart body for excitation-contraction, material and particles */
	HeartBody mechanics_body(system, "ContractionHeart");
	SharedPtr<ParticleGenerator> mechanics_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		mechanics_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, mechanics_body.getBodyName());
	SharedPtr<ActiveMuscle<LocallyOrthotropicMuscle>> myocardium_muscle =
		makeShared<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
	ActiveMuscleParticles mechanics_particles(mechanics_body, myocardium_muscle, mechanics_body_particle_generator);

	/** check whether reload material properties. */
	if (!system.run_particle_relaxation_ && system.reload_particles_)
	{
		ReloadMaterialParameterIO read_muscle_fiber_and_sheet(in_output, myocardium_muscle);
		ReloadMaterialParameterIO read_myocardium_excitation_fiber(in_output, myocardium_excitation, myocardium_muscle->LocalParametersName());
		read_muscle_fiber_and_sheet.readFromFile();
		read_myocardium_excitation_fiber.readFromFile();
	}

	/** Creat a Purkinje network for fast diffusion, material and particles */
	PurkinjeBody pkj_body(system, "Purkinje");
	AlievPanfilowModel pkj_reaction_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	SharedPtr<MonoFieldElectroPhysiology> pkj_myocardium_muscle =
		makeShared<MonoFieldElectroPhysiology>(pkj_reaction_model, diffusion_coff * acceleration_factor, bias_coff, fiber_direction);
	ElectroPhysiologyReducedParticles pkj_muscle_particles(pkj_body, pkj_myocardium_muscle,
		makeShared<NetworkGeneratorWithExtraCheck>(starting_point, second_point, 50, 1.0));
	TreeTerminates pkj_leaves(pkj_body);
	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		HeartBody relax_body(system, "RelaxationHeart");
		StdVec<std::string> species_name_list{"Phi"};
		SharedPtr<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>> relax_body_material =
			makeShared<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>>(
				species_name_list, rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
		relax_body_material->initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
		DiffusionReactionParticles<ElasticSolidParticles, LocallyOrthotropicMuscle> diffusion_particles(relax_body, relax_body_material);
		/** topology */
		BodyRelationInner relax_body_inner(relax_body);
		/** Random reset the relax solid particle position. */
		RandomizePartilePosition random_particles(relax_body);
		/**Algorithms for particle relaxation.*/
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(relax_body_inner);
		/** Diffusion process.*/
		/** Time step for diffusion. */
		GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(relax_body);
		/** Diffusion process for diffusion body. */
		DiffusionRelaxation diffusion_relaxation(relax_body_inner);
		/** Compute the fiber and sheet after diffusion. */
		ComputeFiberandSheetDirections compute_fiber_sheet(relax_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_relax_body_state(in_output, {relax_body});
		/** Write the particle reload files. */
		ReloadParticleIO
			write_particle_reload_files(in_output, {&relax_body, &relax_body},
										{physiology_body.getBodyName(), mechanics_body.getBodyName()});
		/** Write material property to xml file. */
		ReloadMaterialParameterIO
			write_material_property(in_output, relax_body_material, myocardium_muscle->LocalParametersName());
		/**Physics relaxation starts here.*/
		/** Relax the elastic structure. */
		random_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_relax_body_state.writeToFile(0.0);
		/**
		 * From here the time stepping begines.
		 * Set the starting time.
		 */
		int ite = 0;
		int relax_step = 1000;
		int diffusion_step = 100;
		while (ite < relax_step)
		{
			relaxation_step_inner.parallel_exec();
			ite++;
			if (ite % 100 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relax_body_state.writeToFile(Real(ite) * 1.0e-4);
			}
		}
		/** constraint boundary condition for diffusion. */
		BodySurface surface_part(relax_body);
		DiffusionBCs impose_diffusion_bc(relax_body, surface_part);
		impose_diffusion_bc.parallel_exec();
		write_relax_body_state.writeToFile(Real(ite) * 1.0e-4);
		Real dt = get_time_step_size.parallel_exec();
		while (ite <= diffusion_step + relax_step)
		{
			diffusion_relaxation.parallel_exec(dt);
			impose_diffusion_bc.parallel_exec();
			if (ite % 10 == 0)
			{
				cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
				write_relax_body_state.writeToFile(Real(ite) * 1.0e-4);
			}
			ite++;
		}
		compute_fiber_sheet.exec();
		ite++;
		write_relax_body_state.writeToFile(Real(ite) * 1.0e-4);
		compute_fiber_sheet.parallel_exec();
		write_material_property.writeToFile(0);
		write_particle_reload_files.writeToFile(0);

		return 0;
	}
	//----------------------------------------------------------------------
	//	SPH Observation section
	//----------------------------------------------------------------------
	ObserverBody voltage_observer(system, "VoltageObserver");
	ObserverParticles observer_particles(voltage_observer, makeShared<ObserverParticleGenerator>());
	ObserverBody myocardium_observer(system, "MyocardiumObserver");
	ObserverParticles disp_observer_particles(myocardium_observer, makeShared<ObserverParticleGenerator>());

	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	/** topology */
	BodyRelationInner physiology_body_inner(physiology_body);
	BodyRelationInner mechanics_body_inner(mechanics_body);
	BodyRelationContact physiology_body_contact(physiology_body, {&mechanics_body});
	BodyRelationContact mechanics_body_contact(mechanics_body, {&physiology_body});
	BodyRelationContact voltage_observer_contact(voltage_observer, {&physiology_body});
	BodyRelationContact myocardium_observer_contact(myocardium_observer, {&mechanics_body});
	ComplexBodyRelation physiology_body_complex(physiology_body, {&pkj_leaves});
	GenerativeBodyRelationInner pkj_inner(pkj_body);

	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration correct_configuration_excitation(physiology_body_inner);
	/** Time step size calculation. */
	electro_physiology::GetElectroPhysiologyTimeStepSize get_myocardium_physiology_time_step(physiology_body);
	/** Diffusion process for diffusion body. */
	electro_physiology::ElectroPhysiologyDiffusionRelaxationComplex myocardium_diffusion_relaxation(physiology_body_complex);
	/** Solvers for ODE system */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward myocardium_reaction_relaxation_forward(physiology_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward myocardium_reaction_relaxation_backward(physiology_body);
	/** Physiology for PKJ*/
	/** Time step size calculation. */
	electro_physiology::GetElectroPhysiologyTimeStepSize get_pkj_physiology_time_step(pkj_body);
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner pkj_diffusion_relaxation(pkj_inner);
	/** Solvers for ODE system */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward pkj_reaction_relaxation_forward(pkj_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward pkj_reaction_relaxation_backward(pkj_body);
	/**IO for observer.*/
	ObservedQuantityRecording<Real> write_voltage("Voltage", in_output, voltage_observer_contact);
	ObservedQuantityRecording<Vecd> write_displacement("Position", in_output, myocardium_observer_contact);
	/**Apply the Iron stimulus.*/
	ApplyStimulusCurrentToMmyocardium apply_stimulus_myocardium(physiology_body);
	ApplyStimulusCurrentToPKJ apply_stimulus_pkj(pkj_body);
	/** Active mechanics. */
	solid_dynamics::CorrectConfiguration correct_configuration_contraction(mechanics_body_inner);
	/** Observer Dynamics */
	observer_dynamics::CorrectInterpolationKernelWeights
		correct_kernel_weights_for_interpolation(mechanics_body_contact);
	/** Interpolate the active contract stress from electrophysiology body. */
	observer_dynamics::InterpolatingAQuantity<Real>
		active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
	/** Interpolate the particle position in physiology_body  from mechanics_body. */
	observer_dynamics::InterpolatingAQuantity<Vecd>
		interpolation_particle_position(physiology_body_contact, "Position", "Position");
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize get_mechanics_time_step(mechanics_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(mechanics_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(mechanics_body_inner);
	/** Constrain region of the inserted body. */
	MuscleBaseShapeParameters muscle_base_parameters;
	TriangleMeshShapeBrick muscle_base_shape(muscle_base_parameters);
	BodyRegionByParticle muscle_base(mechanics_body, "Holder",  muscle_base_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(mechanics_body, muscle_base);
	/** 
	 * Pre-simultion. 
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	correct_configuration_excitation.parallel_exec();
	correct_configuration_contraction.parallel_exec();
	correct_kernel_weights_for_interpolation.parallel_exec();
	/**Output global basic parameters. */
	write_states.writeToFile(0);
	write_voltage.writeToFile(0);
	write_displacement.writeToFile(0);
	write_states.writeToFile(0);
	/**
	 * main loop. 
	 */
	int screen_output_interval = 10;
	int ite = 0;
	int reaction_step = 2;
	Real End_Time = 80;
	Real Ouput_T = End_Time / 200.0;
	Real Observer_time = 0.01 * Ouput_T;
	Real dt_myocardium = 0.0;
	Real dt_pkj = 0.0;
	Real dt_muscle = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	cout << "Main Loop Starts Here : "
		 << "\n";
	/** Main loop starts here. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < Ouput_T)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observer_time)
			{
				if (ite % screen_output_interval == 0)
				{
					cout << fixed << setprecision(9) << "N=" << ite << "	Time = "
						 << GlobalStaticVariables::physical_time_
						 << "	dt_pkj = " << dt_pkj
						 << "	dt_myocardium = " << dt_myocardium
						 << "	dt_muscle = " << dt_muscle << "\n";
				}
				/** Apply stimulus excitation. */
				// if( 0 <= GlobalStaticVariables::physical_time_
				// 	&&  GlobalStaticVariables::physical_time_ <= 0.5)
				// {
				// 	apply_stimulus_myocardium.parallel_exec(dt_myocardium);
				// }

				Real dt_pkj_sum = 0.0;
				while (dt_pkj_sum < dt_myocardium)
				{
					/**
					 * When network generates particles, the final particle spacing, which is after particle projected in to 
					 * complex geometry, may small than the reference one, therefore, a smaller time step size is required. 
					 */
					dt_pkj = 0.5 * get_pkj_physiology_time_step.parallel_exec();
					if (dt_myocardium - dt_pkj_sum < dt_pkj)
						dt_pkj = dt_myocardium - dt_pkj_sum;

					if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= 0.5)
					{
						apply_stimulus_pkj.parallel_exec(dt_pkj);
					}
					/**Strang splitting method. */
					int ite_pkj_forward = 0;
					while (ite_pkj_forward < reaction_step)
					{
						pkj_reaction_relaxation_forward.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						ite_pkj_forward++;
					}
					/** 2nd Runge-Kutta scheme for diffusion. */
					pkj_diffusion_relaxation.parallel_exec(dt_pkj);
					//backward reaction
					int ite_pkj_backward = 0;
					while (ite_pkj_backward < reaction_step)
					{
						pkj_reaction_relaxation_backward.parallel_exec(0.5 * dt_pkj / Real(reaction_step));
						ite_pkj_backward++;
					}

					dt_pkj_sum += dt_pkj;
				}

				/**Strang splitting method. */
				int ite_forward = 0;
				while (ite_forward < reaction_step)
				{
					myocardium_reaction_relaxation_forward.parallel_exec(0.5 * dt_myocardium / Real(reaction_step));
					ite_forward++;
				}
				/** 2nd Runge-Kutta scheme for diffusion. */
				myocardium_diffusion_relaxation.parallel_exec(dt_myocardium);

				//backward reaction
				int ite_backward = 0;
				while (ite_backward < reaction_step)
				{
					myocardium_reaction_relaxation_backward.parallel_exec(0.5 * dt_myocardium / Real(reaction_step));
					ite_backward++;
				}

				active_stress_interpolation.parallel_exec();
				Real dt_muscle_sum = 0.0;
				while (dt_muscle_sum < dt_myocardium)
				{
					dt_muscle = get_mechanics_time_step.parallel_exec();
					if (dt_myocardium - dt_muscle_sum < dt_muscle)
						dt_muscle = dt_myocardium - dt_muscle_sum;
					stress_relaxation_first_half.parallel_exec(dt_muscle);
					constrain_holder.parallel_exec(dt_muscle);
					stress_relaxation_second_half.parallel_exec(dt_muscle);
					dt_muscle_sum += dt_muscle;
				}

				ite++;
				dt_myocardium = get_myocardium_physiology_time_step.parallel_exec();

				relaxation_time += dt_myocardium;
				integration_time += dt_myocardium;
				GlobalStaticVariables::physical_time_ += dt_myocardium;
			}
			write_voltage.writeToFile(ite);
			write_displacement.writeToFile(ite);
		}
		tick_count t2 = tick_count::now();
		interpolation_particle_position.parallel_exec();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}