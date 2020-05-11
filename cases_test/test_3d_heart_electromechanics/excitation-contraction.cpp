/**
 * @file 	biventricular_electromechanics.cpp
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
 */
/** 
 * SPHinXsys Library. 
 */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/heart-new.stl";
Real length_scale = 1.0;
Real time_scale  = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
/** Paremeters and physical properties. */
Vec3d domain_lower_bound(-55.0 * length_scale,-75.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);	
Real dp_0 	= (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0;	/**< Initial particle spacing. */
/** Material properties. */
Real rho_0 = 1.06e-3; 	
/** Active stress factor */
Real k_a = 100 * stress_scale;
Real a_0[4] = {496.0 * stress_scale, 15196.0 * stress_scale, 3283.0 * stress_scale, 662.0 * stress_scale};
Real b_0[4] = {7.209, 20.417, 11.176, 9.466};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a_0[0] *(1.0 + poisson) / (3.0 * (1.0 - 2.0*poisson));
/** Electrophysiology prameters. */
Real diffusion_coff_ = 0.8;
Real bias_diffusion_coff_ = 0.0;
/** Electrophysiology prameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
/** Define the geometry. */
Geometry *CreateHeart()
{
	Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
	Geometry *geometry_myocardium = new Geometry(full_path_to_stl_file, translation, length_scale);

	return geometry_myocardium;
}
Geometry *CreateBaseShape()
{
	Real l = domain_upper_bound[0] - domain_lower_bound[0];
	Real w = domain_upper_bound[2] - domain_lower_bound[2];
	Vecd halfsize_shape(0.5 * l, 1.0 * dp_0, 0.5 * w);
	Vecd translation_shape(-10.0 * length_scale, -1.0 * dp_0, 0.0);
	Geometry *geometry = new Geometry(halfsize_shape, 20, translation_shape);

	return geometry;
}
/**
 * Setup eletro_physiology reation properties
 */
class MuscleReactionModel : public AlievPanfilowModel
{
public:
	MuscleReactionModel() : AlievPanfilowModel()
	{
		/** Basic reaction parameters*/
		k_a_ 		= k_a;
		c_m_ 		= c_m;
		k_ 			= k;
		a_ 			= a;
		b_ 			= b;
		mu_1_ 		= mu_1;
		mu_2_ 		= mu_2;
		epsilon_ 	= epsilon;
		/** Compute the derived material parameters*/
		assignDerivedReactionParameters();
	}
};
/**
 * Setup material properties of myocardium
 */
class MyocardiumPhysiology
 	: public LocalMonoFieldElectroPhysiology
{
 public:
 	MyocardiumPhysiology(ElectroPhysiologyReaction* electro_physiology_reaction)
		: LocalMonoFieldElectroPhysiology(electro_physiology_reaction)
	{
		/** Basic material parameters*/
		diff_cf_ 		= diffusion_coff_;
		bias_diff_cf_ 	= bias_diffusion_coff_;
		bias_direction_ = fiber_direction;
		/** Compute the derived material parameters. */
		assignDerivedMaterialParameters();
		/** Create the vector of diffusion properties. */
		initializeDiffusion();
	}
};
 /**
 * Setup material properties of myocardium
 */
class MyocardiumMuscle
 	: public LocallyOrthotropicMuscle
{
 public:
 	MyocardiumMuscle() : LocallyOrthotropicMuscle()
	{
		rho_0_ = rho_0;
		bulk_modulus_ = bulk_modulus;
		f0_ = fiber_direction;
		s0_ = sheet_direction;
		std::copy(a_0, a_0 + 4, a_0_);
		std::copy(b_0, b_0 + 4, b_0_);

		assignDerivedMaterialParameters();
	}
};
/** 
 * Define geometry and initial conditions of SPH bodies. 
 */
class HeartBody : public SolidBody
{
public:
	HeartBody(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		body_region_.add_geometry(CreateHeart(), RegionBooleanOps::add);
		body_region_.done_modeling();
	}
};
/**
* @brief define the beam base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class MuscleBase : public BodyPartByParticle
{
public:
	 MuscleBase(SolidBody *solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		body_part_region_.add_geometry(CreateBaseShape(), RegionBooleanOps::add);
		body_part_region_.done_modeling();
		/** Tag the constrained particles to the base for constraint. */
		TagBodyPartParticles();
	}
};
 /**
 * application dependent initial condition 
 */
class ApplyStimulusCurrentSI
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_particle_i, Real dt) override
	{
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		DiffusionReactionData& electro_physiology_data_i = particles_->diffusion_reaction_data_[index_particle_i];

		if( -30.0  * length_scale <= base_particle_data_i.pos_n_[0] && base_particle_data_i.pos_n_[0] <= -15.0  * length_scale)
		{
			if( -2.0  * length_scale <= base_particle_data_i.pos_n_[1] && base_particle_data_i.pos_n_[1] <= 0.0)
			{
				if( -3.0  * length_scale <= base_particle_data_i.pos_n_[2] && base_particle_data_i.pos_n_[2] <= 3.0  * length_scale)
				{
					electro_physiology_data_i.species_n_[voltage_] = 0.92;
				}
			}
		}
	};
public:
	ApplyStimulusCurrentSI(SolidBody* muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle) 
	{
		voltage_ = material_->getSpeciesIndexMap()["Voltage"];
	};
};
 /**
 * application dependent initial condition 
 */
class ApplyStimulusCurrentSII
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_particle_i, Real dt) override
	{
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		DiffusionReactionData& electro_physiology_data_i = particles_->diffusion_reaction_data_[index_particle_i];

		if( 0.0 <= base_particle_data_i.pos_n_[0] && base_particle_data_i.pos_n_[0] <= 6.0 * length_scale) 
		{
			if( -6.0 * length_scale <= base_particle_data_i.pos_n_[1])
			{
				if( 12.0 * length_scale <= base_particle_data_i.pos_n_[2])
				{
					electro_physiology_data_i.species_n_[voltage_] = 0.95;
				}
			}
		}
	};
public:
	ApplyStimulusCurrentSII(SolidBody* muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle) 
	{
		voltage_ = material_->getSpeciesIndexMap()["Voltage"];
	};
};
/**
 * Voltage observer body definition.
 */
class VoltageObserver : public FictitiousBody
{
public:
	VoltageObserver(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(make_pair(Point(-45.0 * length_scale, -30.0 * length_scale, 0.0),  0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -30.0 * length_scale, 26.0 * length_scale), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(-30.0 * length_scale, -50.0 * length_scale, 0.0),  0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -50.0 * length_scale, 20.0 * length_scale), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -70.0 * length_scale, 0.0),  0.0));
	}
};
/**
 * Muscle observer.
 */
class MyocardiumObserver : public FictitiousBody
{
public:
	MyocardiumObserver(SPHSystem &system, string body_name, 
			   int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(make_pair(Point(-45.0 * length_scale, -30.0 * length_scale, 0.0),  0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -30.0 * length_scale, 26.0 * length_scale), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(-30.0 * length_scale, -50.0 * length_scale, 0.0),  0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -50.0 * length_scale, 20.0 * length_scale), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(0.0,   -70.0 * length_scale, 0.0),  0.0));
	}
};
/** 
 * The main program. 
 */
int main()
{
	/** 
	 * Build up context -- a SPHSystem. 
	 */
	SPHSystem system(domain_lower_bound, domain_upper_bound, dp_0, 6);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	/** Tag for reload initially repaxed particles. */
	system.reload_particles_ = true;
	/** Configuration of materials, crate particle container and electrophysiology body. */
	HeartBody *physiology_body  					
		=  new HeartBody(system, "ExcitationHeart", 0, ParticlesGeneratorOps::lattice);
	MuscleReactionModel *muscle_reaction_model 	= new MuscleReactionModel();
	MyocardiumPhysiology 	*myocardium_excitation 	= new MyocardiumPhysiology(muscle_reaction_model);
	ElectroPhysiologyParticles 	physiology_articles(physiology_body, myocardium_excitation);

	/** Configuration of materials, crate particle container and electrophysiology body. */
	HeartBody *mechanics_body  					
		=  new HeartBody(system, "ContractionHeart", 0, ParticlesGeneratorOps::lattice);
	MyocardiumMuscle *myocardium_muscle = new MyocardiumMuscle();
	ActiveMuscle 	*myocardium_contraction = new ActiveMuscle(myocardium_muscle);
	ActiveMuscleParticles 	mechanics_particles(mechanics_body, myocardium_contraction);

	/**
	 * Particle and body creation of fluid observer.
	 */
	VoltageObserver *voltage_observer 
		= new VoltageObserver(system, "VoltageObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles 		observer_particles(voltage_observer);
	/** Define muscle Observer. */
	MyocardiumObserver *myocardium_observer 
		= new MyocardiumObserver(system, "MyocardiumObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles 	disp_observer_particles(myocardium_observer);

	/** Outputs. */
	In_Output 					in_output(system);
	WriteBodyStatesToVtu 		write_states(in_output, system.real_bodies_);
	ReadReloadParticle			excitation_reload_particles(in_output, {physiology_body}, { "ExcitationHeart" });;
	ReadReloadParticle			contraction_reload_particles(in_output, {mechanics_body}, { "ContractionHeart" });;
	/** Read material property, e.g., sheet and fiber, from xml file. */
	ReadReloadMaterialProperty  read_material_property(in_output, myocardium_muscle);
	/** Set body contact map. */
	SPHBodyTopology body_topology = {{physiology_body, {mechanics_body}}, {mechanics_body, {physiology_body}},
									 {voltage_observer, {physiology_body}}, {myocardium_observer, {mechanics_body}} };
	system.SetBodyTopology(&body_topology);
	/**
	 * Simulation data structure set up.
	 */
	/** check whether reload particles. */
	if (system.reload_particles_)
	{
		excitation_reload_particles.ReadFromFile();
		contraction_reload_particles.ReadFromFile();
		read_material_property.ReadFromFile();
		myocardium_excitation->assignFiberProperties(myocardium_muscle->local_f0_);
	}
	system.SetupSPHSimulation();
	/** 
	 * Corrected strong configuration. 
	 */	
	solid_dynamics::CorrectConfiguration 					
		correct_configuration_excitation(physiology_body);
	/** 
	 * Time step size caclutation. 
	 */
	electro_physiology::GetElectroPhysiologyTimeStepSize 
		get_physiology_time_step(physiology_body);
	/** 
	 * Diffusion process for diffusion body. 
	 */
	electro_physiology::ElectroPhysiologyDiffusionRelaxation 	
		diffusion_relaxation(physiology_body);
	/** 
	 * Sovlers for ODE system 
	 */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward 		
		reaction_relaxation_forward(physiology_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward 	
		reaction_relaxation_backward(physiology_body);
	/**
	 * IO for observer.
	 */
	WriteObservedDiffusionReactionQuantity<ElectroPhysiologyParticles>
		write_voltage("Voltage", in_output, voltage_observer, physiology_body);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_displacement("Displacement", in_output, myocardium_observer, mechanics_body);
	/**
	 * Apply the Iron stimulus.
	 */
	ApplyStimulusCurrentSI		
		apply_stimulus_s1(physiology_body);
	ApplyStimulusCurrentSII		
		apply_stimulus_s2(physiology_body);
	/**
	 * Active mechanics. */
	solid_dynamics::CorrectConfiguration 
		correct_configuration_contraction(mechanics_body);
	/** */
	observer_dynamics::CorrectKenelWeightsForInterpolation
		correct_kernel_weights_for_interpolation(mechanics_body, {physiology_body});
	/** Interpolate the active contract stress from eletrophyisology body. */
	observer_dynamics::InterpolatingADiffusionReactionQuantity<ActiveMuscleParticles, ActiveMuscleData,
		ElectroPhysiologyParticles, &ActiveMuscleParticles::active_muscle_data_, &ActiveMuscleData::active_contraction_stress_>
		active_stress_interpolation("ActiveContractionStress", mechanics_body, { physiology_body });
	/** Interpolate the particle position in physiology_body  from mechanics_body. */
	observer_dynamics::InterpolatingAQuantity<Vecd, BaseParticles, BaseParticleData, BaseParticles, BaseParticleData, 
		&BaseParticles::base_particle_data_, &BaseParticles::base_particle_data_,
		&BaseParticleData::pos_n_, &BaseParticleData::pos_n_>
		interpolation_particle_position(physiology_body, { mechanics_body });
	/** Time step size caclutation. */
	solid_dynamics::GetAcousticTimeStepSize 
		get_mechanics_time_step(mechanics_body);
	/** active-pative stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(mechanics_body);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(mechanics_body);
	/** Constrain region of the inserted body. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(mechanics_body, new MuscleBase(mechanics_body, "Holder"));
	/** 
	 * Pre-simultion. 
	 */
	correct_configuration_excitation.parallel_exec();
	correct_configuration_contraction.parallel_exec();
	correct_kernel_weights_for_interpolation.parallel_exec();
	/** 
	 * Output global basic parameters. 
	 */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_voltage.WriteToFile(GlobalStaticVariables::physical_time_);
	write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * Physical parameters for main loop. 
	 */
	int screen_output_interval 	= 10;
	int ite 					= 0;
	int reaction_step 			= 2;
	Real End_Time 				= 100;
	Real Ouput_T 				= End_Time / 200.0;
	Real Observer_time 			= 0.01 * Ouput_T;	
	Real dt 					= 0.0; 				/**< Default acoustic time step sizes for physiology. */
	Real dt_s 					= 0.0;				/**< Default acoustic time step sizes for mechanics. */
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	cout << "Main Loop Starts Here : " << "\n";
	/** Main loop starts here. */ 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < Ouput_T) 
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observer_time) 
			{
				if (ite % screen_output_interval == 0) 
				{
					cout << fixed << setprecision(9) << "N=" << ite << "	Time = "
						<< GlobalStaticVariables::physical_time_
						<< "	dt = " << dt 
						<< "	dt_s = " << dt_s << "\n";
				}
				/** Apply stimulus excitation. */
				if( 0 <= GlobalStaticVariables::physical_time_ 
					&&  GlobalStaticVariables::physical_time_ <= 0.5)
				{
					apply_stimulus_s1.parallel_exec(dt);
				}
				/** Single spiral wave. */
				// if( 105 <= GlobalStaticVariables::physical_time_ 
				// 	&&  GlobalStaticVariables::physical_time_ <= 105.2)	
				// {
				// 	apply_stimulus_s2.parallel_exec(dt);
				// }
				/**Strang's splitting method. */
				//forward reaction 
				int ite_forward = 0;
				while (ite_forward < reaction_step )
				{
					reaction_relaxation_forward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_forward ++;
				}
				/** 2nd Runge-Kutta scheme for diffusion. */
				diffusion_relaxation.parallel_exec(dt);

				//backward reaction
				int ite_backward = 0;
				while (ite_backward < reaction_step)
				{
					reaction_relaxation_backward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_backward ++;
				}

				active_stress_interpolation.parallel_exec();

				Real dt_s_sum = 0.0;
				while (dt_s_sum < dt) 
				{
					dt_s = get_mechanics_time_step.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_holder.parallel_exec(dt_s);
					stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}

				ite++;
				dt = get_physiology_time_step.parallel_exec();

				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_voltage.WriteToFile(GlobalStaticVariables::physical_time_);
			write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}
		tick_count t2 = tick_count::now();
		interpolation_particle_position.parallel_exec();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}