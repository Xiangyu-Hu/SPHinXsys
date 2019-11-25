/**
 * @file electro_activation.cpp
 * @brief This is the first example of elctro activation of myocaridum
 * @author Chi Zhang and Xiangyu Hu
 * @version 0.1.0
 */
#include "sphinxsys.h"

using namespace SPH;

Real PL = 10.0; 		/**< Lenght. */
Real PH = 2.0; 		/**< Thickness for thick plate. */
Real PW = 2.0;		/**< Width. */				
Real SL = 2.0; 		/**< Extensioin for holder. */

Real particle_spacing_ref = PL / 50.0;		/**< Initial particle spacing. */
Real BW = particle_spacing_ref * 4; 		/**< Boundary width. */

int resolution(20);		/**< SimTK geometric modeling resolution. */
/** For material properties of the solid. */
Real rho_0 = 1.0; 
Real poisson = 0.49; 
Real Youngs_modulus = 1.0;
Real a_0[4] = {496.0, 15196.0, 3283.0, 662.0};
Real b_0[4] = {7.209, 20.417, 11.176, 9.466};
// Real a_0[4] = {496.0, 0.0, 0.0, 0.0};
// Real b_0[4] = {7.209, 0.0, 0.0, 0.0};
Vec3d d_0(1.0, 0.0, 0.0);
/** Electrophysiology prameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.04;
/** Active stress factor */
Real k_a = 250.0;
/** Define the geometry. */
Geometry *CreateMyocardium()
{
	Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_myocardium(0.5 * (PL - SL), 0.0, 0.0);
	Geometry *geometry_myocardium = new Geometry(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
Geometry *CreateHolder()
{
	Vecd halfsize_holder(0.5*SL, 0.5*PH , 0.5*PW);
	Vecd translation_holder(-0.5*SL, 0.0, 0.0);
	Geometry *geometry_holder = new Geometry(halfsize_holder,
		resolution, translation_holder);

	return geometry_holder;
}
/** Define the myocardium body. */
class Myocardium : public SolidBody
{
public:
	Myocardium(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		body_region_.add_geometry(CreateMyocardium(), RegionBooleanOps::add);
		body_region_.add_geometry(CreateHolder(),	RegionBooleanOps::add);

		body_region_.done_modeling();
	}
};
/**
* @brief define the Holder base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class Holder : public SolidBodyPart
{
public:
	Holder(SolidBody *solid_body, string constrianed_region_name)
		: SolidBodyPart(solid_body, constrianed_region_name)
	{
		body_part_region_.add_geometry(CreateHolder(), RegionBooleanOps::add);
		body_part_region_.done_modeling();

		TagBodyPartParticles();
	}
};
 /**
 * application dependent initial condition 
 */
class ElectroActivationInitialization
	: public electro_mechanics::ElectroMechanicsInitialCondition
{
public:
	ElectroActivationInitialization(SolidBody *myocardium)
		: electro_mechanics::ElectroMechanicsInitialCondition(myocardium) {};
protected:
	void Update(size_t index_particle_i, Real dt) override 
	{
		electro_mechanics::ElectroMechanicsInitialCondition::Update(index_particle_i, dt);

		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];
		if(base_particle_data_i.pos_n_[0] <= 0)
		{
			muscle_particle_data_i.voltage_n_ = 0.0;
		}else{
			muscle_particle_data_i.voltage_n_ = exp(-4.0 * (base_particle_data_i.pos_n_[0] - 5.0) 
				* (base_particle_data_i.pos_n_[0] - 5.0));
		}	
	};
};
 /**
 * Setup local properties of myocardium
 */
class MyocardiumMuscle
 	: public Muscle
{
 public:
 	MyocardiumMuscle(string myocardium_name,SPHBody *body, Real a0[4], Real b0[4], Vecd d0, 
		Real rho0, Real nu, Real eta0)
		: Muscle(myocardium_name, body, a0, b0, d0, rho0, nu, eta0){}

 	void SetupLocalProperties(SPHBody *body) override
 	{
	
 		Vecd e_1(1.0,0.0,0.0);
 		Vecd e_2(0.0,1.0,0.0);
 		for(size_t i = 0; i < body->number_of_particles_; i++)
 		{
 			f0_.push_back(e_1);
 			s0_.push_back(e_2);
 			f0f0_.push_back(SimTK::outer(e_1, e_1));
 			s0s0_.push_back(SimTK::outer(e_2, e_2));
 			f0s0_.push_back(SimTK::outer(e_1, e_2));
			Matd diff_i = d_0_[0] * Matd(1.0) + d_0_[1] * SimTK::outer(f0_[i], f0_[i]);
			diff_cd_0.push_back(inverseCholeskyDecomposition(diff_i));
		}
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(Vecd(-SL - BW, -0.5 * (PL + BW), -0.5 * (SL + BW)),
		Vecd(PL, 0.5 * (PL + BW), 0.5 * (SL + BW)), particle_spacing_ref);
	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body =
		new Myocardium(system, "MyocardiumBody", 0, ParticlesGeneratorOps::lattice);
	MyocardiumMuscle 	material("Muscle", myocardium_body, a_0, b_0,d_0, rho_0, poisson, 1.0);
	MuscleParticles 	particles(myocardium_body);
	AlievPanfilowModel 	reaction_model("TwoVariableModel", myocardium_body, c_m, k, a, mu_1, mu_2, epsilon, k_a);
	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bidies to build neighbor particle lists
	 */
	SPHBodyTopology body_topology = { { myocardium_body,{} } };
	system.SetBodyTopology(&body_topology);
	/** Setting up the simulation. */
	system.SetupSPHSimulation();
	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** initial condition for the elastic solid bodies */
	ElectroActivationInitialization 
		setup_initial_condition(myocardium_body);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_body, {});
	/** Time step size caclutation. */
	solid_dynamics::GetAcousticTimeStepSize 
		computing_time_step_size(myocardium_body);
	/** Compute the active contraction stress */
	electro_mechanics::computeActiveContractionStress
		compute_active_contraction_stress(myocardium_body);
	/** active-pative stress relaxation. */
	electro_mechanics::ActivePassiveStressRelaxationFirstStep
		stress_relaxation_first_step(myocardium_body);
	electro_mechanics::ActivePassiveStressRelaxationSecondStep
		stress_relaxation_second_step(myocardium_body);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** apply initial condition */
	setup_initial_condition.exec();
	material.SetupLocalProperties(myocardium_body);
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 10.0;
	Real output_period = end_time / 100.0;		
	Real dt = 0.0; 
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < output_period) 
		{
			if (ite % 100 == 0) {
				cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			compute_active_contraction_stress.parallel_exec(dt);
			stress_relaxation_first_step.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_step.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
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

	return 0;
}
