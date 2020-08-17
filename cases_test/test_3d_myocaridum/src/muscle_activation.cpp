/**
 * @file muscle_activation.cpp
 * @brief This is the first example of electro activation of myocardium
 * @author Chi Zhang and Xiangyu Hu
 * @version 0.1.0
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 1.0; 		/**< Length. */
Real PH = 1.0; 		/**< Thickness for thick plate. */
Real PW = 1.0;		/**< Width. */				
/**< Initial particle spacing. */
Real particle_spacing_ref = PH / 25.0;
Real SL = 4.0 * particle_spacing_ref; 		/**< Extension for holder. */
/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 1.0; 
Real a_0[4] = {0.059, 0.0, 0.0, 0.0};
Real b_0[4] = {8.023, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real reference_voltage = 30.0;
Real linear_active_stress_factor = - 0.5;
/** reference stress to achieve weakly compressible condition */
Real bulk_modulus = 30.0 * reference_voltage * fabs(linear_active_stress_factor);

/** Define the geometry. */
TriangleMeshShape* CreateMyocardium()
{
	Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5*PW);
	TriangleMeshShape* geometry_myocardium = new TriangleMeshShape(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
TriangleMeshShape* CreateHolder()
{
	Vecd halfsize_shape(0.5 * SL, 0.5 * PH, 0.5 * PW);
	Vecd translation_shape(-0.5 * SL, 0.5 * PH, 0.5 * PW);
	TriangleMeshShape* geometry = new TriangleMeshShape(halfsize_shape, resolution,
		translation_shape);

	return geometry;
}
/** Define the myocardium body. */
class Myocardium : public SolidBody
{
public:
	Myocardium(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		body_shape_.addTriangleMeshShape(CreateMyocardium(), ShapeBooleanOps::add);
	}
};
/**
* @brief define the Holder base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class Holder : public BodyPartByParticle
{
public:
	Holder(SolidBody *solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		body_part_shape_.addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);

		TagBodyPart();
	}
};
 /**
 * Assign case dependent muscle activation histroy  
 */
class MuscleActivation
	: public active_muscle_dynamics::ActiveMuscleSimple
{
public:
	MuscleActivation(SolidBody *myocardium)
		: active_muscle_dynamics
		::ActiveMuscleSimple(myocardium) {};
protected:
	void Update(size_t index_particle_i, Real dt) override 
	{
		ActiveMuscleParticleData& active_muscle_data_i 	= particles_->active_muscle_data_[index_particle_i];
		SolidParticleData &solid_data_i			= particles_->solid_body_data_[index_particle_i];
		BaseParticleData& base_particle_data_i	= particles_->base_particle_data_[index_particle_i];

		Real voltage = base_particle_data_i.pos_0_[0] <= 0 ? 0.0 : reference_voltage * base_particle_data_i.pos_0_[0] / PL;
		active_muscle_data_i.active_contraction_stress_ += GlobalStaticVariables::physical_time_ <= 1.0
			? linear_active_stress_factor * voltage * dt : 0.0;
	};
};
 /**
 * Setup material properties of myocardium
 */
 /**
 * Setup material properties of myocardium
 */
class MyocardiumMuscle
	: public Muscle
{
public:
	MyocardiumMuscle() : Muscle()
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
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(Vecd(-SL, -SL, -SL),
		Vecd(PL + SL, PH + SL, PW + SL), particle_spacing_ref);

	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_muscle_body =
		new Myocardium(system, "MyocardiumMuscleBody", 0, ParticlesGeneratorOps::lattice);
	ActiveMuscleParticles 	myocardium_muscle_particles(myocardium_muscle_body, new ActiveMuscle(new MyocardiumMuscle()));

	/** topology */
	SPHBodyInnerRelation* myocardium_muscle_body_inner = new SPHBodyInnerRelation(myocardium_muscle_body);

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_muscle_body_inner);
	/** Time step size calculation. */
	solid_dynamics::GetAcousticTimeStepSize 
		computing_time_step_size(myocardium_muscle_body);
	/** Compute the active contraction stress */
	MuscleActivation muscle_activation(myocardium_muscle_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(myocardium_muscle_body_inner);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(myocardium_muscle_body_inner);
	/** Constrain region of the inserted body. */
	solid_dynamics::constrainNormDirichletBoundary
		constrain_holder(myocardium_muscle_body, new Holder(myocardium_muscle_body, "Holder"), 0);
	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 1.2;
	Real output_period = end_time / 60.0;		
	Real dt = 0.0; 
	/** Statistics for computing time. */
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
				cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			muscle_activation.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
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
