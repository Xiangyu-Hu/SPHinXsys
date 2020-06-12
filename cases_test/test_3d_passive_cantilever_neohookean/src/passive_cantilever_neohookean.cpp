/**
 * @file pssive_cantilever_neohookean.cpp
 * @brief This is the example of myocaridum with simple neohookean tissue model 
 * @author Bence Rochlitz, Chi Zhang  and Xiangyu Hu
 * @version 0.1.0
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 0.1;
Real PH = 0.04;
Real PW = 0.04;
Real SL = 0.02;
Real particle_spacing_ref = PH / 6.0;		/**< Initial particle spacing. */
Real BW = particle_spacing_ref * 4; 		/**< Boundary width. */
/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 1265.0; // Gheorghe 2019 
Real poisson = 0.45; // nearly incompressible
Real Youngs_modulus = 5e4; // Sommer 2015
Real physical_viscosity = 0.0; //physical damping
Real gravity_g = 9.8; 					/**< Value of gravity. */
Real time_to_full_gravity = 4.0;

Vec3d d_0(1.0e-4, 0.0, 0.0);
/** Define the geometry. */
Geometry *CreateMyocardium()
{
	Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
	Geometry *geometry_myocardium = new Geometry(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
Geometry *CreateHolder()
{
	Vecd halfsize_holder(0.5*SL, 0.5*PH , 0.5*PW);
	Vecd translation_holder(-0.5*SL, 0.5 * PH, 0.5 * PW);
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
class Holder : public BodyPartByParticle
{
public:
	Holder(SolidBody *solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		body_part_region_.add_geometry(CreateHolder(), RegionBooleanOps::add);
		body_part_region_.done_modeling();

		TagBodyPartParticles();
	}
};
/**
 * Setup material properties of myocardium
 */
class MyocardiumMuscle : public NeoHookeanSolid
{
public:
	MyocardiumMuscle() : NeoHookeanSolid()
	{
		rho_0_ 	= rho_0;
		E_0_ = Youngs_modulus;
		nu_ = poisson;
		eta_0_ = physical_viscosity;

		assignDerivedMaterialParameters();
	}
};
/**
 * define time dependent gravity
 */
class TimeDependentGravity : public Gravity
{
public:
	TimeDependentGravity(Vecd gravity_vector) 
		: Gravity(gravity_vector) {}
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_gravity ?
			current_time * global_acceleration_ / time_to_full_gravity : global_acceleration_;
	}
};
 //define an observer body
class MyocardiumObserver : public FictitiousBody
{
public:
	MyocardiumObserver(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		body_input_points_volumes_.push_back(make_pair(Point(PL, PH, PW), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(Vecd(-SL - BW, BW, BW),
		Vecd(PL + BW, PH + BW, PH + BW), particle_spacing_ref, 6);

	/** Define the external force. */
	TimeDependentGravity gravity(Vec3d(0.0, -gravity_g, 0.0));

	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body =
		new Myocardium(system, "MyocardiumBody", 0, ParticlesGeneratorOps::lattice);
	MyocardiumMuscle 	*muscle_material = new MyocardiumMuscle();
	ElasticSolidParticles 	particles(myocardium_body, muscle_material);
	/** Define Observer. */
	MyocardiumObserver *myocardium_observer 
		= new MyocardiumObserver(system, "MyocardiumObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles observer_particles(myocardium_observer);
	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bidies to build neighbor particle lists
	 */
	SPHBodyTopology body_topology = { { myocardium_body,{} }, {myocardium_observer, {myocardium_body}} };
	system.SetBodyTopology(&body_topology);

	//-------- common paritcle dynamics ----------------------------------------
	InitializeATimeStep 	initialize_gravity(myocardium_body, &gravity);

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_body);
	/** Time step size caclutation. */
	solid_dynamics::GetAcousticTimeStepSize 
		computing_time_step_size(myocardium_body);
	/** active-pative stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(myocardium_body);
	/** Setup the damping stress, if you know what you are doing. */
	//stress_relaxation_first_step.setupDampingStressFactor(1.0);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(myocardium_body);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));

	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_displacement("Displacement", in_output, myocardium_observer, myocardium_body);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 8.0;
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

			initialize_gravity.parallel_exec(); // gravity force
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
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
