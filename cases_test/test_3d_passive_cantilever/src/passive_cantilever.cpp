/**
 * @file passive_cantilever.cpp
 * @brief This is the first example of myocardium 
 * @author Chi Zhang and Xiangyu Hu
 * @version 0.1.0
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 6.0; 
Real PH = 1.0;
Real PW = 1.0;		
Real SL = 0.5; 
Real particle_spacing_ref = PH / 12.0;		/**< Initial particle spacing. */
Real BW = particle_spacing_ref * 4; 		/**< Boundary width. */
/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 1100.0; 
Real poisson = 0.45; 
Real Youngs_modulus = 1.7e7;
Real a = Youngs_modulus/(2.0 *(1.0 + poisson));
Real a_f = 0.0 * a;
Real a_0[4] = {a, a_f, 0.0, 0.0};
Real b_0[4] = {1.0, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real bulk_modulus = Youngs_modulus / 3.0 / (1.0 - 2.0 * poisson);

Vec3d d_0(1.0e-4, 0.0, 0.0);
/** Define the geometry. */
TriangleMeshShape* CreateMyocardium()
{
	Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
	TriangleMeshShape* geometry_myocardium = new TriangleMeshShape(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
TriangleMeshShape* CreateHolder()
{
	Vecd halfsize_holder(0.5*SL, 0.5*PH , 0.5*PW);
	Vecd translation_holder(-0.5*SL, 0.5 * PH, 0.5 * PW);
	TriangleMeshShape* geometry_holder = new TriangleMeshShape(halfsize_holder,
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
		body_shape_.addTriangleMeshShape(CreateMyocardium(), ShapeBooleanOps::add);
		body_shape_.addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);
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
	Holder(SolidBody *solid_body, string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		body_part_shape_.addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);

		TagBodyPart();
	}
};
/**
 * Setup material properties of myocardium
 */
class MyocardiumMuscle : public Muscle
{
public:
	MyocardiumMuscle() : Muscle() 
	{
		rho_0_ 	= rho_0;
		bulk_modulus_ = bulk_modulus;
		f0_ 	= fiber_direction;
		s0_ 	= sheet_direction;
		std::copy(a_0, a_0 + 4, a_0_);
		std::copy(b_0, b_0 + 4, b_0_);

		assignDerivedMaterialParameters();
	}
};
/**
 * application dependent initial condition 
 */
class MyocardiumInitialCondition
	: public solid_dynamics::ElasticSolidDynamicsInitialCondition
{
public:
	MyocardiumInitialCondition(SolidBody *myocardium)
		: solid_dynamics::ElasticSolidDynamicsInitialCondition(myocardium) {};
protected:
	void Update(size_t index_particle_i, Real dt) override 
	{
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		if (base_particle_data_i.pos_n_[0] > 0.0) 
		{
			base_particle_data_i.vel_n_[1] = 5.0 * sqrt(3.0);
			base_particle_data_i.vel_n_[2] = 5.0;
		}
	};
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
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem system(Vecd(- SL - BW, -BW, -BW),
		Vecd(PL + BW, PH + BW, PH + BW), particle_spacing_ref);
	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body =
		new Myocardium(system, "MyocardiumBody", 0, ParticlesGeneratorOps::lattice);
	MyocardiumMuscle 	*muscle_material = new MyocardiumMuscle();
	ElasticSolidParticles 	particles(myocardium_body, muscle_material);
	/** Define Observer. */
	MyocardiumObserver *myocardium_observer 
		= new MyocardiumObserver(system, "MyocardiumObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles observer_particles(myocardium_observer);

	/** topology */
	SPHBodyInnerRelation* myocardium_body_inner = new SPHBodyInnerRelation(myocardium_body);
	SPHBodyContactRelation* myocardium_observer_contact = new SPHBodyContactRelation(myocardium_observer, { myocardium_body });

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Initialization. */
	MyocardiumInitialCondition initialization(myocardium_body);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_body_inner);
	/** Time step size calculation. */
	solid_dynamics::GetAcousticTimeStepSize 
		computing_time_step_size(myocardium_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_displacement("Displacement", in_output, myocardium_observer_contact);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	/** apply initial condition */
	initialization.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 3.0;
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
		Real integration_time = 0.0;
		while (integration_time < output_period) 
		{
			if (ite % 100 == 0) {
				cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
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
