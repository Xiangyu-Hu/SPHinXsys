/**
 * @file passive_cantilever.cpp
 * @brief This is the first example of myocardium 
 * @author Chi Zhang and Xiangyu Hu
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
Real resolution_ref = PH / 12.0;		/**< Initial particle spacing. */
Real BW = resolution_ref * 4; 		/**< Boundary width. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-SL - BW, -BW, -BW),
	Vecd(PL + BW, PH + BW, PH + BW));
											
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
TriangleMeshShape* createMyocardium()
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
	Myocardium(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		body_shape_ = new ComplexShape(mesh);
		mesh->addTriangleMeshShape(createMyocardium(), ShapeBooleanOps::add);
		mesh->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);
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
	Holder(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		body_part_shape_ = new ComplexShape(mesh);
		mesh->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);

		tagBodyPart();
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
		rho0_ 	= rho_0;
		bulk_modulus_ = bulk_modulus;
		f0_ 	= fiber_direction;
		s0_ 	= sheet_direction;
		std::copy(a_0, a_0 + 4, a0_);
		std::copy(b_0, b_0 + 4, b0_);

		assignDerivedMaterialParameters();
	}
};
/**
 * application dependent initial condition 
 */
class MyocardiumInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	MyocardiumInitialCondition(SolidBody *myocardium)
		: solid_dynamics::ElasticDynamicsInitialCondition(myocardium) {};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		if (pos_n_[index_i][0] > 0.0)
		{
			vel_n_[index_i][1] = 5.0 * sqrt(3.0);
			vel_n_[index_i][2] = 5.0;
		}
	};
};
 //define an observer body
class MyocardiumObserver : public FictitiousBody
{
public:
	MyocardiumObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		body_input_points_volumes_.push_back(std::make_pair(Vecd(PL, PH, PW), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body = new Myocardium(system, "MyocardiumBody");
	MyocardiumMuscle 	*muscle_material = new MyocardiumMuscle();
	ElasticSolidParticles 	particles(myocardium_body, muscle_material);
	/** Define Observer. */
	MyocardiumObserver *myocardium_observer = new MyocardiumObserver(system, "MyocardiumObserver");
	BaseParticles observer_particles(myocardium_observer);

	/** topology */
	BodyRelationInner* myocardium_body_inner = new BodyRelationInner(myocardium_body);
	BodyRelationContact* myocardium_observer_contact = new BodyRelationContact(myocardium_observer, { myocardium_body });

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Initialization. */
	MyocardiumInitialCondition initialization(myocardium_body);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize 
		computing_time_step_size(myocardium_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd> 
		write_displacement("Position", in_output, myocardium_observer_contact);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	/** apply initial condition */
	initialization.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.writeToFile(0);
	write_displacement.writeToFile(0);
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
				std::cout << "N=" << ite << " Time: "
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
		write_displacement.writeToFile(ite);
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
