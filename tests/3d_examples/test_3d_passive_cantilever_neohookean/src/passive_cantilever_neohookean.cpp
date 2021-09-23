/**
 * @file passive_cantilever_neohookean.cpp
 * @brief This is the example of myocardium with simple neohookean tissue model 
 * @author Bence Rochlitz, Chi Zhang  and Xiangyu Hu
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
Real resolution_ref = PH / 6.0;		/**< Initial particle spacing. */
Real BW = resolution_ref * 4; 		/**< Boundary width. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-SL, 0, 0), Vecd(PL, PH, PH));

/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 1265.0; // Gheorghe 2019 
Real poisson = 0.45; // nearly incompressible
Real Youngs_modulus = 5e4; // Sommer 2015
Real physical_viscosity = 50.0; //physical damping, here we choose the same value as numerical viscosity
Real gravity_g = 9.8; 					/**< Value of gravity. */
Real time_to_full_gravity = 0.0;

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
class MyocardiumMuscle : public NeoHookeanSolid
{
public:
	MyocardiumMuscle() : NeoHookeanSolid()
	{
		rho0_ 	= rho_0;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

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
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** Define the external force. */
	TimeDependentGravity gravity(Vec3d(0.0, -gravity_g, 0.0));

	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body = new Myocardium(system, "MyocardiumBody");
	MyocardiumMuscle 	*muscle_material = new MyocardiumMuscle();
	ElasticSolidParticles 	myocardium_particles(myocardium_body, muscle_material);
	/** Define Observer. */
	MyocardiumObserver *myocardium_observer = new MyocardiumObserver(system, "MyocardiumObserver");
	BaseParticles observer_particles(myocardium_observer);

	/** topology */
	BodyRelationInner* myocardium_body_inner = new BodyRelationInner(myocardium_body);
	BodyRelationContact* myocardium_observer_contact = new BodyRelationContact(myocardium_observer, { myocardium_body });

	//-------- common particle dynamics ----------------------------------------
	TimeStepInitialization 	initialize_gravity(myocardium_body, &gravity);

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration 
		corrected_configuration_in_strong_form(myocardium_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize 
		computing_time_step_size(myocardium_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(myocardium_body_inner);
	/** Setup the damping stress, if you know what you are doing. */
	//stress_relaxation_first_step.setupDampingStressFactor(1.0);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(myocardium_body_inner);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	DampingWithRandomChoice<DampingBySplittingInner<indexVector, Vec3d>>
		muscle_damping(myocardium_body_inner, 0.1, "Velocity", physical_viscosity);
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
	corrected_configuration_in_strong_form.parallel_exec();
	write_states.writeToFile(0);
	write_displacement.writeToFile(0);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 1.0;
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

			initialize_gravity.parallel_exec(); // gravity force
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
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
