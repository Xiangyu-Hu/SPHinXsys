/**
 * @file passive_cantilever_neohookean.cpp
 * @brief This is the example of myocardium with simple neohookean tissue model 
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
Real physical_viscosity = 50.0; //physical damping, here we choose the same value as numerical viscosity
Real gravity_g = 9.8; 					/**< Value of gravity. */
Real time_to_full_gravity = 0.0;

Vec3d d_0(1.0e-4, 0.0, 0.0);

/**
 * define time dependent gravity
 */
class TimeDependentGravity : public Gravity //*** move it to source
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

/** Define the geometry. */
TriangleMeshShape* CreateBeamShape() //*** move it to source
{
	Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
	Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
	TriangleMeshShape* geometry_myocardium = new TriangleMeshShape(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}

/** Define the holder geometry. */
TriangleMeshShape* CreateHolder() //*** move it to source
{
	Vecd halfsize_holder(0.5*SL, 0.5*PH , 0.5*PW);
	Vecd translation_holder(-0.5*SL, 0.5 * PH, 0.5 * PW);
	TriangleMeshShape* geometry_holder = new TriangleMeshShape(halfsize_holder,
		resolution, translation_holder);

	return geometry_holder;
}

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
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);

		tagBodyPart();
	}
};

/** define an observer body */
class MyocardiumObserver : public FictitiousBody
{
public:
	MyocardiumObserver(SPHSystem &system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		body_input_points_volumes_.push_back(make_pair(Point(PL, PH, PW), 0.0));
	}
};

class PassiveCantilever
{
public:
	PassiveCantilever(){}
	~PassiveCantilever(){}
	/** Setup the system. */
    SPHSystem* system;
	/** Define the external force. */
    TimeDependentGravity* gravity;
	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
    SolidBody *myocardium_body;
    NeoHookeanSolid *muscle_material;
    ElasticSolidParticles *myocardium_particles;
	/** Define Observer. */
    MyocardiumObserver *myocardium_observer;
    BaseParticles *observer_particles;
	/** topology */
    SPHBodyInnerRelation* myocardium_body_inner;
    SPHBodyContactRelation* myocardium_observer_contact;

	//-------- common particle dynamics ----------------------------------------
	InitializeATimeStep *initialize_gravity;
	/** This section define all numerical methods will be used in this case. */
	solid_dynamics::CorrectConfiguration *corrected_configuration_in_strong_form;
	solid_dynamics::AcousticTimeStepSize *computing_time_step_size;
	solid_dynamics::StressRelaxationFirstHalf *stress_relaxation_first_half;
	solid_dynamics::StressRelaxationSecondHalf *stress_relaxation_second_half;
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion *constrain_holder;
	DampingBySplittingWithRandomChoice<SPHBodyInnerRelation, DampingBySplittingPairwise<Vec3d>, Vec3d> *muscle_damping;

	/** Output */
	/*In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_displacement("Displacement", in_output, myocardium_observer_contact);*/

	////-----------main functions-----------////
	Real dt;
    void initialize_simulation();
	void perform_simulation_step();

};

void PassiveCantilever::initialize_simulation(){
    /** Setup the system. */
    system = new SPHSystem (Vecd(-SL, 0, 0), Vecd(PL, PH, PH), particle_spacing_ref);
	/** Define the external force. */
    gravity = new TimeDependentGravity (Vec3d(0.0, -gravity_g, 0.0));
    /** Define the myocardium body. */
    myocardium_body = new SolidBody(system, "MyocardiumBody", 0);
    myocardium_body->body_shape_ = new ComplexShape(body_name);
	myocardium_body->body_shape_->addTriangleMeshShape(CreateBeamShape(), ShapeBooleanOps::add); // soft body
	myocardium_body->body_shape_->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add); // wall fixation
    /** Setup material properties of myocardium **/
    muscle_material =  new NeoHookeanSolid(system, "MyocardiumBody", 0); //*** need new constructor for NeoHookeanSolid with material paremeters
    muscle_material->E_0_ = Youngs_modulus;
    muscle_material->nu_ = poisson;
    muscle_material->eta_0_ = physical_viscosity;
    muscle_material->assignDerivedMaterialParameters();
    /** myocardium_particles */
    myocardium_particles = new ElasticSolidParticles(myocardium_body, muscle_material);
    /** Define Observer. */
    myocardium_observer = new MyocardiumObserver(system, "MyocardiumObserver", 0);
    observer_particles = new BaseParticles(myocardium_observer);
    /** topology */
    myocardium_body_inner = new SPHBodyInnerRelation(myocardium_body);
    myocardium_observer_contact = new SPHBodyContactRelation(myocardium_observer, { myocardium_body });

	//--------common particle dynamics--------//
	initialize_gravity = new InitializeATimeStep(myocardium_body, &gravity);

	/** This section define all numerical methods will be used in this case. */
	/** Corrected strong configuration. */	
	corrected_configuration_in_strong_form = new solid_dynamics::CorrectConfiguration(myocardium_body_inner);
	/** Time step size calculation. */
	computing_time_step_size = new solid_dynamics::AcousticTimeStepSize(myocardium_body);
	/** active and passive stress relaxation. */
	stress_relaxation_first_half = new solid_dynamics::StressRelaxationFirstHalf(myocardium_body_inner);
	/** Setup the damping stress, if you know what you are doing. */
	stress_relaxation_second_half = new solid_dynamics::StressRelaxationSecondHalf (myocardium_body_inner);
	/** Constrain the holder. */
	constrain_holder = new solid_dynamics::ConstrainSolidBodyRegion(myocardium_body, new Holder(myocardium_body, "Holder"));
	muscle_damping = new DampingBySplittingWithRandomChoice<SPHBodyInnerRelation, DampingBySplittingPairwise<Vec3d>, Vec3d>(myocardium_body_inner, 0.1, myocardium_particles.vel_n_, physical_viscosity);
	
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration_in_strong_form.parallel_exec();
	/*write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_displacement.WriteToFile(GlobalStaticVariables::physical_time_);*/

	dt = 0.0; 
}

void PassiveCantilever::perform_simulation_step(){

	initialize_gravity.parallel_exec(); // gravity force
	stress_relaxation_first_half.parallel_exec(dt);
	constrain_holder.parallel_exec(dt);
	muscle_damping.parallel_exec(dt);
	constrain_holder.parallel_exec(dt);
	stress_relaxation_second_half.parallel_exec(dt);

	dt = computing_time_step_size.parallel_exec();
}

