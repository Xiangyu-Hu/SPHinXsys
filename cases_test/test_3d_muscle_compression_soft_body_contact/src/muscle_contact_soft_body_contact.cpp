/**
 * @file 	muscle_contact_soft_body_contact.cpp
 * @brief 	This is the test for muscle compression with our new contact model. Different particle resolutions are used for the two soft bodies that are in contact.
 * @author 	Chi Zhang and Xiangyu Hu, Bence Rochlitz
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real L  = 0.04;
Real PL = 0.1; 
Real resolution_ref = L / 12.0;
Real BW = resolution_ref * 4;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -0.5 * PL, -0.5 * PL),
	Vecd(2.0 * L + BW, 0.5 * PL, 0.5 * PL));

/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho_0 = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0; 
/** Define the geometry. */
TriangleMeshShape* CreateMyocardium()
{
	Vecd halfsize_myocardium(0.5 * L, 0.5 * L, 0.5 * L);
	Vecd translation_myocardium(0.5 * L, 0.0, 0.0);
	TriangleMeshShape* geometry_myocardium = new TriangleMeshShape(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
TriangleMeshShape* CreateStationaryPlate()
{
	Vecd halfsize_plate(0.5 * BW, 0.5 * L + BW , 0.5 * L + BW);
	Vecd translation_plate(-0.5 * BW, 0.0, 0.0);
	TriangleMeshShape* geometry_plate = new TriangleMeshShape(halfsize_plate,
		resolution, translation_plate);

	return geometry_plate;
}
/** Define the holder geometry. */
TriangleMeshShape* CreateMovingPlate()
{
	Vecd halfsize_plate(0.5 * BW,0.5*PL, 0.5*PL);
	Vecd translation_plate(L + BW, 0.0, 0.0);
	TriangleMeshShape* geometry_plate = new TriangleMeshShape(halfsize_plate,
		resolution, translation_plate);

	return geometry_plate;
}
/** Define the myocardium body. */
class Myocardium : public SolidBody
{
public:
	Myocardium(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addTriangleMeshShape(CreateMyocardium(), ShapeBooleanOps::add);
		body_shape_->addTriangleMeshShape(CreateStationaryPlate(), ShapeBooleanOps::add);
	}
};
/**
* @brief define the moving plate
*/
class MovingPlate : public SolidBody
{
public:
	MovingPlate(SPHSystem &system, std::string body_name, Real resolution_ref)
		: SolidBody(system, body_name, resolution_ref)
	{
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addTriangleMeshShape(CreateMovingPlate(), ShapeBooleanOps::add);
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
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateStationaryPlate(), ShapeBooleanOps::add);

		tagBodyPart();
	}
};
class HolderSpring : public BodyPartByParticle
{
public:
	HolderSpring(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateMovingPlate(), ShapeBooleanOps::add);
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
		E0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief Define moving plate material.
 */
class MovingPlateMaterial : public LinearElasticSolid
{
public:
	MovingPlateMaterial() : LinearElasticSolid()
	{
		rho0_ = rho_0;
		E0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
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
	ElasticSolidParticles 	myocardium_particles(myocardium_body, muscle_material);
	/** Plate. */
	MovingPlate *moving_plate = new MovingPlate(system, "MovingPlate", resolution_ref * 0.5);
	MovingPlateMaterial* moving_plate_material = new MovingPlateMaterial();
	ElasticSolidParticles 	moving_plate_particles(moving_plate, moving_plate_material);
	/** topology */
	InnerBodyRelation*   myocardium_body_inner = new InnerBodyRelation(myocardium_body);
	InnerBodyRelation*   moving_plate_inner = new InnerBodyRelation(moving_plate);
	
	SolidContactBodyRelation* myocardium_plate_contact = new SolidContactBodyRelation(myocardium_body, {moving_plate});
	SolidContactBodyRelation* plate_myocardium_contact = new SolidContactBodyRelation(moving_plate, {myocardium_body});
	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** initialize a time step */
	InitializeATimeStep 	myocardium_initialize_gravity(myocardium_body);
	Gravity* gravity = new Gravity(Vecd(-100.0, 0.0, 0.0));
	InitializeATimeStep plate_initialize_gravity(moving_plate, gravity);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form(myocardium_body_inner);
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form_2(moving_plate_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(myocardium_body);
	solid_dynamics::AcousticTimeStepSize computing_time_step_size_2(moving_plate);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half_2(moving_plate_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half_2(moving_plate_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation myocardium_update_contact_density(myocardium_plate_contact);
	solid_dynamics::ContactForce myocardium_compute_solid_contact_forces(myocardium_plate_contact);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation plate_update_contact_density(plate_myocardium_contact);
	solid_dynamics::ContactForce plate_compute_solid_contact_forces(plate_myocardium_contact);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Add spring contraint on the plate. */
	//solid_dynamics::ParticleWiseAcceleration spring_contraint(moving_plate, Vecd(1e6, 1e6, 1e6));
	// active_muscle_dynamics::SpringConstrainMuscleRegion
	// 	spring_contraint(moving_plate, new HolderSpring(moving_plate, "HolderSpring"));
	// spring_contraint.setUpSpringStiffness(Vec3d(1e10, 1e10, 1e10));
	/** Damping with the solid body*/
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		muscle_damping(myocardium_body_inner, 0.1, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		plate_damping(moving_plate_inner, 0.1, "Velocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);

	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	/** apply initial condition */
	corrected_configuration_in_strong_form.parallel_exec();
	corrected_configuration_in_strong_form_2.parallel_exec();
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.1;
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
			/** Gravity. */
			myocardium_initialize_gravity.parallel_exec();
			plate_initialize_gravity.parallel_exec();

			/** Contact model for myocardium. */
			myocardium_update_contact_density.parallel_exec();
			myocardium_compute_solid_contact_forces.parallel_exec();
			/** Contact model for plate. */
			plate_update_contact_density.parallel_exec();
			plate_compute_solid_contact_forces.parallel_exec();

			/** Stress relaxation and damping. */
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			stress_relaxation_first_half_2.parallel_exec(dt);
			//spring_contraint.parallel_exec(dt);
			plate_damping.parallel_exec(dt);
			//spring_contraint.parallel_exec(dt);
			stress_relaxation_second_half_2.parallel_exec(dt);

			ite++;
			dt = system.getSmallestTimeStepAmongSolidBodies();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			myocardium_body->updateCellLinkedList();
			moving_plate->updateCellLinkedList();

			myocardium_plate_contact->updateConfiguration();
			plate_myocardium_contact->updateConfiguration();
		}
		tick_count t2 = tick_count::now();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
