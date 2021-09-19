/**
 * @file 	muscle_contact_soft_body_contact.cpp
 * @brief 	This is the test for muscle compression with our new contact model. 
 * Different particle resolutions are used for the two soft bodies that are in contact.
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
TriangleMeshShape* createMyocardium()
{
	Vecd halfsize_myocardium(0.5 * L, 0.5 * L, 0.5 * L);
	Vecd translation_myocardium(0.5 * L, 0.0, 0.0);
	TriangleMeshShape* geometry_myocardium = new TriangleMeshShape(halfsize_myocardium, resolution,
		translation_myocardium);

	return geometry_myocardium;
}
/** Define the holder geometry. */
TriangleMeshShape* createStationaryPlate()
{
	Vecd halfsize_plate(0.5 * BW, 0.5 * L + BW , 0.5 * L + BW);
	Vecd translation_plate(-0.5 * BW, 0.0, 0.0);
	TriangleMeshShape* geometry_plate = new TriangleMeshShape(halfsize_plate,
		resolution, translation_plate);

	return geometry_plate;
}
/** Define the holder geometry. */
TriangleMeshShape* createMovingPlate()
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
		mesh_.reset(new ComplexShapeTriangleMesh());
		body_shape_ = new ComplexShape(mesh_.get());
		mesh_->addTriangleMeshShape(createMyocardium(), ShapeBooleanOps::add);
		mesh_->addTriangleMeshShape(createStationaryPlate(), ShapeBooleanOps::add);
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
};
/**
* @brief define the moving plate
*/
class MovingPlate : public SolidBody
{
public:
	MovingPlate(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 1.5))
	{
		mesh_.reset(new ComplexShapeTriangleMesh());
		body_shape_ = new ComplexShape(mesh_.get());
		mesh_->addTriangleMeshShape(createMovingPlate(), ShapeBooleanOps::add);
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
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
		mesh_.reset(new ComplexShapeTriangleMesh());
		body_part_shape_ = new ComplexShape(mesh_.get());
		mesh_->addTriangleMeshShape(createStationaryPlate(), ShapeBooleanOps::add);

		tagBodyPart();
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
};
class HolderSpring : public BodyPartByParticle
{
public:
	HolderSpring(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		mesh_.reset(new ComplexShapeTriangleMesh());
		body_part_shape_ = new ComplexShape(mesh_.get());
		mesh_->addTriangleMeshShape(createMovingPlate(), ShapeBooleanOps::add);
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
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
 * @brief Define moving plate material.
 */
class MovingPlateMaterial : public NeoHookeanSolid
{
public:
	MovingPlateMaterial() : NeoHookeanSolid()
	{
		rho0_ = rho_0;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

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
	MovingPlate *moving_plate = new MovingPlate(system, "MovingPlate");
	MovingPlateMaterial* moving_plate_material = new MovingPlateMaterial();
	ElasticSolidParticles 	moving_plate_particles(moving_plate, moving_plate_material);
	/** topology */
	BodyRelationInner*   myocardium_body_inner = new BodyRelationInner(myocardium_body);
	BodyRelationInner*   moving_plate_inner = new BodyRelationInner(moving_plate);
	
	SolidBodyRelationContact* myocardium_plate_contact = new SolidBodyRelationContact(myocardium_body, {moving_plate});
	SolidBodyRelationContact* plate_myocardium_contact = new SolidBodyRelationContact(moving_plate, {myocardium_body});
	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** initialize a time step */
	TimeStepInitialization 	myocardium_initialize_gravity(myocardium_body);
	Gravity* gravity = new Gravity(Vecd(-100.0, 0.0, 0.0));
	TimeStepInitialization plate_initialize_gravity(moving_plate, gravity);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form(myocardium_body_inner);
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form_2(moving_plate_inner);
	/** active and passive stress relaxation. */
	solid_dynamics::KirchhoffStressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	solid_dynamics::KirchhoffStressRelaxationFirstHalf stress_relaxation_first_half_2(moving_plate_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half_2(moving_plate_inner);
	//stress_relaxation_first_half_2.post_processes_(spring_constraint);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation myocardium_update_contact_density(myocardium_plate_contact);
	solid_dynamics::ContactForce myocardium_compute_solid_contact_forces(myocardium_plate_contact);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ContactDensitySummation plate_update_contact_density(plate_myocardium_contact);
	solid_dynamics::ContactForce plate_compute_solid_contact_forces(plate_myocardium_contact);

	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Add spring constraint on the plate. */
	solid_dynamics::SpringDamperConstraintParticleWise spring_constraint(moving_plate, Vecd(0.2, 0, 0), 0.01);

	/** Damping with the solid body*/
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		muscle_damping(myocardium_body_inner, 0.2, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		plate_damping(moving_plate_inner, 0.2, "Velocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_states(in_output, system.real_bodies_);

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
	write_states.writeToFile(0);
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

			spring_constraint.parallel_exec();

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
			plate_damping.parallel_exec(dt);
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
