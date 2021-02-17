/**
 * @file 	muscle_compression.cpp
 * @brief 	This is the test for muscle compression with our new contact model. 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.3.0
 */

#include "base_data_package.h"
#include "all_kernels.h"
#include "all_particles.h"
#include "all_geometries.h"
#include "all_types_of_bodies.h"
#include "sph_system.h"
#include "all_materials.h"
#include "all_physical_dynamics.h"

#include "xml_engine.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/Force_BuiltIns.h"

#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/CableTrackerSubsystem.h"


#include "in_output.h"
/** Standrad c++ libraries. */
#include <iostream>



/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real L  = 0.04;
Real PL = 0.1; 
Real particle_spacing_ref = L / 12.0;
Real BW = particle_spacing_ref * 4;
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
	Myocardium(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
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
	MovingPlate(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
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
	Holder(SolidBody *solid_body, string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateStationaryPlate(), ShapeBooleanOps::add);

		tagBodyPart();
	}
};
/**
* @brief body part constrainted by multi-body dynamics
*/
class SimbodyPlate : public SolidBodyPartForSimbody
{
public:
	SimbodyPlate(SolidBody* solid_body,
		string constrained_region_name, Real solid_body_density)
		: SolidBodyPartForSimbody(solid_body,constrained_region_name)
	{
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateMovingPlate(), ShapeBooleanOps::add);
		/** tag the constrained particle. */
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
		rho_0_ 	= rho_0;
		E_0_ = Youngs_modulus;
		nu_ = poisson;
		eta_0_ = physical_viscosity;

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
		rho_0_ = rho_0;
		E_0_ = Youngs_modulus;
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
	SPHSystem system(Vecd(-BW, -0.5 * PL, - 0.5 * PL),
		Vecd(2.0 * L + BW, 0.5 * PL, 0.5 * PL), particle_spacing_ref);
	/** Creat a Myocardium body, corresponding material, particles and reaction model. */
	Myocardium *myocardium_body = new Myocardium(system, "MyocardiumBody", 0);
	MyocardiumMuscle 	*muscle_material = new MyocardiumMuscle();
	ElasticSolidParticles 	myocardium_particles(myocardium_body, muscle_material);
	/** Plate. */
	MovingPlate *moving_plate = new MovingPlate(system, "MovingPlate",	0);
	MovingPlateMaterial* moving_plate_material = new MovingPlateMaterial();
	SolidParticles 	moving_plate_particles(moving_plate, moving_plate_material);
	/** topology */
	SPHBodyInnerRelation*   myocardium_body_inner = new SPHBodyInnerRelation(myocardium_body);
	SolidBodyContactRelation* myocardium_plate_contact = new SolidBodyContactRelation(myocardium_body, {moving_plate});
	SolidBodyContactRelation* plate_myocardium_contact = new SolidBodyContactRelation(moving_plate, {myocardium_body});
	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** initialize a time step */
	InitializeATimeStep 	myocardium_initialize_gravity(myocardium_body);
	/** Corrected strong configuration. */	
	solid_dynamics::CorrectConfiguration corrected_configuration_in_strong_form(myocardium_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(myocardium_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(myocardium_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_body_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::SummationContactDensity myocardium_update_contact_density(myocardium_plate_contact);
	solid_dynamics::ContactForce myocardium_compute_solid_contact_forces(myocardium_plate_contact);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::SummationContactDensity plate_update_contact_density(plate_myocardium_contact);
	solid_dynamics::ContactForce plate_compute_solid_contact_forces(plate_myocardium_contact);
	/** Constrain the holder. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(myocardium_body, new Holder(myocardium_body, "Holder"));
	/** Damping with the solid body*/
	DampingBySplittingWithRandomChoice<SPHBodyInnerRelation, DampingBySplittingPairwise<Vec3d>, Vec3d>
		muscle_damping(myocardium_body_inner, 0.1, myocardium_particles.vel_n_, physical_viscosity);
	/** Output */
	In_Output in_output(system);
	WriteBodyStatesToVtu write_states(in_output, system.real_bodies_);
	/** Simbody interface. */
	/**
	* The multi body system from simbody.
	*/
	SimTK::MultibodySystem          MBsystem;
	/** The bodies or matter of the MBsystem. */
	SimTK::SimbodyMatterSubsystem   matter(MBsystem);
	/** The forces of the MBsystem.*/
	SimTK::GeneralForceSubsystem    forces(MBsystem);
	SimTK::CableTrackerSubsystem    cables(MBsystem);
	/** mass proeprties of the fixed spot. */
	SimbodyPlate*  	plate_multibody = new SimbodyPlate(moving_plate, "Plate", rho_0);
	/** Mass properties of the consrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid      rigid_info(*plate_multibody->body_part_mass_properties_);
	SimTK::MobilizedBody::Slider 
		plateMBody(matter.Ground(), SimTK::Transform(SimTK::Vec3(0)), rigid_info, SimTK::Transform(SimTK::Vec3(0)));
	/** Gravity. */
	SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(Real(-100.0), 0.0, 0.0));
	/** discreted forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	/** Damper. */
	SimTK::Force::MobilityLinearDamper linear_damper(forces, plateMBody, SimTK::MobilizerUIndex(0), 20.0);
	/** Time steping method for multibody system.*/
	SimTK::State state = MBsystem.realizeTopology();
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);
	/** Coupling between SimBody and SPH.*/
	solid_dynamics::TotalForceOnSolidBodyPartForSimBody
		force_on_plate(moving_plate, plate_multibody, MBsystem, plateMBody, force_on_bodies, integ);
	solid_dynamics::ConstrainSolidBodyPartBySimBody
		constraint_plate(moving_plate, plate_multibody, MBsystem, plateMBody, force_on_bodies, integ);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	/** apply initial condition */
	corrected_configuration_in_strong_form.parallel_exec();
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
				cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			/** Gravity. */
			myocardium_initialize_gravity.parallel_exec();
			/** Contact model for myocardium. */
			myocardium_update_contact_density.parallel_exec();
			myocardium_compute_solid_contact_forces.parallel_exec();
			/** Contact model for plate. */
			plate_update_contact_density.parallel_exec();
			plate_compute_solid_contact_forces.parallel_exec();
			{
				SimTK::State& state_for_update = integ.updAdvancedState();
				force_on_bodies.clearAllBodyForces(state_for_update);
				force_on_bodies.setOneBodyForce(state_for_update, plateMBody, force_on_plate.parallel_exec());
				integ.stepBy(dt);
				constraint_plate.parallel_exec();
			}
			/** Stress relaxation and damping. */
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
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
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
