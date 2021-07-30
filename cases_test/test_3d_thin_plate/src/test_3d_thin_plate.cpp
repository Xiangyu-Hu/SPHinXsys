/**
* @file 	test_3d_thin_plate.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider the body force applied on a quasi-static square plate.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
*/
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 10.0;                                          /** Length of the square plate. */
Real PH = 10.0;                                          /** Width of the square plate. */
Real PT = 1.0;                                           /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                        /** Pseudo-normal. */
int particle_number = 14;								 /** Particle number in the direction of the length */
Real resolution_ref = PL / (Real)particle_number;        /** Initial reference particle spacing. */
int BWD = 1;
Real BW = resolution_ref * (Real)BWD;                    /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
	Vec3d(PL + BW, PH + BW, 0.5 * resolution_ref));

/** For material properties of the solid. */
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 1.3024653e6;	                     /** Normalized Youngs Modulus. */
Real poisson = 0.3; 			                         /** Poisson ratio. */
Real physical_viscosity = 200.0;                         /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 100.0 * Youngs_modulus * 1.0e-4;                /** Total distributed load. */
Real time_to_full_external_force = 0.1;

Real gravitational_acceleration = 0.009646;

class Plate : public ThinStructure
{
public:
	Plate(SPHSystem &system, std::string body_name, ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: ThinStructure(system, body_name, particle_adaptation, particle_generator)
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			for (int j = 0; j < (particle_number + 2 * BWD); j++)
			{
				Real x = resolution_ref * i - BW + resolution_ref * 0.5;
				Real y = resolution_ref * j - BW + resolution_ref * 0.5;
				body_input_points_volumes_.push_back(std::make_pair(Vecd(x, y, 0.0), resolution_ref * resolution_ref));
			}
		}
	}
};

/**
 * application dependent initial condition
 */
class PlateDynamicsInitialCondition
	: public thin_structure_dynamics::ShellDynamicsInitialCondition
{
public:
	PlateDynamicsInitialCondition(SolidBody *plate)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(plate) {};
protected:
	void Update(size_t index_i, Real dt) override {
		/** initial pseudo-normal. */
		n_0_[index_i] = n_0;
		n_[index_i] = n_0;
		pseudo_n_[index_i] = n_0_[index_i];
	};
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
protected:
	virtual void tagBodyPart() override
	{
		BaseParticles* base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			if (base_particles->pos_n_[i][0] < 0.0 || base_particles->pos_n_[i][1] < 0.0
				|| base_particles->pos_n_[i][0] > PL || base_particles->pos_n_[i][1] > PH)
			{
				tagAParticle(i);
			}
		}
	};
public:
	BoundaryGeometry(SPHBody *body, std::string body_part_name)
		: BodyPartByParticle(body, body_part_name) {
		tagBodyPart();
	};
	virtual ~BoundaryGeometry() {};
};


/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force ?
			current_time * global_acceleration_ / time_to_full_external_force : global_acceleration_;
	}
};


/** Define an observer body. */
class PlateObserver : public FictitiousBody
{
public:
	PlateObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(std::make_pair(Vecd(0.5 * PL, 0.5 * PH, 0.0), 0.0));
	}
};

class PlateMaterial : public LinearElasticSolid
{
public:
	PlateMaterial(): LinearElasticSolid()
	{
		rho0_ = rho0_s;
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
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** Define the external force. */
	TimeDependentExternalForce external_force(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration));

	/** Creat a plate body. */
	Plate *plate_body = new Plate(system, "PlateBody", new ParticleAdaptation(1.15, 1.0), new ParticleGeneratorDirect());
	/** elastic solid material properties */
	PlateMaterial *plate_material = new PlateMaterial();
	/** Creat particles for the elastic body. */
	ShellParticles plate_body_particles(plate_body, plate_material, PT);

	/** Define Observer. */
	PlateObserver *plate_observer = new PlateObserver(system, "PlateObserver");
	BaseParticles observer_particles(plate_observer);

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner* plate_body_inner = new BodyRelationInner(plate_body);
	BodyRelationContact* plate_observer_contact = new BodyRelationContact(plate_observer, { plate_body });

	/** Common particle dynamics. */
	TimeStepInitialization 	initialize_external_force(plate_body, &external_force);

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** initial condition */
	PlateDynamicsInitialCondition plate_initial_pseudo_normal(plate_body);
	 /** Corrected strong configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration_in_strong_form(plate_body_inner);
	/** Time step size calculation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	/** active-passive stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_holder(plate_body, new BoundaryGeometry(plate_body, "BoundaryGeometry"));
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		plate_position_damping(plate_body_inner, 0.5, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		plate_rotation_damping(plate_body_inner, 0.5, "AngularVelocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		write_plate_max_displacement("Position", in_output, plate_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	plate_initial_pseudo_normal.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_plate_max_displacement.writeToFile(0);
	
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.8;
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
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			plate_position_damping.parallel_exec(dt);
			plate_rotation_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

		}
		write_plate_max_displacement.writeToFile(ite);
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

