/**
* @file 	2d_plate.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider the body force apply on a 2D plate.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 10.0;									  /** Length of the square plate. */
Real PT = 1.0;									  /** Thickness of the square plate. */
Vec2d n_0 = Vec2d(0.0, 1.0);					  /** Pseudo-normal. */
int particle_number = 20;						  /** Particle number in the direction of the length */
Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 4;
Real BW = resolution_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -0.5 * resolution_ref),
								 Vec2d(PL + BW, 0.5 * resolution_ref));

/** For material properties of the solid. */
Real rho0_s = 1.0;				   /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;				   /** Poisson ratio. */
Real physical_viscosity = 200.0;   /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 100.0 * Youngs_modulus * 1.0e-4; /** Total distributed load. */
Real time_to_full_external_force = 0.1;

Real gravitational_acceleration = 0.0;

/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public ParticleGeneratorDirect
{
public:
	PlateParticleGenerator() : ParticleGeneratorDirect()
	{
		// the plate and boundary
		for (int i = 0; i < (particle_number + 2 * BWD); i++)
		{
			Real x = resolution_ref * i - BW + resolution_ref * 0.5;
			positions_volumes_.push_back(std::make_pair(Vecd(x, 0.0), resolution_ref));
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
	explicit PlateDynamicsInitialCondition(SolidBody &plate)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(plate){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial pseudo-normal. */
		n_0_[index_i] = n_0;
		n_[index_i] = n_0;
		pseudo_n_[index_i] = n_0_[index_i];
	};
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_->pos_n_[index_i][0] < 0.0 || base_particles_->pos_n_[index_i][0] > PL)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * define time dependent external force
 */
class TimeDependentExternalForce : public Gravity
{
public:
	explicit TimeDependentExternalForce(Vecd external_force)
		: Gravity(external_force) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_external_force ? current_time * global_acceleration_ / time_to_full_external_force : global_acceleration_;
	}
};

/** Define an observer particle generator. */
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(0.5 * PL, 0.0), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** Creat a plate body. */
	ThinStructure plate_body(system, "PlateBody", makeShared<SPHAdaptation>(1.15, 1.0));
	/** Creat particles for the elastic body. */
	ShellParticles plate_body_particles(plate_body,
										makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
										makeShared<PlateParticleGenerator>(), PT);

	/** Define Observer. */
	ObserverBody plate_observer(system, "PlateObserver");
	ObserverParticles observer_particles(plate_observer, makeShared<ObserverParticleGenerator>());

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner plate_body_inner(plate_body);
	BodyRelationContact plate_observer_contact(plate_observer, {&plate_body});

	/** Common particle dynamics. */
	TimeDependentExternalForce external_force(Vec2d(0.0, q / (PT * rho0_s) - gravitational_acceleration));
	TimeStepInitialization initialize_external_force(plate_body, external_force);

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** initial condition */
	PlateDynamicsInitialCondition plate_initial_velocity(plate_body);
	/** Corrected configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration(plate_body_inner);
	/** Time step size. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(plate_body);
	/** stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf stress_relaxation_first_half(plate_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf stress_relaxation_second_half(plate_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	thin_structure_dynamics::FixedFreeRotateShellBoundary constrain_holder(plate_body_inner, boundary_geometry);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		plate_position_damping(plate_body_inner, 0.5, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		plate_rotation_damping(plate_body_inner, 0.5, "AngularVelocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<indexVector, Vecd>>
		write_plate_max_displacement("Position", in_output, plate_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	plate_initial_velocity.parallel_exec();
	corrected_configuration.parallel_exec();

	/**
	* From here the time stepping begins.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_plate_max_displacement.writeToFile(0);

	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.5;
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
			if (ite % 1000 == 0)
			{
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

	write_plate_max_displacement.newResultTest();

	return 0;
}
