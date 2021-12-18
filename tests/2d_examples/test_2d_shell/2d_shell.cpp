/**
* @file 	dw_2d_shell.cpp
* @brief 	This is the benchmark test of the 2d shell.
* @details  We consider large deformation of a 2D shell.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 24.5;											   /** Radius of the inner boundary of the cylinder. */
Real thickness = 1.0;										   /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;			   /** Radius of the mid surface. */
int particle_number = 4;									   /** Particle number in the thickness. */
Real particle_spacing_ref = thickness / (Real)particle_number; /** Initial reference particle spacing. */
int particle_number_mid_surface = 2 * radius_mid_surface * Pi * 80.0 / 360.0 / particle_spacing_ref;
int BWD = 4;
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-radius - thickness, 0.0),
								 Vec2d(radius + thickness, radius + thickness));

/** For material properties of the solid. */
Real rho0_s = 3.67346939;		  /** Normalized density. */
Real Youngs_modulus = 4.32e7;	  /** Normalized Youngs Modulus. */
Real poisson = 0.3;				  /** Poisson ratio. */
Real physical_viscosity = 1000.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.1;
Real gravitational_acceleration = -1500.0;
/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public ParticleGeneratorDirect
{
public:
	CylinderParticleGenerator() : ParticleGeneratorDirect()
	{
		// the cylinder and boundary
		for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
		{
			Real x = radius_mid_surface *
					 cos(50.0 / 180.0 * Pi + (i + 0.5 - BWD) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
			Real y = radius_mid_surface *
					 sin(50.0 / 180.0 * Pi + (i + 0.5 - BWD) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
			positions_volumes_.push_back(std::make_pair(Vecd(x, y), particle_spacing_ref));
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
		n_0_[index_i] = Vec2d(pos_0_[index_i][0] / radius_mid_surface, pos_0_[index_i][1] / radius_mid_surface);
		n_[index_i] = n_0_[index_i];
		pseudo_n_[index_i] = n_0_[index_i];
		transformation_matrix_[index_i] = getTransformationMatrix(n_0_[index_i]);
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
		if (base_particles_->pos_n_[index_i][0] < -radius_mid_surface * cos(50.0 / 180.0 * Pi) || base_particles_->pos_n_[index_i][0] > radius_mid_surface * cos(50.0 / 180.0 * Pi))
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

/** Define an observer body. */
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(0.0, radius_mid_surface), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Create a Cylinder body. */
	ThinStructure cylinder_body(system, "CylinderBody", makeShared<SPHAdaptation>(1.15, 1.0));
	/** Create particles for the elastic body. */
	ShellParticles cylinder_body_particles(cylinder_body,
										   makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
										   makeShared<CylinderParticleGenerator>(), thickness);
	cylinder_body_particles.addAVariableToWrite<indexVector, Vecd>("PseudoNormal");

	/** Define Observer. */
	ObserverBody cylinder_observer(system, "CylinderObserver");
	ObserverParticles observer_particles(cylinder_observer, makeShared<ObserverParticleGenerator>());

	/** Set body contact map
	 *  The contact map gives the data connections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner cylinder_body_inner(cylinder_body);
	BodyRelationContact cylinder_observer_contact(cylinder_observer, {&cylinder_body});

	/** Common particle dynamics. */
	TimeDependentExternalForce external_force(Vec2d(0.0, gravitational_acceleration));
	TimeStepInitialization initialize_external_force(cylinder_body, external_force);

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/** initial condition */
	PlateDynamicsInitialCondition cylinder_initial_pseudo_normal(cylinder_body);
	/** Corrected configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration(cylinder_body_inner);
	/** Time step size calculation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(cylinder_body);
	/** stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(cylinder_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(cylinder_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(cylinder_body, "BoundaryGeometry");
	thin_structure_dynamics::ClampConstrainShellBodyRegion
		fixed_free_rotate_shell_boundary(cylinder_body_inner, boundary_geometry);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		cylinder_position_damping(cylinder_body_inner, 0.5, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		cylinder_rotation_damping(cylinder_body_inner, 0.5, "AngularVelocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<indexVector, Vecd>>
		write_cylinder_max_displacement("Position", in_output, cylinder_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	cylinder_initial_pseudo_normal.parallel_exec();
	corrected_configuration.parallel_exec();

	/**
	* From here the time stepping begins.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_cylinder_max_displacement.writeToFile(0);

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
		Real integral_time = 0.0;
		while (integral_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			fixed_free_rotate_shell_boundary.parallel_exec(dt);
			cylinder_position_damping.parallel_exec(dt);
			cylinder_rotation_damping.parallel_exec(dt);
			fixed_free_rotate_shell_boundary.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_cylinder_max_displacement.writeToFile(ite);
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_cylinder_max_displacement.newResultTest();

	return 0;
}
