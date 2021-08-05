/**
* @file 	dw_3d_cylindrical_surface.cpp
* @brief 	This is the benchmark test of the shell.
* @details  We consider large deformation of a cylindrical surface.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 24.875;                                     /** Radius of the inner boudary of the cylinder. */
Real height = 50.0;                                     /** Height of the cylinder. */
Real thickness = 0.25;                                  /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
int particle_number = 40;								 /** Particle number in the peripheral direction. */
Real particle_spacing_ref = 2 * radius_mid_surface * Pi * 80.0 / 360.0 / (Real)particle_number;     /** Initial reference particle spacing. */
int BWD = 1;
Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, 0.0, -radius - thickness),
	Vec3d(radius + thickness + BW, height, radius + thickness));

/** For material properties of the solid. */
Real rho0_s = 36.0; 			                     /** Normalized density. */
Real Youngs_modulus = 4.32e8;	                         /** Normalized Youngs Modulus. */
Real poisson = 0.0; 			                         /** Poisson ratio. */
Real physical_viscosity = 1.0e4;                         /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.1;
Real gravitational_acceleration = -10.0;

class Cylinder : public ThinStructure
{
public:
	Cylinder(SPHSystem &system, std::string body_name, ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: ThinStructure(system, body_name, particle_adaptation, particle_generator)
	{
		// the cylinder and boundary
		for (int i = 0; i < particle_number; i++)
		{
			for (int j = 0; j < (height / particle_spacing_ref + 2 * BWD - 1); j++)
			{
				Real x = radius_mid_surface * cos(50.0 / 180.0 * Pi + (i + 0.5) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number);
				Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
				Real z = radius_mid_surface * sin(50.0 / 180.0 * Pi + (i + 0.5) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number);
				body_input_points_volumes_.push_back(std::make_pair(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref));
			}
		}
	}
};

/**
 * application dependent initial condition
 */
class CylinderDynamicsInitialCondition
	: public thin_structure_dynamics::ShellDynamicsInitialCondition
{
public:
	CylinderDynamicsInitialCondition(SolidBody *plate)
		: thin_structure_dynamics::ShellDynamicsInitialCondition(plate) {};
protected:
	void Update(size_t index_i, Real dt) override {
		/** initial pseudo-normal. */
		n_0_[index_i] = Vec3d(pos_0_[index_i][0] / radius_mid_surface, 0.0, pos_0_[index_i][2] / radius_mid_surface);
		n_[index_i] = n_0_[index_i];
		pseudo_n_[index_i] = n_0_[index_i];
		transformation_matrix_[index_i] = getTransformationMatrix(n_0_[index_i]);
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
			if (base_particles->pos_n_[i][1] < 0.0 || base_particles->pos_n_[i][1] > height - 0.5 * particle_spacing_ref)
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
class CylinderObserver : public FictitiousBody
{
public:
	CylinderObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(std::make_pair(Vecd(radius_mid_surface * cos(5.0 / 18.0 * Pi), 0.5 * height, radius_mid_surface * sin(5.0 / 18.0 * Pi)), 0.0));
	}
};

class CylinderMaterial : public LinearElasticSolid
{
public:
	CylinderMaterial(): LinearElasticSolid()
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
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Define the external force. */
	TimeDependentExternalForce external_force(Vec3d(0.0, 0.0, gravitational_acceleration));

	/** Creat a Cylinder body. */
	Cylinder *cylinder_body = new Cylinder(system, "CylinderBody", new ParticleAdaptation(1.15, 1.0), new ParticleGeneratorDirect());
	/** elastic soild material properties */
	CylinderMaterial *cylinder_material = new CylinderMaterial();
	/** Creat particles for the elastic body. */
	ShellParticles cylinder_body_particles(cylinder_body, cylinder_material, thickness);
	cylinder_body_particles.addAVariableToWrite<indexMatrix, Mat3d>("FirstPiolaKirchhoffStress");

	/** Define Observer. */
	CylinderObserver *cylinder_observer = new CylinderObserver(system, "CylinderObserver");
	BaseParticles observer_particles(cylinder_observer);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner* cylinder_body_inner = new BodyRelationInner(cylinder_body);
	BodyRelationContact* cylinder_observer_contact = new BodyRelationContact(cylinder_observer, { cylinder_body });

	/** Common particle dynamics. */
	TimeStepInitialization 	initialize_external_force(cylinder_body, &external_force);

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** initial condition */
	CylinderDynamicsInitialCondition cylinder_initial_pseudo_normal(cylinder_body);
	 /** Corrected strong configuration. */
	thin_structure_dynamics::ShellCorrectConfiguration
		corrected_configuration_in_strong_form(cylinder_body_inner);
	/** Time step size caclutation. */
	thin_structure_dynamics::ShellAcousticTimeStepSize computing_time_step_size(cylinder_body);
	/** active-pative stress relaxation. */
	thin_structure_dynamics::ShellStressRelaxationFirstHalf
		stress_relaxation_first_half(cylinder_body_inner);
	thin_structure_dynamics::ShellStressRelaxationSecondHalf
		stress_relaxation_second_half(cylinder_body_inner);
	solid_dynamics::ConstrainSolidBodyRegionVelocity
		constrain_holder(cylinder_body, new BoundaryGeometry(cylinder_body, "BoundaryGeometry"), Vecd(0.0, 1.0, 0.0));
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vecd>>
		cylinder_position_damping(cylinder_body_inner, 0.1, "Velocity", physical_viscosity);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vecd>>
		cylinder_rotation_damping(cylinder_body_inner, 0.1, "AngularVelocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		write_cylinder_max_displacement("Position", in_output, cylinder_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	cylinder_initial_pseudo_normal.parallel_exec();
	corrected_configuration_in_strong_form.parallel_exec();

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_cylinder_max_displacement.writeToFile(0);
	
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
		Real integeral_time = 0.0;
		while (integeral_time < output_period)
		{
			if (ite % 100 == 0) {
				std::cout << "N=" << ite << " Time: "
					<< GlobalStaticVariables::physical_time_ << "	dt: "
					<< dt << "\n";
				write_states.writeToFile(100);
			}
			dt = 0.3 * computing_time_step_size.parallel_exec();
			initialize_external_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);

			constrain_holder.parallel_exec();
			cylinder_position_damping.parallel_exec(dt);
			cylinder_rotation_damping.parallel_exec(dt);
			constrain_holder.parallel_exec();

			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			integeral_time += dt;
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

	return 0;
}

