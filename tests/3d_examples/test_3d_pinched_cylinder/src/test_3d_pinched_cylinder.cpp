/**
* @file 	dw_3d_pinched_cylinder.cpp
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
Real radius = 4.953;                                     /** Radius of the inner wall of the cylinder. */
Real height = 10.35;                                     /** Height of the cylinder. */
Real thickness = 0.094;                                  /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;      /** Radius of the mid surface. */
int particle_number = 80;								 /** Particle number in the peripheral direction. */
Real particle_spacing_ref = 2 * radius_mid_surface * Pi / (Real)particle_number;     /** Initial reference particle spacing. */
int particle_number_height = height / particle_spacing_ref;
int BWD = 0;
Real BW = particle_spacing_ref * (Real)BWD;              /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, -BW, -radius - thickness),
	Vec3d(radius + thickness, height + BW, radius + thickness));

// Observer location
StdVec<Vecd> observation_location = { Vecd(0.0, 0.5 * height, radius_mid_surface),
									  Vecd(radius_mid_surface, 0.5 * height, 0.0),
									  Vecd(radius_mid_surface, 0.0, 0.0) };

/** For material properties of the solid. */
Real rho0_s = 1.0; 			                             /** Normalized density. */
Real Youngs_modulus = 1.05e7;	                         /** Normalized Youngs Modulus. */
Real poisson = 0.3125; 			                         /** Poisson ratio. */
Real physical_viscosity = 2.0e2;                         /** physical damping, here we choose the same value as numerical viscosity. */

/** Difine point forces. */
Real F_full = 50.0e3;
std::vector<Vecd> point_forces{ Vec3d(0.0, 0.0, F_full), Vec3d(0.0, 0.0, -F_full) };
std::vector<Vecd> reference_positions{ Vec3d(0.0, height / 2.0, radius_mid_surface),
							   Vec3d(0.0, height / 2.0, -radius_mid_surface) };
Real time_to_full_external_force = 0.5;

Real gravitational_acceleration = 0.0;

/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit CylinderParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body) {};
	virtual void initializeGeometricVariables() override
	{
		// the cylinder and boundary
		for (int i = 0; i < particle_number; i++)
		{
			for (int j = 0; j < (height / particle_spacing_ref + 2 * BWD - 1); j++)
			{
				Real x = radius_mid_surface * cos(i * 2 * Pi / (Real)particle_number);
				Real y = particle_spacing_ref * j - BW;
				Real z = radius_mid_surface * sin(i * 2 * Pi / (Real)particle_number);
				initializePositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
				Vecd n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
				initializeSurfaceProperties(n_0, thickness);
			}
		}
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
	SolidBody cylinder_body(system, makeShared<DefaultShape>("CylinderBody"));
	cylinder_body.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	cylinder_body.generateParticles<CylinderParticleGenerator>();

	cylinder_body.addBodyStateForRecording<Vec3d>("PriorAcceleration");
	cylinder_body.addBodyStateForRecording<Vec3d>("NormalDirection");
	cylinder_body.addBodyStateForRecording<Vec3d>("InitialNormalDirection");
	cylinder_body.addBodyStateForRecording<Mat3d>("CorrectionMatrix");
	cylinder_body.addBodyStateForRecording<Mat3d>("TransformationMatrix");
	cylinder_body.addBodyStateForRecording<Mat3d>("DeformationGradient");
	cylinder_body.addBodyStateForRecording<Real>("Density");

	/** Define Observer. */
	ObserverBody cylinder_observer(system, "CylinderObserver");
	cylinder_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation cylinder_body_inner(cylinder_body);
	ContactRelation cylinder_observer_contact(cylinder_observer, { &cylinder_body });

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	 /** Corrected configuration. */
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
		corrected_configuration(cylinder_body_inner);
	/** Time step size calculation. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(cylinder_body);
	/** stress relaxation. */
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf>
		stress_relaxation_first_half(cylinder_body_inner, 3, true);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf>
		stress_relaxation_second_half(cylinder_body_inner);
	SimpleDynamics<thin_structure_dynamics::DistributingPointForcesToShell>
		apply_point_forces(cylinder_body, point_forces, reference_positions,
			time_to_full_external_force, particle_spacing_ref, 1.15);
	SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> 
		constrain_mass_center(cylinder_body, Vecd(1.0, 1.0, 1.0));
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_position_damping(0.1, cylinder_body_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
		cylinder_rotation_damping(0.1, cylinder_body_inner, "AngularVelocity", physical_viscosity);
	/** Output */
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_cylinder_max_displacement("Position", io_environment, cylinder_observer_contact);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();

	/**
	* From here the time stepping begines.
	* Set the starting time.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	write_states.writeToFile(0);
	write_cylinder_max_displacement.writeToFile(0);

	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 2.0;
	Real output_period = end_time / 200.0;
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
			}
			apply_point_forces.parallel_exec(dt);

			stress_relaxation_first_half.parallel_exec(dt);
			constrain_mass_center.parallel_exec(dt);
			cylinder_position_damping.parallel_exec(dt);
			cylinder_rotation_damping.parallel_exec(dt);
			constrain_mass_center.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = 0.2 * computing_time_step_size.parallel_exec();
			integeral_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			Real check_time = GlobalStaticVariables::physical_time_;
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

