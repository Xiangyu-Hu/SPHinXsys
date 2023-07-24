/**
* @file 	dw_3d_volum_oscillating_plate.cpp
* @brief 	This is the test case for the hourglass manuscript.
* @details  We consider vibration deformation of a square plate under initial vertical velocity field.
* @author 	Dong Wu, Chi Zhang and Xiangyu Hu
* @version  0.1
 */
#include "sphinxsys.h"
#include "all_continuum.h"
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 0.4;                                              /** Length of the square plate. */
Real PH = 0.4;                                              /** Width of the square plate. */
Real PT = 0.01;                                             /** Thickness of the square plate. */
int particle_number = 3;                                    /** Particle number in the direction of the thickness. */
Real particle_spacing_ref = PT / (Real)particle_number;     /** Initial reference particle spacing. */
int particle_nuber_PL = PL / particle_spacing_ref;
int particle_nuber_PH = PH / particle_spacing_ref;
int BWD = 1;
Real BW = particle_spacing_ref * (Real)BWD;                 /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-BW, -BW, -PT / 2), Vec3d(PL + BW, PH + BW, PT / 2));
// Observer location
StdVec<Vecd> observation_location = { Vecd(0.5 * PL, 0.5 * PH, 0.0),  Vecd(-BW, -BW, 0.0) };
/** For material properties of the solid. */
Real rho0_s = 1000.0; 			                               /** Normalized density. */
Real Youngs_modulus = 100.0e6;	                       /** Normalized Youngs Modulus. */
Real poisson = 0.3; 			                           /** Poisson ratio. */
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
Real gravity_g = 0.0;

Real governing_vibration_integer_x = 1.0;
Real governing_vibration_integer_y = 1.0;
Real U_max = 1.0;  //Maximum velocity
/** Define application dependent particle generator for thin structure. */
class PlateParticleGenerator : public ParticleGenerator
{
public:
	explicit PlateParticleGenerator(SPHBody& sph_body) : ParticleGenerator(sph_body) {};
	virtual void initializeGeometricVariables() override
	{
		for (int k = 0; k < particle_number; k++)
		{
			for (int i = 0; i < particle_nuber_PL + 2 * BWD; i++)
			{
				for (int j = 0; j < particle_nuber_PH + 2 * BWD; j++)
				{
					Real x = particle_spacing_ref * i - BW + particle_spacing_ref * 0.5;
					Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
					Real z = particle_spacing_ref * (k - ((particle_number - 1.0) / 2.0));
					initializePositionAndVolumetricMeasure(Vecd(x, y, z),
						particle_spacing_ref * particle_spacing_ref * particle_spacing_ref);
				}
			}
		}
	}
};

/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody& body, const std::string& body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry() {};

private:
	void tagManually(size_t index_i)
	{
		if ((base_particles_.pos_[index_i][2] < 0.5 * particle_spacing_ref)
			& (base_particles_.pos_[index_i][2] > -0.5 * particle_spacing_ref)
			& (base_particles_.pos_[index_i][0] < 0.0 || base_particles_.pos_[index_i][0] > PL
				|| base_particles_.pos_[index_i][1] < 0.0 || base_particles_.pos_[index_i][1] > PH))
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/** Define the initial condition. */
class BeamInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit BeamInitialCondition(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body) {};

	void update(size_t index_i, Real dt)
	{
		/** initial velocity profile */
		vel_[index_i][2] = sin(governing_vibration_integer_x * Pi * pos_[index_i][0] / PL)
			* sin(governing_vibration_integer_y * Pi * pos_[index_i][1] / PH);
	};
};

/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);

	/** create a plate body. */
	SolidBody plate_body(system, makeShared<DefaultShape>("PlateBody"));
	plate_body.defineParticlesAndMaterial<ContinuumParticles,  GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
	plate_body.generateParticles<PlateParticleGenerator>();
	//plate_body.addBodyStateForRecording<Real>("VolumetricStress");
	plate_body.addBodyStateForRecording<Real>("VonMisesStress");
	plate_body.addBodyStateForRecording<Real>("VonMisesStrain");
	plate_body.addBodyStateForRecording<Real>("Pressure");
	plate_body.addBodyStateForRecording<Real>("Density");

	/** Define Observer. */
	ObserverBody plate_observer(system, "PlateObserver");
	plate_observer.defineParticlesAndMaterial();
	plate_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	InnerRelation plate_body_inner(plate_body);
	ContactRelation plate_observer_contact(plate_observer, { &plate_body });
	/**
	 * This section define all numerical methods will be used in this case.
	 */
	SimpleDynamics<BeamInitialCondition> initial_velocity(plate_body);
	/** Time step size calculation. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(plate_body, U_max, 0.2);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(plate_body, 0.4);
	/** stress relaxation. */
	Dynamics1Level<continuum_dynamics::Integration1stHalf> plate_pressure_relaxation(plate_body_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfDissipativeRiemann> plate_density_relaxation(plate_body_inner);
	InteractionDynamics<continuum_dynamics::AngularConservativeShearAccelerationRelaxation>
		plate_shear_acceleration_angular_conservative(plate_body_inner);
	/** Corrected configuration. */
	InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(plate_body_inner);
	Dynamics1Level<continuum_dynamics::ShearStressRelaxation> plate_shear_stress_relaxation(plate_body_inner);
	/** Constrain the Boundary. */
	BoundaryGeometry boundary_geometry(plate_body, "BoundaryGeometry");
	SimpleDynamics<continuum_dynamics::FixedInAxisDirection> constrain_holder(boundary_geometry, Vecd(1.0, 1.0, 0.0));
	SimpleDynamics<continuum_dynamics::ConstrainSolidBodyMassCenter>
		constrain_mass_center(plate_body, Vecd(1.0, 1.0, 0.0));
	/** Output */
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RestartIO restart_io(io_environment, system.real_bodies_);
	// ObservedQuantityRecording<Vecd>
	// 	write_plate_displacement("Position", io_environment, plate_observer_contact);
	RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
		write_plate_displacement("Position", io_environment, plate_observer_contact);
	ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>
		write_kinetic_energy(io_environment, plate_body);

	/** Apply initial condition. */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	initial_velocity.exec();
	constrain_holder.exec();
	corrected_configuration.exec();

	/**
	* @brief The time stepping starts here.
	*/
	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
	}
	GlobalStaticVariables::physical_time_ = 0.0;
	/** first output. */
	write_states.writeToFile(0);
	write_plate_displacement.writeToFile(0);
	write_kinetic_energy.writeToFile(0);
	//write_elastic_strain_energy.writeToFile(0);

	/** Setup physical parameters. */
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 500;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 0.1;
	Real output_period = end_time / 100.0;
	/** Statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			Real relaxation_time = 0.0;
			Real advection_dt = fluid_advection_time_step.exec();

			while (relaxation_time < advection_dt)
			{
				Real acoustic_dt = fluid_acoustic_time_step.exec();
				plate_shear_stress_relaxation.exec(acoustic_dt);
				plate_pressure_relaxation.exec(acoustic_dt);
				constrain_holder.exec(acoustic_dt);
				plate_density_relaxation.exec(acoustic_dt);
				//shear acceleration with angular conservative
				plate_shear_acceleration_angular_conservative.exec(acoustic_dt);
				number_of_iterations++;
				relaxation_time += acoustic_dt;
				integration_time += acoustic_dt;
				GlobalStaticVariables::physical_time_ += acoustic_dt;
				if (number_of_iterations % screen_output_interval == 0)
				{
					std::cout << "N=" << number_of_iterations << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	advection_dt: "
						<< advection_dt << "	acoustic_dt: "
						<< acoustic_dt << "\n";
					if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.RestartStep())
						restart_io.writeToFile(Real(number_of_iterations));
				}
			}
			plate_body.updateCellLinkedList();
			plate_body_inner.updateConfiguration();
			corrected_configuration.exec(); //put correct matrix in the end
		}
		write_plate_displacement.writeToFile(number_of_iterations);
		write_kinetic_energy.writeToFile(number_of_iterations);
		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}

	TickCount t4 = TickCount::now();
	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	//system.generate_regression_data_ = true;
    if (system.generate_regression_data_)
    {
        write_plate_displacement.generateDataBase(Vecd(1.0e-2, 1.0e-2, 1.0e-2), Vecd(1.0e-2, 1.0e-2, 1.0e-2));
    }
    else
    {
        write_plate_displacement.testResult();
    }
	return 0;
}