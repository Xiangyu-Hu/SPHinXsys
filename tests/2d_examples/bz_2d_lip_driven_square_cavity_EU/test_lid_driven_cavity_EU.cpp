/**
 * @file 	Lid_driven_cavity.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_Lid_driven_cavity example.
 * @author 	Zhentong Wang
 */
#include "eulerian_fluid_dynamics.hpp" // eulerian classes for fluid.
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 50.0; /**< Global reference resolution. */
Real BW = resolution_ref * 4;	 /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 1000.0;						/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	/** Geometry definition. */
	std::vector<Vecd> water_body_shape;
	water_body_shape.push_back(Vecd(0.0, 0.0));
	water_body_shape.push_back(Vecd(0.0, DH));
	water_body_shape.push_back(Vecd(DL, DH));
	water_body_shape.push_back(Vecd(DL, 0.0));
	water_body_shape.push_back(Vecd(0.0, 0.0));

	return water_body_shape;
}
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon outer_boundary(createWaterBlockShape());
		add<MultiPolygonShape>(outer_boundary, "InnerWaterBlock");
	}
};

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};

//----------------------------------------------------------------------
//	Application dependent initial condition
//----------------------------------------------------------------------
class MovingWallInitialCondition
	: public LocalDynamics, public SolidDataSimple
{
public:
	explicit MovingWallInitialCondition(SolidBody& solid_body)
		: LocalDynamics(solid_body), SolidDataSimple(solid_body),
		Vol_(particles_->Vol_), vel_(particles_->vel_), n_(particles_->n_), pos_(particles_->pos_) {};

	void update(size_t index_i, Real dt)
	{
		if (pos_[index_i][1] > DH)
		{
			vel_[index_i][0] = 1.0;
			vel_[index_i][1] = 0.0;
		}
	}
protected:
	StdLargeVec<Vecd>& vel_, & n_, & pos_;
	StdLargeVec<Real>& Vol_;
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class VelocityXObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit VelocityXObserverParticleGenerator(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_point = DL / resolution_ref;
		Real range_of_measure = DL - resolution_ref;
		Real start_of_measure = 0.5 * resolution_ref;

		for (size_t i = 0; i < number_of_observation_point; ++i)
		{
			Vec2d point_corrdinate(range_of_measure * (Real)i / (Real)(number_of_observation_point - 1) + start_of_measure, 0.5 * DL);
			positions_.push_back(point_corrdinate);
		}
	}
};
class VelocityYObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit VelocityYObserverParticleGenerator(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_point = DL / resolution_ref;
		Real range_of_measure = DL - resolution_ref;
		Real start_of_measure = 0.5 * resolution_ref;

		for (size_t i = 0; i < number_of_observation_point; ++i)
		{
			Vec2d point_corrdinate(0.5 * DH, range_of_measure * (Real)i / (Real)(number_of_observation_point - 1) + start_of_measure);
			positions_.push_back(point_corrdinate);
		}
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	sph_system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	sph_system.setReloadParticles(true);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_body.defineComponentLevelSetShape("InnerWaterBlock")->writeLevelSet(io_environment);
	water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
		: water_body.generateParticles<ParticleGeneratorLattice>();
	water_body.addBodyStateForRecording<Real>("Density");

	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
	//----------------------------------------------------------------------
	//	Particle and body creation of fluid observers.
	//----------------------------------------------------------------------
	ObserverBody horizontal_observer(sph_system, "HorizontalVelocity");
	horizontal_observer.generateParticles<VelocityXObserverParticleGenerator>();
	ObserverBody vertical_observer(sph_system, "VerticalVelocity");
	vertical_observer.generateParticles<VelocityYObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_block_complex(water_body, { &wall_boundary });
	ComplexRelation water_block_complex_corrected(water_body, { &wall_boundary });
	ContactRelation horizontal_observer_contact(horizontal_observer, { &water_body });
	ContactRelation vertical_observer_contact(vertical_observer, { &water_body });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		/** body topology only for particle relaxation */
		InnerRelation water_body_inner(water_body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(water_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_water_body_to_vtp(io_environment, { &water_body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, water_body);

		/* Relaxation method: including based on the 0th and 1st order consistency. */
		InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex(water_block_complex);

		//relax_dynamics::RelaxationStepComplexImplicit relaxation_step_complex(water_block_complex, "InnerWaterBlock");
		relax_dynamics::RelaxationStepComplexImplicit<CorrectionMatrixRelaxation> relaxation_step_complex(water_block_complex, "InnerWaterBlock");
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_water_block_kinetic_energy(water_body_inner);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_water_block_average_kinetic_energy(water_body, "ParticleKineticEnergy");
		ReduceDynamics<QuantityMaximum<Real>> calculate_water_block_maximum_kinetic_energy(water_body, "ParticleKineticEnergy");
		water_body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		relaxation_step_complex.SurfaceBounding().exec();
		write_water_body_to_vtp.writeToFile(0);
		Real water_block_average_kinetic_energy = 100.0;
		Real water_block_maximum_kinetic_energy = 100.0;
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.
		GlobalStaticVariables::physical_time_ = ite;
		while (ite < 200)
		{
			kernel_correction_complex.exec();
			relaxation_step_complex.exec();

			ite += 1;
			if (ite % 200 == 0)
			{
				update_water_block_kinetic_energy.exec();
				water_block_average_kinetic_energy = calculate_water_block_average_kinetic_energy.exec();
				water_block_maximum_kinetic_energy = calculate_water_block_maximum_kinetic_energy.exec();
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				std::cout << "Body: "
					      << " Average: " << water_block_average_kinetic_energy
					      << " Maximum: " << water_block_maximum_kinetic_energy << std::endl;
				write_water_body_to_vtp.writeToFile(ite);
			}
		}

		std::cout << "Residual: "
			      << " maximum: " << calculate_water_block_maximum_kinetic_energy.exec()
			      << " average: " << calculate_water_block_average_kinetic_energy.exec() << std::endl;

		std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

		ite++;
		write_water_body_to_vtp.writeToFile(ite);
		write_particle_reload_files.writeToFile(0);

		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initial condition with momentum and energy field */
	SimpleDynamics<MovingWallInitialCondition>  solid_initial_condition(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> boundary_normal_direction(wall_boundary);
	/** Initialize a time step. */
	SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_maxtrix(water_block_complex);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_maxtrix_corrected(water_block_complex_corrected);
	InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(kernel_correction_maxtrix_corrected);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_corrected);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfAcousticRiemannWithWall> pressure_relaxation(water_block_complex_corrected);
	InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfAcousticRiemannWithWall> density_and_energy_relaxation(water_block_complex_corrected);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	ObservedQuantityRecording<Vecd> write_horizontal_velocity("Velocity", io_environment, horizontal_observer_contact);
	ObservedQuantityRecording<Vecd> write_vertical_velocity("Velocity", io_environment, vertical_observer_contact);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	solid_initial_condition.exec();
	boundary_normal_direction.exec();
	body_states_recording.writeToFile();
	kernel_correction_maxtrix.exec();
	kernel_correction_maxtrix_corrected.exec();
	kernel_gradient_update.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 50.0; /**< End time. */
	Real D_Time = 1.0;	 /**< Time stamps for output of body states. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
    TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force. */
			time_step_initialization.exec();
			Real dt = get_fluid_time_step_size.exec();
			viscous_acceleration.exec();
			/** Dynamics including pressure relaxation. */
			integration_time += dt;
			pressure_relaxation.exec(dt);
			density_and_energy_relaxation.exec(dt);
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
				{
					restart_io.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;
		}

		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		write_horizontal_velocity.writeToFile(number_of_iterations);
		write_vertical_velocity.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
        interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
