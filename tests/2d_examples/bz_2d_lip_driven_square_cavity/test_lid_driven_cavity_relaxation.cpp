/**
 * @file 	Lid_driven_square_cavity.cpp
 * @brief 	2d lip driven square cavity example
 * @details This is the one of the basic test cases.
 * @author 	Bo Zhang, Zhentong Wang
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 25.0; /**< Global reference resolution. */
Real BW = resolution_ref * 6;	  /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 1000.0;					/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_body_shape;
		water_body_shape.push_back(Vecd(0.0, 0.0));
		water_body_shape.push_back(Vecd(0.0, DH));
		water_body_shape.push_back(Vecd(DL, DH));
		water_body_shape.push_back(Vecd(DL, 0.0));
		water_body_shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);

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
		  vel_(particles_->vel_), pos_(particles_->pos_) {};

	void update(size_t index_i, Real dt)
	{
		if (pos_[index_i][1] > DH)
		{
			vel_[index_i][0] = 1.0;
			vel_[index_i][1] = 0.0;
		}
	}
protected:
	StdLargeVec<Vecd>& vel_, & pos_;
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
		size_t number_of_observation_point = 51;
		Real range_of_measure = 1.0;
		Real start_of_measure = 0.0;

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
		size_t number_of_observation_point = 51;
		Real range_of_measure = 1.0;
		Real start_of_measure = 0.0;

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
	sph_system.setReloadParticles(false);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_body.defineBodyLevelSetShape();
	water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
		: water_body.generateParticles<ParticleGeneratorLattice>();
	water_body.addBodyStateForRecording<Real>("Density");
	water_body.addBodyStateForRecording<Real>("Pressure");

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
	ComplexRelation wall_boundary_complex(wall_boundary, { &water_body });
	ComplexRelation water_block_complex(water_body, { &wall_boundary });
	ContactRelation horizontal_observer_contact(horizontal_observer, { &water_body });
	ContactRelation vertical_observer_contact(vertical_observer, { &water_body });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		/** body topology only for particle relaxation */
		InnerRelation water_block_inner(water_body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(water_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, { &water_body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, water_body);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(water_block_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		write_inserted_body_to_vtp.writeToFile(0);

		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_inserted_body_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of the water block finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/** Initialize a time step. */
	SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
	/** Initial condition with momentum field */
	SimpleDynamics<MovingWallInitialCondition>  solid_initial_condition(wall_boundary);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_complex);
	/** Pressure and density relaxation algorithm by using verlet time stepping. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannCorrectWithWall> pressure_relaxation_with_correction(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannConsistencyWithWall> pressure_relaxation_with_consistency(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
	/** Kernel correction matrix and transport velocity formulation. */
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex(water_block_complex);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex_wall(wall_boundary_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(water_block_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityConsistencyComplex<AllParticles>> transport_velocity_consistency(water_block_complex, 0.25);
	InteractionSplit<fluid_dynamics::TransportVelocityConsistencyComplexImplicit<AllParticles>> transport_velocity_consistency_implicit(water_block_complex, 10);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_body, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	/** Computing viscous acceleration with wall. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
	water_body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToPlt body_states_recording(io_environment, sph_system.real_bodies_);
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
	body_states_recording.writeToFile();
	kernel_correction_complex.exec();
	kernel_correction_complex_wall.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	Real End_Time = 100.0; /**< End time. */
	Real output_interval = 0.1;
	Real dt = 1.0;	      /**< Time stamps for output of body states. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
    TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			time_step_initialization.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();
			update_density_by_summation.exec();
			viscous_acceleration.exec();
			
			kernel_correction_complex.exec();
			transport_velocity_correction.exec();
			//transport_velocity_consistency.exec();
			//transport_velocity_consistency_implicit.exec();
			
			Real relaxation_time = 0.0;
			while(relaxation_time < Dt)
			{
				// avoid possible smaller acoustic time step size for viscous flow
				dt = SMIN(get_fluid_time_step_size.exec(), Dt);
				relaxation_time += dt;
				integration_time += dt;
				//pressure_relaxation.exec(dt);
				//pressure_relaxation_with_correction.exec(dt);
				pressure_relaxation_with_consistency.exec(dt);
				density_relaxation.exec(dt);
				GlobalStaticVariables::physical_time_ += dt;
			}
			
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	dt = " << dt << "\n";
			}
			number_of_iterations++;
			water_body.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			horizontal_observer_contact.updateConfiguration();
			vertical_observer_contact.updateConfiguration();
		}
		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		if (GlobalStaticVariables::physical_time_ > 90)
		{
			write_horizontal_velocity.writeToFile(number_of_iterations);
			write_vertical_velocity.writeToFile(number_of_iterations);
		}
		
		TickCount t3 = TickCount::now();
        interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
