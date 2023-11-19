/**
 * @file 	Lid_driven_cavity.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_Lid_driven_cavity example.
 * @author 	Zhentong Wang
 */
#include "general_eulerian_fluid_dynamics.hpp" // eulerian classes for fluid.
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 129.0; /**< Global reference resolution. */
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
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class WeaklyCompressibleFluidInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit WeaklyCompressibleFluidInitialCondition(SPHBody& sph_body)
		: FluidInitialCondition(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
	{
		particles_->registerVariable(mom_, "Momentum");
		particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
		particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
	};

	void update(size_t index_i, Real dt) {};

protected:
	StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
	StdLargeVec<Real>& rho_, & p_;
};
using namespace SPH;   //	Namespace cite here.
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
	EulerianFluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_body.defineBodyLevelSetShape();
	water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
		: water_body.generateParticles<ParticleGeneratorLattice>();
	//water_body.generateParticles<ParticleGeneratorLattice>();
	water_body.addBodyStateForRecording<Real>("Density");
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_block_complex(water_body, { &wall_boundary });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		/** body topology only for particle relaxation */
		InnerRelation cylinder_inner(water_body);
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
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner, true);
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
		std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/** Initial condition with momentum and energy field */
	SimpleDynamics<MovingWallInitialCondition>  solid_initial_condition(wall_boundary);
	SimpleDynamics<WeaklyCompressibleFluidInitialCondition> fluid_initial_condition(water_body);
	/** Initialize a time step. */
	SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
	InteractionWithUpdate<KernalGredientWithCorrectionComplex> kernel_gredient_update(water_block_complex);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_body_correct_inner);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	InteractionWithUpdate<ICEIntegration1stHalfHLLRiemannWithWall> pressure_relaxation(water_block_complex, 0.0);
	InteractionWithUpdate<ICEIntegration1stHalfHLLRiemannWithWall> density_and_energy_relaxation(water_block_complex, 0.0);
	//InteractionDynamics<eulerian_weakly_compressible_fluid_dynamics::VorticityInner> vorticity_inner(water_body_inner);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	//wall_boundary_normal_direction.parallel_exec();
	fluid_initial_condition.exec();
	solid_initial_condition.exec();
	body_states_recording.writeToFile();
	kernel_gredient_update.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 30.0; /**< End time. */
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
		TickCount t3 = TickCount::now();
        interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
