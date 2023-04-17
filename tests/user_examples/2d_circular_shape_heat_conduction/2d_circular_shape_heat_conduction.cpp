/**
 * @file 	2d_circular_shape_heat_conduction.cpp
 * @brief 	2D circular shape model of heat conduction problem.
 * @details This is a case to implement three kinds of boundary conditions related to heat conduction problem.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
#include "2d_circular_shape_heat_conduction.h"


using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec2d insert_circle_center(0.0, 0.0);
Real insert_inner_wall_circle_radius = 0.65;
Real insert_in_circle_radius = 0.75;
Real insert_out_circle_radius = 1.0;
Real insert_outer_wall_circle_radius = 1.1;
Real resolution_ref = 0.02;
BoundingBox system_domain_bounds(Vec2d(-1.1, -1.1), Vec2d(1.1, 1.1));
// observer location
StdVec<Vecd> observation_location = { Vecd(0.8, 0.5) };

//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1;
Real bias_coff = 0.0;
std::array<std::string, 1> species_name_list{ "Phi" };

Real rho0_f = 1.0;					/**< Density. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Speed of sound. */
Real Re = 100.0;					/**< Reynolds number100. */
Real mu_f = rho0_f * U_f * 2 * (insert_out_circle_radius - insert_in_circle_radius) / Re; /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Global parameters on the initial condition.
//----------------------------------------------------------------------
//Real phi_outer_wall = 40.0;
Real phi_inner_wall = 350.0;
Real phi_fluid_initial = 200.0;
Real T_infinity = 100.0;
Real convection = 100.0;

//----------------------------------------------------------------------
//	Definition of the fluid block shape.
//----------------------------------------------------------------------
class FluidBlock : public ComplexShape
{
public:
	explicit FluidBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_out_circle_radius, 100, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_in_circle_radius, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};

//----------------------------------------------------------------------
//	Definition of the cylinder.
//----------------------------------------------------------------------
class InnerWall : public ComplexShape
{
public:
	explicit InnerWall(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_in_circle_radius, 100, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_inner_wall_circle_radius, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Definition of the outer wall.
//----------------------------------------------------------------------
class OuterWall : public ComplexShape
{
public:
	explicit OuterWall(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_outer_wall_circle_radius, 100, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_out_circle_radius, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};
////----------------------------------------------------------------------
////	Definition of the wall.
////----------------------------------------------------------------------
//class SolidWall : public ComplexShape
//{
//public:
//	explicit SolidWall(const std::string& shape_name) : ComplexShape(shape_name)
//	{
//		MultiPolygon multi_polygon;
//		multi_polygon.addACircle(insert_circle_center, insert_outer_wall_circle_radius, 100, ShapeBooleanOps::add);
//		multi_polygon.addACircle(insert_circle_center, insert_out_circle_radius, 100, ShapeBooleanOps::sub);
//		multi_polygon.addACircle(insert_circle_center, insert_in_circle_radius, 100, ShapeBooleanOps::add);
//		multi_polygon.addACircle(insert_circle_center, insert_inner_wall_circle_radius, 100, ShapeBooleanOps::sub);
//		add<MultiPolygonShape>(multi_polygon);
//	}
//};

//----------------------------------------------------------------------
//	Setup fluid material properties.
//----------------------------------------------------------------------
class FluidMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	FluidMaterial() : DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, SharedPtr<NoReaction>(), rho0_f, c_f, mu_f)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	};
};
using DiffusionFluidParticles = DiffusionReactionParticles<FluidParticles, FluidMaterial>;
//----------------------------------------------------------------------
//	Setup solid wall material properties.
//----------------------------------------------------------------------
class WallMaterial : public DiffusionReaction<Solid>
{
public:
	WallMaterial() : DiffusionReaction<Solid>({ "Phi" }, SharedPtr<NoReaction>())
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
	};
};
using DiffusionSolidParticles = DiffusionReactionParticles<SolidParticles, WallMaterial>;

//----------------------------------------------------------------------
//	Solid wall initial condition.
//----------------------------------------------------------------------
class WallInitialCondition
	: public DiffusionReactionInitialCondition<DiffusionSolidParticles>
{
protected:
	size_t phi_;

public:
	explicit WallInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<DiffusionSolidParticles>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		if (pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2) <= pow(insert_in_circle_radius, 2)
			&& pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2) >= pow(insert_inner_wall_circle_radius, 2))
		{

			all_species_[phi_][index_i] = phi_inner_wall;
		}

		if (pow(insert_out_circle_radius, 2) <= pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2)
			&& pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2) <= pow(insert_outer_wall_circle_radius, 2)
			&& pos_[index_i][0] < 0.1 * pos_[index_i][1] && pos_[index_i][0] > 0)
		{
			convection_[index_i] = convection;
			T_infinity_[index_i] = T_infinity;
		}

	};
};

//----------------------------------------------------------------------
//	Fluid initial condition.
//----------------------------------------------------------------------
class FluidInitialCondition
	: public DiffusionReactionInitialCondition<DiffusionFluidParticles>
{
protected:
	size_t phi_;

public:
	explicit FluidInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<DiffusionFluidParticles>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		if (pow(insert_in_circle_radius, 2) <= pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2)
			&& pow(pos_[index_i][0], 2) + pow(pos_[index_i][1], 2) <= pow(insert_out_circle_radius, 2))
		{
			all_species_[phi_][index_i] = phi_fluid_initial;
		}
	};
};
//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
class ThermalRelaxationWithBoundary
	: public RelaxationOfAllDiffusionSpeciesRK2<
	RelaxationOfAllDiffusionSpeciesWithBoundary<DiffusionFluidParticles, DiffusionSolidParticles>>
{
public:
	explicit ThermalRelaxationWithBoundary(ComplexRelation& complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(complex_relation) {};
	virtual ~ThermalRelaxationWithBoundary() {};
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.setRunParticleRelaxation(false);
	/** Tag for computation start with relaxed body fitted particles distribution. */
	sph_system.setReloadParticles(true);
	//handle command line arguments
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av);
#endif

	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody fluid_body(sph_system, makeShared<FluidBlock>("FluidBlock"));
	fluid_body.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, FluidMaterial>();
	fluid_body.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_body(sph_system, makeShared<SolidWall>("SolidWall"));
	wall_body.defineBodyLevelSetShape();
	wall_body.defineParticlesAndMaterial<DiffusionReactionParticles<SolidParticles, Solid>, WallMaterial>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? wall_body.generateParticles<ParticleGeneratorReload>(io_environment, wall_body.getName())
		: wall_body.generateParticles<ParticleGeneratorLattice>();

	//----------------------------------------------------------------------
	//	Particle and body creation of fluid observers.
	//----------------------------------------------------------------------
	ObserverBody temperature_observer(sph_system, "FluidObserver");
	temperature_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation fluid_body_inner(fluid_body);
	//InnerRelation wall_body_inner(wall_body);
	//ComplexRelation fluid_body_complex(fluid_body_inner, { &wall_body });
	ComplexRelation fluid_body_complex(fluid_body, { &wall_body });
	ComplexRelation wall_boundary_complex(wall_body, { &fluid_body });
	//ContactRelation wall_contact(wall_body, { &fluid_body });
	ContactRelation temperature_observer_contact(temperature_observer, { &fluid_body });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		/** body topology only for particle relaxation */
		InnerRelation wall_body_inner(wall_body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, { &wall_body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, wall_body);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_body_inner);
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
		std::cout << "The physics relaxation process of the wall finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}

	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<FluidInitialCondition> fluid_initial_condition(fluid_body);
	SimpleDynamics<WallInitialCondition> wall_condition(wall_body);
	//SimpleDynamics<NormalDirectionFromBodyShape> wall_body_normal_direction(wall_body);
	/** Initialize particle acceleration. */
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(fluid_body);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_body_complex);

	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step(fluid_body, U_f);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step(fluid_body);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<DiffusionFluidParticles> get_thermal_time_step(fluid_body);

	/** Diffusion process between two diffusion bodies. */
	ThermalRelaxationWithBoundary thermal_relaxation_with_boundary(fluid_body_complex);

	InteractionDynamics<UpdateUnitVectorNormalToBoundary<FluidParticles, WeaklyCompressibleFluid, SolidParticles, Solid>> update_fluid_body_normal_vector(fluid_body_complex);
	InteractionDynamics<UpdateUnitVectorNormalToBoundary<SolidParticles, Solid, FluidParticles, WeaklyCompressibleFluid>> update_wall_boundary_normal_vector(wall_boundary_complex);

	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(fluid_body_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(fluid_body_complex);
	/** Computing viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_body_complex);
	/** Apply transport velocity formulation. */
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction(fluid_body_complex);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_body_inner);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
	//RegressionTestEnsembleAveraged<ObservedQuantityRecording<Real>> write_fluid_phi("Phi", io_environment, temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

	/** computing surface normal direction for the wall. */
	wall_condition.exec();
	fluid_initial_condition.exec();
	update_fluid_body_normal_vector.exec();
	update_wall_boundary_normal_vector.exec();

	Real dt_thermal = get_thermal_time_step.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	Real end_time = 10;
	Real output_interval = end_time / 100.0; /**< time stamps for output,WriteToFile*/
	int number_of_iterations = 0;
	int screen_output_interval = 40;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.exec();
			Real Dt = get_fluid_advection_time_step.exec();
			update_density_by_summation.exec();
			viscous_acceleration.exec();
			transport_velocity_correction.exec();

			size_t inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(SMIN(dt_thermal, get_fluid_time_step.exec()), Dt);
				pressure_relaxation.exec(dt);
				density_relaxation.exec(dt);
				thermal_relaxation_with_boundary.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
			}
			number_of_iterations++;

			/** Water block configuration. */
			fluid_body.updateCellLinkedListWithParticleSort(100);
			fluid_body_complex.updateConfiguration();
		}
		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
		compute_vorticity.exec();
		temperature_observer_contact.updateConfiguration();
		write_states.writeToFile();
		//write_fluid_phi.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}

	TickCount t4 = TickCount::now();
	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	//write_fluid_phi.newResultTest();

	return 0;
}
