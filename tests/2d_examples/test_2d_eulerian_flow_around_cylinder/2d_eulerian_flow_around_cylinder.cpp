/**
 * @file 	2d_eulerian_flow_around_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;						   /**< Channel length. */
Real DH = 10.0;						   /**< Channel height. */
Real resolution_ref = 1.0 / 5.0;	   /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Vec2d cylinder_center(4.0, 5.0);	   /**< Location of the cylinder center. */
Real cylinder_radius = 1.0;			   /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;										 /**< Density. */
Real U_f = 1.0;											 /**< freestream velocity. */
Real c_f = 10.0 * U_f;									 /**< Speed of sound. */
Real Re = 100.0;										 /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	return water_block_shape;
}
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon outer_boundary(createWaterBlockShape());
		add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
		MultiPolygon circle(cylinder_center, cylinder_radius, 100);
		subtract<MultiPolygonShape>(circle);
	}
};
class Cylinder : public MultiPolygonShape
{
public:
	explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Far field boundary condition
//----------------------------------------------------------------------
class FarFieldBoundary : public eulerian_weakly_compressible_fluid_dynamics::NonReflectiveBoundaryVariableCorrection
{
public:
	explicit FarFieldBoundary(BaseInnerRelation &inner_relation)
		: eulerian_weakly_compressible_fluid_dynamics::NonReflectiveBoundaryVariableCorrection(inner_relation)
	{
		rho_farfield_ = rho0_f;
		vel_farfield_ = Vecd(U_f, 0.0);
		sound_speed_ = c_f;
	};
	virtual ~FarFieldBoundary(){};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	sph_system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	sph_system.setReloadParticles(true);
	// Handle command line arguments and override the tags for particle relaxation and reload.
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	EulerianFluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
	water_block.defineComponentLevelSetShape("OuterBoundary");
	water_block.defineParticlesAndMaterial<WeaklyCompressibleFluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<int>("SurfaceIndicator");

	SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
	cylinder.defineAdaptationRatios(1.15, 2.0);
	cylinder.defineBodyLevelSetShape();
	cylinder.defineParticlesAndMaterial<SolidParticles, Solid>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
		: cylinder.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//	Note that the same relation should be defined only once.
	//----------------------------------------------------------------------
	ComplexRelation water_block_complex(water_block, {&cylinder});
	ContactRelation cylinder_contact(cylinder, {&water_block});
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		InnerRelation cylinder_inner(cylinder); // extra body topology only for particle relaxation
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
		SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
		BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
		ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner, true);
		relax_dynamics::RelaxationStepComplex relaxation_step_complex(water_block_complex, "OuterBoundary", true);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.exec(0.25);
		random_water_body_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		relaxation_step_complex.SurfaceBounding().exec();
		write_real_body_states.writeToFile(0);

		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			relaxation_step_complex.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
				write_real_body_states.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process finish !" << std::endl;

		write_real_body_particle_reload_files.writeToFile(0);

		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<eulerian_weakly_compressible_fluid_dynamics::EulerianFlowTimeStepInitialization> initialize_a_fluid_step(water_block);
	SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
	Dynamics1Level<eulerian_weakly_compressible_fluid_dynamics::Integration1stHalfHLLCRiemannWithLimiterWithWall> pressure_relaxation(water_block_complex);
	InteractionWithUpdate<eulerian_weakly_compressible_fluid_dynamics::Integration2ndHalfHLLCRiemannWithLimiterWithWall> density_relaxation(water_block_complex);
	InteractionDynamics<eulerian_weakly_compressible_fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
	ReduceDynamics<eulerian_weakly_compressible_fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationComplex> surface_indicator(water_block_complex.getInnerRelation(), water_block_complex.getContactRelation());
	InteractionDynamics<FarFieldBoundary> variable_reset_in_boundary_condition(water_block_complex.getInnerRelation());
	//----------------------------------------------------------------------
	//	Compute the force exerted on solid body due to fluid pressure and viscosity
	//----------------------------------------------------------------------
	InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(cylinder_contact);
	InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_solid_update(cylinder_contact, viscous_force_on_solid);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
	RegressionTestTimeAveraged<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
		write_total_viscous_force_on_inserted_body(io_environment, viscous_force_on_solid, "TotalViscousForceOnSolid");
	ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>
		write_total_force_on_inserted_body(io_environment, fluid_force_on_solid_update, "TotalPressureForceOnSolid");
	ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>> write_maximum_speed(io_environment, water_block);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	cylinder_normal_direction.exec();
	surface_indicator.exec();
	variable_reset_in_boundary_condition.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 1000;
	Real end_time = 80.0;
	Real output_interval = 5.0; /**< time stamps for output. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.exec();
			Real dt = get_fluid_time_step_size.exec();
			viscous_acceleration.exec();
			pressure_relaxation.exec(dt);
			density_relaxation.exec(dt);

			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			variable_reset_in_boundary_condition.exec();

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					 << GlobalStaticVariables::physical_time_
					 << "	dt = " << dt << "\n";
			}
			number_of_iterations++;
		}

		TickCount t2 = TickCount::now();
		write_real_body_states.writeToFile();

		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		write_total_force_on_inserted_body.writeToFile(number_of_iterations);

		write_maximum_speed.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	if (sph_system.generate_regression_data_)
	{
		// The lift force at the cylinder is very small and not important in this case.
		write_total_viscous_force_on_inserted_body.generateDataBase(Vec2d{1.0e-2, 1.0e-2}, Vec2d{1.0e-2, 1.0e-2});
	}
	else
	{
		write_total_viscous_force_on_inserted_body.newResultTest();
	}

	return 0;
}
