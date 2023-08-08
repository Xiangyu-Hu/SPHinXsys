/**
 * @file 	column_collapse.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for soil simulation.
 */
#include "sphinxsys.h" //SPHinXsys Library.
#include "all_continuum.h"
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
//unit system - 1
Real DL = 0.5;					/**< Tank length. */
Real DH = 0.15;					/**< Tank height. */
Real LL = 0.2;						/**< Liquid column length. */
Real LH = 0.1;						/**< Liquid column height. */

Real particle_spacing_ref = LH / 50;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = { Vecd(DL, 0.2) };
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
/*
 * Dilatancy angle is always zero for non-associate flow rule
 */
Real rho0_s = 2650;						 /**< Reference density of soil. */
Real gravity_g = 9.8;					 /**< Gravity force of soil. */
Real Youngs_modulus = 0.84e6; //reference Youngs modulus
Real poisson = 0.3;		 //Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s*3*(1-2* poisson)));
Real friction_angle = 19.8 * Pi / 180;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d soil_block_translation = soil_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
		subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
	}
};
std::vector<Vecd> soil_shape{
		Vecd(0, 0), Vecd(0, LH), Vecd(LL, LH), Vecd(LL, 0), Vecd(0, 0) };

class Soil : public MultiPolygonShape
{
public:
	explicit Soil(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
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
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
	// handle command line arguments
	sph_system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	RealBody soil_block(
		sph_system, makeShared<Soil>("GranularBody"));
	soil_block.defineParticlesAndMaterial<PlasticContinuumParticles, PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle);
	soil_block.generateParticles<ParticleGeneratorLattice>();
	soil_block.addBodyStateForRecording<Real>("Pressure");
	soil_block.addBodyStateForRecording<Real>("Density");
	soil_block.addBodyStateForRecording<Real>("VerticalStress");
	soil_block.addBodyStateForRecording<Real>("AccDeviatoricPlasticStrain");

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation soil_block_complex(soil_block, { &wall_boundary });
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> soil_step_initialization(soil_block, gravity_ptr);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> soil_acoustic_time_step(soil_block, 0.2);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> soil_density_by_summation(soil_block_complex);
	InteractionDynamics<continuum_dynamics::StressDiffusion> stress_diffusion(soil_block_complex.getInnerRelation());
	//stress relaxation with Riemann solver
	Dynamics1Level<continuum_dynamics::StressRelaxation1stHalfRiemannWithWall> granular_stress_relaxation_1st(soil_block_complex);
	Dynamics1Level<continuum_dynamics::StressRelaxation2ndHalfRiemannWithWall> granular_stress_relaxation_2nd(soil_block_complex);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	IOEnvironment io_environment(sph_system);
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_mechanical_energy(io_environment, soil_block, gravity_ptr);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
		soil_block.updateCellLinkedList();
		soil_block_complex.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 500;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1.0; /**< End time. */
	Real D_Time = End_Time / 50;	  /**< Time stamps for output of body states. */
	Real Dt = 0.1*D_Time;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	TimeInterval interval_computing_time_step;
	TimeInterval interval_computing_soil_stress_relaxation;
	TimeInterval interval_updating_configuration;
	TickCount time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	write_mechanical_energy.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			time_instance = TickCount::now();
			soil_step_initialization.exec();

			soil_density_by_summation.exec();
			interval_computing_time_step += TickCount::now() - time_instance;

			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = soil_acoustic_time_step.exec();

				granular_stress_relaxation_1st.exec(dt);
				stress_diffusion.exec();
				granular_stress_relaxation_2nd.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				interval_computing_soil_stress_relaxation += TickCount::now() - time_instance;

				/** screen output, write body reduced values and restart files  */
				if (number_of_iterations % screen_output_interval == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4) << "	Time = "
						<< GlobalStaticVariables::physical_time_
						<< std::scientific << "	dt = " << dt << "\n";

					if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
					{
						write_mechanical_energy.writeToFile(number_of_iterations);
					}
					if (number_of_iterations % restart_output_interval == 0)
						restart_io.writeToFile(number_of_iterations);
				}
				number_of_iterations++;
				/** Update cell linked list and configuration. */
				soil_block.updateCellLinkedList();
				soil_block_complex.updateConfiguration();
			}
			time_instance = TickCount::now();

			interval_updating_configuration += TickCount::now() - time_instance;
		}
		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << std::fixed << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_soil_stress_relaxation = "
		<< interval_computing_soil_stress_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";

	//sph_system.generate_regression_data_ = true;
	if (sph_system.generate_regression_data_)
	{
		write_mechanical_energy.generateDataBase(1.0e-3);
	}
	else if (sph_system.RestartStep() == 0)
	{
		write_mechanical_energy.testResult();
	}

	return 0;
};