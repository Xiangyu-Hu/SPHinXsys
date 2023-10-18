/**
 * @file	dambreak.cpp
 * @brief	2D dambreak example.
 * @details	This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
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
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(
        sph_system, makeShared<TransformShape<GeometricShapeBox>>(
                        Transform(water_block_translation), water_block_halfsize, "WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.getBaseParticles().registerDeviceMemory();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
    wall_boundary.getBaseParticles().registerDeviceMemory();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    fluid_observer.getBaseParticles().registerDeviceMemory();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    water_block_complex.getInnerRelation().allocateInnerConfigurationDevice();
    water_block_complex.getContactRelation().allocateContactConfiguration();
    fluid_observer_contact.allocateContactConfiguration();

    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall, ParallelSYCLDevicePolicy> fluid_pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall, ParallelSYCLDevicePolicy> fluid_density_relaxation(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex, ParallelSYCLDevicePolicy> fluid_density_by_summation(water_block_complex);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SharedPtr<Gravity> gravity_ptr = makeSharedDevice<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization, ParallelSYCLDevicePolicy> fluid_step_initialization(water_block, gravity_ptr);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize, ParallelSYCLDevicePolicy> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize, ParallelSYCLDevicePolicy> fluid_acoustic_time_step(water_block);

    water_block.getBaseParticles().copyToDeviceMemory();
    executionQueue.setWorkGroupSize(16);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RestartIO restart_io(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy, ParallelSYCLDevicePolicy>>>
        write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real, ParallelSYCLDevicePolicy>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists(execution::par_sycl).wait();
    auto system_configuration_update_event = sph_system.initializeSystemDeviceConfigurations();

    wall_boundary_normal_direction.exec();
    wall_boundary.getBaseParticles().copyToDeviceMemory();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList(execution::par_sycl).wait();
        system_configuration_update_event.wait();
        system_configuration_update_event.add(water_block_complex.updateDeviceConfiguration())
            .add(fluid_observer_contact.updateDeviceConfiguration());
    }

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 20.0;
    Real output_interval = 0.1;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TimeInterval interval_writing_files;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    auto async_body_states_copy_event = body_states_recording.copyDeviceData();
    async_body_states_copy_event.then([&] { body_states_recording.writeToFile(); }).wait();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    write_recorded_water_pressure.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_step_initialization.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            /** Wait on system configuration to finish updating */
            time_instance = TickCount::now();
            system_configuration_update_event.wait();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** Evaluation of density by summation */
            time_instance = TickCount::now();
            fluid_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_pressure_relaxation.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                GlobalStaticVariables::physical_time_ += acoustic_dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedListWithParticleSort(100, execution::par_sycl)
                .then([&, number_of_iterations=number_of_iterations, advection_dt=advection_dt, acoustic_dt=acoustic_dt]{
                          /** screen output, write body reduced values and restart files  */
                          if (number_of_iterations % screen_output_interval == 0)
                          {
                              std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                                        << GlobalStaticVariables::physical_time_
                                        << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                              if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                              {
                                  time_instance = TickCount::now();
                                  write_water_mechanical_energy.writeToFile(number_of_iterations);
                                  write_recorded_water_pressure.writeToFile(number_of_iterations);
                                  interval_writing_files += TickCount::now() - time_instance;
                              }
                              if (number_of_iterations % restart_output_interval == 0)
                              {
                                  time_instance = TickCount::now();
                                  restart_io.writeToFile(number_of_iterations);
                                  interval_writing_files += TickCount::now() - time_instance;
                              }
                          }
                }).wait();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** Submit task to copy data from device to host in preparation of next file output */
            if (integration_time >= output_interval)
                async_body_states_copy_event = body_states_recording.copyDeviceData();

            /** Submit task for system configuration update */
            system_configuration_update_event = water_block_complex.updateDeviceConfiguration().add(
                fluid_observer_contact.updateDeviceConfiguration());

            number_of_iterations++;
        }

        time_instance = TickCount::now();
        async_body_states_copy_event.then([&] { body_states_recording.writeToFile(); }).wait();
        interval_writing_files += TickCount::now() - time_instance;
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    executionQueue.getQueue().wait_and_throw();

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_writing_files = "
              << interval_writing_files.seconds() << "\n";

    if (sph_system.generate_regression_data_)
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        write_recorded_water_pressure.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_water_mechanical_energy.testResult();
        write_recorded_water_pressure.testResult();
    }

    return 0;
};
