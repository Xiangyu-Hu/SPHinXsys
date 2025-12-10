/**
 * @file dambreak_sycl.cpp
 * @brief 2D dambreak example using SYCL.
 * @author Xiangyu Hu
 */
#include "sphinxsys.h"

using namespace SPH; // Namespace cite here.
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
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Artificial sound speed. */
//----------------------------------------------------------------------
//	Parameters for the shapes.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with inital shape, materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    ComplexShape wall_complex_shape("WallBoundary");
    wall_complex_shape.add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
    wall_complex_shape.subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    SolidBody wall_boundary(sph_system, wall_complex_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The relations give the topological connections within a body
    //  or with other bodies within interaction range.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall_boundary});
    Contact<> fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &water_block_update_complex_relation = main_methods.addRelationDynamics(water_block_inner, water_wall_contact);
    auto &fluid_observer_contact_relation = main_methods.addRelationDynamics(fluid_observer_contact);
    auto &particle_sort = main_methods.addSortDynamics(water_block);

    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(water_block, gravity);
    auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary); // run on CPU
    auto &water_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);

    auto &fluid_linear_correction_matrix =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(water_block_inner, 0.5)
            .addPostContactInteraction(water_wall_contact);
    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact);
    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact);
    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::DensityRegularization, WithUpdate, FreeSurface, AllParticles>(water_block_inner)
            .addPostContactInteraction(water_wall_contact);

    auto &fluid_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_ref);
    auto &fluid_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_state_recorder.addToWrite<Real>(water_block, "Density");
    auto &restart_io = main_methods.addIODynamics<RestartIOCK>(sph_system);
    auto &record_water_mechanical_energy = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, TotalMechanicalEnergyCK>(water_block, gravity);
    auto &fluid_observer_pressure = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Real>("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(20.0);
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        time_stepper.setPhysicalTime(restart_io.readRestartFiles(sph_system.RestartStep()));
    }
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep() + 1;
    int screening_interval = 100;
    int observation_interval = screening_interval * 2;
    int restart_output_interval = screening_interval * 10;
    auto &state_recording = time_stepper.addTriggerByInterval(0.1);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec(); // run particle dynamics with host kernels first
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();

    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    fluid_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    record_water_mechanical_energy.writeToFile(advection_steps);
    fluid_observer_pressure.writeToFile(advection_steps);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            water_update_particle_position.exec();

            /** Output body state during the simulation according output_interval. */
            time_instance = TickCount::now();
            /** screen output, write body observables and restart files  */
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                record_water_mechanical_energy.writeToFile(advection_steps);
                fluid_observer_contact_relation.exec();
                fluid_observer_pressure.writeToFile(advection_steps);
            }

            if (advection_steps % restart_output_interval == 0)
            {
                restart_io.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (advection_steps % 100)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    //----------------------------------------------------------------------
    // Post-run regression test to ensure that the case is validated
    //----------------------------------------------------------------------
    if (sph_system.GenerateRegressionData())
    {
        record_water_mechanical_energy.generateDataBase(1.0e-3);
        fluid_observer_pressure.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        record_water_mechanical_energy.testResult();
        fluid_observer_pressure.testResult();
    }

    return 0;
};
