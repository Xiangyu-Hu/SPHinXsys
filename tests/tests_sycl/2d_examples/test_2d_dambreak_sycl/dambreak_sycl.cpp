/**
 * @file dambreak_sycl.cpp
 * @brief 2D dambreak example using SYCL.
 * @author Xiangyu Hu
 */
#include "sphinxsys_sycl.h"
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
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
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
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    TransformShape<GeometricShapeBox> initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>(water_block_shape);

    WallBoundary wall_boundary_shape("WallBoundary");
    SolidBody wall_boundary(sph_system, wall_boundary_shape.getName());
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>(wall_boundary_shape);

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<execution::ParallelDevicePolicy, CellLinkedList> water_cell_linked_list(water_block);
    UpdateCellLinkedList<execution::ParallelDevicePolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);

    Relation<Inner<>> water_block_inner(water_block);
    Relation<Contact<>> water_wall_contact(water_block, {&wall_boundary});
    Relation<Contact<>> fluid_observer_contact(fluid_observer, {&water_block});

    SequencedCombination<UpdateRelation<execution::ParallelDevicePolicy, BodyRelationUpdate<Inner<>, Contact<>>>>
        water_block_update_complex_relation(water_block_inner, water_wall_contact);
    UpdateRelation<execution::ParallelDevicePolicy, BodyRelationUpdate<Contact<>>>
        fluid_observer_contact_relation(fluid_observer_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxillary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    ParticleSortCK<execution::ParallelDevicePolicy, RadixSort> particle_sort(water_block);
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<execution::ParallelDevicePolicy, GravityForceCK<Gravity>> constant_gravity(water_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_block);

    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        fluid_density_regularization(water_block_inner, water_wall_contact);

    ReduceDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_block, U_ref);
    ReduceDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticTimeStepCK> fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    RestartIO restart_io(sph_system);

    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<execution::ParallelDevicePolicy, TotalMechanicalEnergyCK>>
        record_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<execution::ParallelDevicePolicy, Real>>
        fluid_observer_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        sv_physical_time->setValue(restart_io.readRestartFiles(sph_system.RestartStep()));
    }

    wall_boundary_normal_direction.exec();
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();

    water_block_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();
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
    TimeInterval interval_acoustic_steps;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(execution::ParallelDevicePolicy{});
    record_water_mechanical_energy.writeToFile(number_of_iterations);
    fluid_observer_pressure.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();

            fluid_density_regularization.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            water_advection_step_setup.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            interval_acoustic_steps += TickCount::now() - time_instance;

            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    record_water_mechanical_energy.writeToFile(number_of_iterations);
                    fluid_observer_pressure.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(execution::ParallelDevicePolicy{}, number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_advection_step_close.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile(execution::ParallelDevicePolicy{});
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_steps = "
              << interval_acoustic_steps.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

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
