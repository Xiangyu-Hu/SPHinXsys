/**
 * @file 	Lid_driven_square_cavity.cpp
 * @brief 	2d lip driven square cavity example
 * @details This is the one of the basic test cases for the RKGC inner flow.
 * @author 	Bo Zhang, Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;                    /**< box length. */
Real DH = 1.0;                    /**< box height. */
Real resolution_ref = 1.0 / 50.0; /**< Global reference resolution. */
Real BW = resolution_ref * 6;     /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                  /**< Reference density of fluid. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Reference sound speed. */
Real Re = 100.0;                    /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
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
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
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
class BoundaryVelocity : public MotionConstraint<SPHBody>
{
  public:
    BoundaryVelocity(SPHBody &body)
        : MotionConstraint<SPHBody>(body) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        if (pos_[index_i][1] > DH)
        {
            vel_[index_i][0] = 1.0;
            vel_[index_i][1] = 0.0;
        }
    };
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
StdVec<Vecd> VelocityXObserverParticle()
{
    StdVec<Vecd> observation_points;
    size_t number_of_observation_point = 5;
    Real range_of_measure = 1.0 - 0.5 * resolution_ref;
    Real start_of_measure = 0.5 * resolution_ref;

    for (size_t i = 0; i < number_of_observation_point; ++i)
    {
        Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_of_observation_point - 1) + start_of_measure, 0.5 * DL);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
}

StdVec<Vecd> VelocityYObserverParticle()
{
    StdVec<Vecd> observation_points;
    size_t number_of_observation_point = 5;
    Real range_of_measure = 1.0 - 0.5 * resolution_ref;
    Real start_of_measure = 0.5 * resolution_ref;
    for (size_t i = 0; i < number_of_observation_point; ++i)
    {
        Vec2d point_coordinate(0.5 * DH, range_of_measure * (Real)i /
                                                 (Real)(number_of_observation_point - 1) +
                                             start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
}
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
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of fluid observers.
    //----------------------------------------------------------------------
    ObserverBody horizontal_observer(sph_system, "HorizontalVelocity");
    horizontal_observer.generateParticles<ObserverParticles>(VelocityXObserverParticle());
    ObserverBody vertical_observer(sph_system, "VerticalVelocity");
    vertical_observer.generateParticles<ObserverParticles>(VelocityYObserverParticle());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_body);
    Contact<> water_wall_contact(water_body, {&wall_boundary});
    Contact<> horizontal_observer_contact(horizontal_observer, {&water_body});
    Contact<> vertical_observer_contact(vertical_observer, {&water_body});
    // ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
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
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_block_update_complex_relation(water_block_inner, water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> horizontal_observer_contact_relation(horizontal_observer_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> vertical_observer_contact_relation(vertical_observer_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_body);
    //----------------------------------------------------------------------

    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary); // run on CPU
    /** Time step size with considering sound wave speed. */
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_body);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_update_particle_position(water_body);

    /** Initial condition with momentum field */
    SimpleDynamics<BoundaryVelocity> solid_initial_condition(wall_boundary);
    /** Kernel correction matrix and transport velocity formulation. */
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_block_inner, 0.5), water_wall_contact);
    /** Evaluation of density by summation approach. */
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplex>
        fluid_density_regularization(water_block_inner, water_wall_contact); /** Pressure and density relaxation algorithm by using Verlet time stepping. */
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_wall_contact);
    /**
     * Free Surface Indicator and BulkParticles for Transport Velocity Correction are not required for this simulation.
     * They are included here solely for testing and verification purposes,
     * and do not contribute to the primary objectives of the current case.
     */
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCK>
        transport_correction_ck(water_block_inner, water_wall_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> fluid_acoustic_time_step(water_body);
    /** Computing viscous acceleration with wall. */
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(water_body, "Velocity");
    body_states_recording.addToWrite<Real>(water_body, "Density");
    body_states_recording.addToWrite<int>(water_body, "Indicator");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "Velocity");

    RestartIOCK<MainExecutionPolicy> restart_io(sph_system);

    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd>> write_horizontal_velocity("Velocity", horizontal_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd>> write_vertical_velocity("Velocity", vertical_observer_contact);
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
    wall_boundary_normal_direction.exec(); // run particle dynamics on CPU first
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    horizontal_observer_contact_relation.exec();
    vertical_observer_contact_relation.exec();
    solid_initial_condition.exec();
    fluid_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 30.0;
    Real output_interval = 1.0;
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval_writing_body_state;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_acoustic_steps;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_horizontal_velocity.writeToFile(number_of_iterations);
    write_vertical_velocity.writeToFile(number_of_iterations);
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

            water_advection_step_setup.exec();
            fluid_viscous_force.exec();
            fluid_linear_correction_matrix.exec();
            fluid_boundary_indicator.exec();
            transport_correction_ck.exec();

            Real advection_dt = fluid_advection_time_step.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = SMIN(fluid_acoustic_time_step.exec(), advection_dt);
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            water_update_particle_position.exec();
            interval_acoustic_steps += TickCount::now() - time_instance;

            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Particle sort, ipdate cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        TickCount t2 = TickCount::now();
        /** Output body state during the simulation according output_interval. */
        horizontal_observer_contact_relation.exec();
        vertical_observer_contact_relation.exec();
        write_horizontal_velocity.writeToFile(number_of_iterations);
        write_vertical_velocity.writeToFile(number_of_iterations);
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval_writing_body_state += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval_writing_body_state;
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
        write_horizontal_velocity.generateDataBase(1.0e-3);
        write_vertical_velocity.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_horizontal_velocity.testResult();
        write_vertical_velocity.testResult();
    }

    return 0;
}
