/**
 * @file 	Lid_driven_square_cavity.cpp
 * @brief 	2d lip driven square cavity example
 * @details This is the one of the basic test cases for fluid dynamics.
 * @author 	Bo Zhang, Xiangyu Hu
 */
#include "sphinxsys_ck.h" //	SPHinXsys Library.
using namespace SPH;      //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
// geometry data
Real height = 1;
Real width = 1;
Real boundary_width = particle_spacing * 4; // boundary width
Real resolution_ref = height / 50.0;        /**< Global reference resolution. */
BoundingBox system_domain_bounds(Vecd(-boundary_width * 2, -boundary_width * 2),
                                 Vecd(width + boundary_width * 2, height + boundary_width * 2));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                      /**< Reference density of fluid. */
Real U_f = 1.0;                         /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                  /**< Reference sound speed. */
Real Re = 100.0;                        /**< Reynolds number. */
Real mu_f = rho0_f * U_f * height / Re; /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class LidBoundary : public ComplexShape
{
  public:
    explicit LidBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width + boundary_width, 0.5 * boundary_width);
        Transform translate_to_origin(scaled_container);
        Vecd transform(-boundary_width, height);
        Transform translate_to_position(transform + scaled_container);
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_position), scaled_container);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container_outer(0.5 * width + boundary_width, 0.5 * height + boundary_width);
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin_outer(Vec2d(-boundary_width, -boundary_width) + scaled_container_outer);
        Transform translate_to_origin_inner(scaled_container);

        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin_outer), scaled_container_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin_inner), scaled_container);
    }
};
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(scaled_container);
        add<TransformShape<GeometricShapeBox>>(Transform(translate_to_origin), scaled_container);
    }
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
        Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_of_observation_point - 1) + start_of_measure, 0.5 * height);
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
        Vec2d point_coordinate(0.5 * width, range_of_measure * (Real)i / (Real)(number_of_observation_point - 1) + start_of_measure);
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
    IOEnvironment io_environment(sph_system);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
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
    using MyExecutionPolicy = execution::ParallelDevicePolicy; // define execution policy for this case

    UpdateCellLinkedList<MyExecutionPolicy, CellLinkedList> water_cell_linked_list(water_body);
    UpdateCellLinkedList<MyExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);

    Relation<Inner<>> water_body_inner(water_body);
    Relation<Contact<>> water_wall_contact(water_body, {&wall_boundary});
    Relation<Contact<>> horizontal_observer_contact(horizontal_observer, {&water_body});
    Relation<Contact<>> vertical_observer_contact(vertical_observer, {&water_body});

    UpdateRelation<MyExecutionPolicy, Inner<>, Contact<>> water_body_update_complex_relation(water_body_inner, water_wall_contact);
    UpdateRelation<MyExecutionPolicy, Contact<>> horizontal_observer_contact_relation(horizontal_observer_contact);
    UpdateRelation<MyExecutionPolicy, Contact<>> vertical_observer_contact_relation(vertical_observer_contact);
    ParticleSortCK<MyExecutionPolicy, QuickSort> particle_sort(water_body);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<MyExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(water_body, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<MyExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_body);
    StateDynamics<MyExecutionPolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_body);
    BodyRegionByParticle lid_boundary(wall_boundary, makeShared<LidBoundary>("LidBoundary"));
    StateDynamics<MyExecutionPolicy, ConstantConstraintCK<BodyRegionByParticle, Vec2d>> lid_velocity(lid_boundary, "Velocity", Vec2d(U_f, 0.0));

    InteractionDynamicsCK<MyExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(ConstructorArgs(water_body_inner, 0.5), water_wall_contact);
    InteractionDynamicsCK<MyExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MyExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MyExecutionPolicy, fluid_dynamics::DensityRegularizationComplex>
        fluid_density_regularization(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MyExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        viscous_acceleration(water_body_inner, water_wall_contact);

    ReduceDynamicsCK<MyExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_ref);
    ReduceDynamicsCK<MyExecutionPolicy, fluid_dynamics::AcousticTimeStepCK> fluid_acoustic_time_step(water_body);

    InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplex<AllParticles>>
        transport_velocity_correction(water_body_inner, water_body_contact);
    /** Computing viscous acceleration with wall. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_body_inner, water_body_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vecd>(water_body, "Velocity");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_horizontal_velocity("Velocity", horizontal_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_vertical_velocity("Velocity", vertical_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    solid_initial_condition.exec();
    write_real_body_states.writeToFile();
    kernel_correction_complex.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real End_Time = 30.0; /**< End time. */
    Real output_interval = 1.0;
    Real dt = 1.0; /**< Time stamps for output of body states. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    /** Output the start states of bodies. */
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();

            kernel_correction_complex.exec();
            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                // avoid possible smaller acoustic time step size for viscous flow
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                relaxation_time += dt;
                integration_time += dt;
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            water_body.updateCellLinkedList();
            water_body_complex.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        horizontal_observer_contact.updateConfiguration();
        vertical_observer_contact.updateConfiguration();
        write_horizontal_velocity.writeToFile(number_of_iterations);
        write_vertical_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

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
