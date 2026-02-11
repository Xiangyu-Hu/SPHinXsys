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
Real global_resolution = 1.0 / 50.0; /**< Global reference resolution. */
Real BW = global_resolution * 6;     /**< Extending width for BCs. */
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
    Real range_of_measure = 1.0 - 0.5 * global_resolution;
    Real start_of_measure = 0.5 * global_resolution;

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
    Real range_of_measure = 1.0 - 0.5 * global_resolution;
    Real start_of_measure = 0.5 * global_resolution;
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
    SPHSystem sph_system(system_domain_bounds, global_resolution);
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
    InnerRelation water_block_inner(water_body);
    ContactRelation water_block_contact(water_body, {&wall_boundary});
    ContactRelation horizontal_observer_contact(horizontal_observer, {&water_body});
    ContactRelation vertical_observer_contact(vertical_observer, {&water_body});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Initial condition with momentum field */
    SimpleDynamics<BoundaryVelocity> solid_initial_condition(wall_boundary);
    /** Kernel correction matrix and transport velocity formulation. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_block_contact);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    /** Pressure and density relaxation algorithm by using Verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplex<AllParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_body);
    /** Computing viscous acceleration with wall. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
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
            water_block_complex.updateConfiguration();
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
