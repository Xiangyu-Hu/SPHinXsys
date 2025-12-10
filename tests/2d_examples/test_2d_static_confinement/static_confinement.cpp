/**
 * @file 	static_confinement.cpp
 * @brief 	2D dambreak example in which the solid wall boundary are static confinement.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;              /**< Tank length. */
Real DH = 5.366;              /**< Tank height. */
Real LL = 2.0;                /**< Liquid column length. */
Real LH = 1.0;                /**< Liquid column height. */
Real resolution_ref = 0.025;  /**< Global reference resolution. */
Real BW = resolution_ref * 4; /**< Extending width for BCs. */
// Observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity force of fluid. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, LH));
    water_block_shape.push_back(Vecd(LL, LH));
    water_block_shape.push_back(Vecd(LL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));
    return water_block_shape;
}
/** create wall shape */
std::vector<Vecd> createWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, DH));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, 0.0));

    return inner_wall_shape;
}
/** create a structure shape */
std::vector<Vecd> createStructureShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
    water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH + 0.5 * LH));
    water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH));
    water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
    return water_block_shape;
}
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Shape for the wall.
//----------------------------------------------------------------------
class WallShape : public MultiPolygonShape
{
  public:
    explicit WallShape(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Shape for a structure.
//----------------------------------------------------------------------
class Triangle : public MultiPolygonShape
{
  public:
    explicit Triangle(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createStructureShape(), ShapeBooleanOps::add);
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
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);

    Dynamics1Level<fluid_dynamics::Integration1stHalfInnerRiemann> pressure_relaxation(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> density_relaxation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    /** Define the confinement condition for wall. */
    NearShapeSurface near_surface_wall(water_block, makeShared<WallShape>("Wall"));
    near_surface_wall.getLevelSetShape().writeLevelSet(sph_system);
    fluid_dynamics::StaticConfinement confinement_condition_wall(near_surface_wall);
    /** Define the confinement condition for structure. */
    NearShapeSurface near_surface_triangle(water_block, makeShared<InverseShape<Triangle>>("Triangle"));
    near_surface_triangle.getLevelSetShape().writeLevelSet(sph_system);
    fluid_dynamics::StaticConfinement confinement_condition_triangle(near_surface_triangle);
    /** Push back the static confinement condition to corresponding dynamics. */
    update_density_by_summation.post_processes_.push_back(&confinement_condition_wall.density_summation_);
    update_density_by_summation.post_processes_.push_back(&confinement_condition_triangle.density_relaxation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_wall.pressure_relaxation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition_triangle.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_triangle.density_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition_wall.surface_bounding_);
    density_relaxation.post_processes_.push_back(&confinement_condition_triangle.surface_bounding_);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 20.0;       /**< End time. */
    Real output_interval = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;              /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_water_mechanical_energy.writeToFile(0);
    write_recorded_water_pressure.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    write_recorded_water_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_inner.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
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
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
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
}
