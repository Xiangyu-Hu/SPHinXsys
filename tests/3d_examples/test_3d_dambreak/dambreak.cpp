/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry
Real resolution_ref = 0.05;   // particle spacing
Real BW = resolution_ref * 4; // boundary width
Real DL = 5.366;              // tank length
Real DH = 2.0;                // tank height
Real DW = 0.5;                // tank width
Real LL = 2.0;                // liquid length
Real LH = 1.0;                // liquid height
Real LW = 0.5;                // liquid width

// for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

//	define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_water(halfsize_water);
        add<GeometricShapeBox>(Transform(translation_water), halfsize_water);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_wall(halfsize_inner);
        add<GeometricShapeBox>(Transform(translation_wall), halfsize_outer);
        subtract<GeometricShapeBox>(Transform(translation_wall), halfsize_inner);
    }
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    observation_points.push_back(Vecd(DL, 0.01, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.1, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.2, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.24, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.252, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.266, 0.5 * DW));
    return observation_points;
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    WaterBlock initial_water_block("WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_water_block_states(sph_system);
    write_water_block_states.addToWrite<Vec3d>(wall_boundary, "NormalDirection");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 20.0;
    Real output_interval = end_time / 20.0;
    Real dt = 0.0; // default acoustic time step sizes
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_water_block_states.writeToFile(0);
    write_water_mechanical_energy.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();

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

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            write_recorded_water_pressure.writeToFile(number_of_iterations);
        }

        write_water_mechanical_energy.writeToFile(number_of_iterations);

        TickCount t2 = TickCount::now();
        write_water_block_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        write_recorded_water_pressure.generateDataBase(1.0e-3);
    }
    else
    {
        write_water_mechanical_energy.testResult();
        write_recorded_water_pressure.testResult();
    }

    return 0;
}
