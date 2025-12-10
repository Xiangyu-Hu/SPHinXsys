/**
 * @file 	stlw_sycl.cpp
 * @brief 	This is the case file for 2D still water.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 10.0; /**< TOTAL SIMULATION TIME*/
Real DL = 3.0;                   /**< Tank length. */
Real DH = 4.0;                   /**< Tank height. */
Real WH = 2.0;                   /**< Water block height. */
Real particle_spacing_ref = 0.05;
Real BW = particle_spacing_ref * 4.0; /**< Extending width for BCs. */
BoundingBoxd system_domain_bounds(Vec2d(-DL - BW, -DH - BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                    /**< Reference density of fluid. */
Real gravity_g = 9.81;                   /**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                   /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_body_halfsize = Vec2d(0.5 * DL, 0.5 * WH);
Vec2d water_body_translation = Vec2d(0.0, -0.5 * WH);
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(0.0, 0.0);
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = Vec2d(0.0, 0.0);
//----------------------------------------------------------------------
//	Geometry of the measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vec2d gauge_halfsize = Vec2d(0.5 * h, 0.5 * DH);
Vec2d gauge_translation = Vec2d(DL / 3, 0.5 * DH);
//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox water_body_shape(Transform(water_body_translation), water_body_halfsize, "WaterBody");
    FluidBody water_body(sph_system, water_body_shape);
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    GeometricShapeBox wave_probe_buffer_shape(Transform(gauge_translation), gauge_halfsize, "FreeSurfaceGauge");
    BodyRegionByCell wave_probe_buffer(water_body, wave_probe_buffer_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_body_inner(water_body);
    Contact<> water_wall_contact(water_body, {&wall});
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
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_body_update_complex_relation(water_body_inner, water_wall_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_body);

    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(water_body, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal_direction(wall); // run on CPU
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_body);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_update_particle_position(water_body);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        fluid_density_regularization(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_body_inner, water_wall_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> fluid_acoustic_time_step(water_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_real_body_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MainExecutionPolicy, UpperFrontInAxisDirectionCK<BodyRegionByCell>>>
        wave_gauge(wave_probe_buffer, "FreeSurfaceHeight");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");

    wall_normal_direction.exec(); // run particle dynamics on CPU first
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_body_update_complex_relation.exec();
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    int number_of_iterations = 0;
    int screen_output_interval = 1000;
    int observation_interval = screen_output_interval / 20;
    Real end_time = total_physical_time;
    Real output_interval = end_time / 100;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    wave_gauge.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            fluid_viscous_force.exec();

            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                total_time += acoustic_dt;
                if (total_time >= relax_time)
                    sv_physical_time->incrementValue(acoustic_dt);
            }
            water_update_particle_position.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << "	Total Time = " << total_time
                          << "	Physical Time = " << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";
            }
            number_of_iterations++;

            if (number_of_iterations % 100 == 0)
            {
                particle_sort.exec();
            }

            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();

            if (total_time >= relax_time && number_of_iterations % observation_interval == 0)
            {
                wave_gauge.writeToFile(number_of_iterations);
            }
        }

        TickCount t2 = TickCount::now();
        if (total_time >= relax_time)
            write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        wave_gauge.generateDataBase(0.1);
    }
    else
    {
        wave_gauge.testResult();
    }

    return 0;
}
