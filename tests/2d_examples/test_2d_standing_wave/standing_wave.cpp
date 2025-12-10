/**
 * @file	standing_wave.cpp
 * @brief	2D standing_wave example.
 * @author	Bo Zhang, Yaru Ren, Chi Zhang and Xiangyu Hu
 * The correction method is RKGC, more details referring to arXiv:2406.0257.
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;                      /**< Tank length. */
Real DH = 2.0;                      /**< Tank height. */
Real LL = 2.0;                      /**< Liquid column length. */
Real LH = 1.0;                      /**< Liquid column height. */
Real particle_spacing_ref = 0.02;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                    /**< Reference density of fluid. */
Real gravity_g = 9.81;                   /**< Gravity force of fluid. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
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
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//	water block
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    int Nh = 100;
    Real Lstep = DL / Nh;

    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    std::vector<Vecd> pnts;
    for (int n = 0; n <= Nh; n++)
    {
        Real x = n * Lstep;
        Real y = 0.1 * cos(PI * x);

        pnts.push_back(Vecd(x, y));
    }
    for (int n = 0; n <= Nh - 0; n++)
    {
        water_block_shape.push_back(Vecd(pnts[n][0], pnts[n][1] + 1.0));
    }
    water_block_shape.push_back(Vecd(LL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));

    return water_block_shape;
}

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
    }
};
//----------------------------------------------------------------------
//	wave gauge
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(1.0 - h, 0.0));
    pnts.push_back(Vecd(1.0 - h, DH));
    pnts.push_back(Vecd(1.0 + h, DH));
    pnts.push_back(Vecd(1.0 + h, 0.0));
    pnts.push_back(Vecd(1.0 - h, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    // Using relaxed particle distribution if needed
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();
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
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(DynamicsArgs(water_block_inner, 0.5), water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> fluid_pressure_relaxation_correct(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> fluid_advection_time_step(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_block, gravity);
    /** WaveProbes. */
    BodyRegionByCell wave_probe_buffer_(water_block, makeShared<MultiPolygonShape>(createWaveProbeShape(), "WaveProbe"));
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>>
        wave_probe(wave_probe_buffer_, "FreeSurfaceHeight");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    if (sph_system.RestartStep() != 0)
    {
        physical_time = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_block_complex.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 10.0;
    Real output_interval = 0.05;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    wave_probe.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real advection_dt = 0.3 * fluid_advection_time_step.exec();
            fluid_density_by_summation.exec();
            corrected_configuration_fluid.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_pressure_relaxation_correct.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                physical_time += acoustic_dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    wave_probe.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile();
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
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        wave_probe.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_water_mechanical_energy.testResult();
        wave_probe.testResult();
    }

    return 0;
};
