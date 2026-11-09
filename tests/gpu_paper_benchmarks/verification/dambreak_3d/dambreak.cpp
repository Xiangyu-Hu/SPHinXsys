/**
 * @file dambreak.cpp
 * @brief 3D dambreak example using computing kernels.
 * @author Xiangyu Hu
 */
#include "sphinxsys.h"
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include <algorithm>
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real global_resolution = 0.05;   // particle spacing
Real BW = global_resolution * 4; // boundary width
Real DL = 5.366;                 // tank length
Real DH = 2.0;                   // tank height
Real DW = 0.5;                   // tank width
Real LL = 2.0;                   // liquid length
Real LH = 1.0;                   // liquid height
Real LW = 0.5;                   // liquid width
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
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

class WallBoundary : public ComplexShape //	define the static solid wall boundary shape
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
//----------------------------------------------------------------------
//	A group of observation points.
//----------------------------------------------------------------------
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
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av);
    TickCount benchmark_start = TickCount::now();
    global_resolution = benchmark_config.dp;
    BW = global_resolution * 4;
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "dambreak_3d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    TimeInterval io_interval = TickCount::now() - environment_io_start;
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    WaterBlock initial_water_block("WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(createObservationPoints());
    const size_t initial_particle_count =
        water_block.getBaseParticles().TotalRealParticles() +
        wall_boundary.getBaseParticles().TotalRealParticles() +
        fluid_observer.getBaseParticles().TotalRealParticles();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall_boundary});
    Contact<> fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_block_update_complex_relation(water_block_inner, water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(fluid_observer_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);

    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(water_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary); // run on CPU
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_update_particle_position(water_block);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_block_inner, 0.5), water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::CompressionSummation<Inner<>, Contact<>>>
        fluid_density_summation(water_block_inner, water_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, WeaklyCompressibleFluid, FreeSurface>>
        fluid_density_regularization(water_block);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_block_inner, water_wall_contact);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_block, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>> fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "PositionDivergence");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>>
        record_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Real>>
        fluid_observer_pressure(fluid_observer_contact, "Pressure");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    wall_boundary_normal_direction.exec(); // run particle dynamics with host kernels first
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();
    const TimeInterval init_interval =
        TickCount::now() - initialization_start;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    size_t number_of_outer_steps = 0;
    size_t number_of_acoustic_steps = 0;
    int screen_output_interval = 100;
    Real end_time = Real(benchmark_config.end_time);
    Real output_interval = benchmark_config.output_interval > 0.0
                               ? Real(benchmark_config.output_interval)
                               : end_time / 20.0;
    const bool verification_mode = !benchmark_config.benchmark_mode;
    const bool write_heavy_output = benchmark_config.output_enabled;
    size_t feature_samples = 0;
    //----------------------------------------------------------------------
    //	Wall and I/O time. Kernel completion follows the synchronization
    //	semantics already required by reductions returning host-visible dt.
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    {
        TickCount initial_io_start = TickCount::now();
        record_water_mechanical_energy.writeToFile(number_of_iterations);
        fluid_observer_pressure.writeToFile(number_of_iterations);
        feature_samples++;
        if (write_heavy_output)
        {
            body_states_recording.writeToFile();
        }
        io_interval += TickCount::now() - initial_io_start;
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval &&
               sv_physical_time->getValue() < end_time)
        {
            fluid_density_summation.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            fluid_boundary_indicator.exec();
            fluid_linear_correction_matrix.exec();

            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt &&
                   integration_time < output_interval &&
                   sv_physical_time->getValue() < end_time)
            {

                acoustic_dt = std::min(
                    {fluid_acoustic_time_step.exec(),
                     advection_dt - relaxation_time,
                     output_interval - integration_time,
                     end_time - sv_physical_time->getValue()});
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
                number_of_acoustic_steps++;
            }
            water_update_particle_position.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                TickCount console_io_start = TickCount::now();
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";
                io_interval += TickCount::now() - console_io_start;
            }
            number_of_iterations++;

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            if (write_heavy_output)
            {
                TickCount observer_io_start = TickCount::now();
                fluid_observer_pressure.writeToFile(number_of_iterations);
                feature_samples++;
                io_interval += TickCount::now() - observer_io_start;
            }
        }

        {
            TickCount periodic_io_start = TickCount::now();
            record_water_mechanical_energy.writeToFile(number_of_iterations);
            fluid_observer_pressure.writeToFile(number_of_iterations);
            feature_samples += 2;
            if (write_heavy_output)
            {
                body_states_recording.writeToFile();
            }
            io_interval += TickCount::now() - periodic_io_start;
        }
        number_of_outer_steps++;
    }
    std::string reference_status = "not_checked_benchmark";
    if (verification_mode && write_heavy_output)
    {
        TickCount regression_io_start = TickCount::now();
        if (sph_system.GenerateRegressionData())
        {
            record_water_mechanical_energy.generateDataBase(1.0e-3);
            fluid_observer_pressure.generateDataBase(1.0e-3);
            reference_status = "generated";
        }
        else
        {
            record_water_mechanical_energy.testResult();
            fluid_observer_pressure.testResult();
            reference_status = "tested";
        }
        io_interval += TickCount::now() - regression_io_start;
    }
    else if (verification_mode)
    {
        reference_status = "not_checked_output_disabled";
    }

    const TimeInterval total_wall = TickCount::now() - benchmark_start;
    const double compute_seconds = std::max(
        0.0, total_wall.seconds() - init_interval.seconds() -
                 io_interval.seconds());
    std::cout << "Total wall time: " << total_wall.seconds()
              << " seconds (compute: " << compute_seconds
              << ", initialization: " << init_interval.seconds()
              << ", I/O: " << io_interval.seconds() << ")." << std::endl;

    paper_bench::BenchmarkSummary summary;
    summary.dp = global_resolution;
    summary.initial_particle_count = initial_particle_count;
    summary.particle_count =
        water_block.getBaseParticles().TotalRealParticles() +
        wall_boundary.getBaseParticles().TotalRealParticles() +
        fluid_observer.getBaseParticles().TotalRealParticles();
    summary.physical_time = sv_physical_time->getValue();
    summary.outer_steps = number_of_outer_steps;
    summary.advection_steps = number_of_iterations;
    summary.acoustic_steps = number_of_acoustic_steps;
    summary.wall_seconds = total_wall.seconds();
    summary.compute_seconds = compute_seconds;
    summary.io_seconds = io_interval.seconds();
    summary.init_seconds = init_interval.seconds();
    summary.time_per_outer_step =
        number_of_outer_steps == 0
            ? 0.0
            : compute_seconds / static_cast<double>(number_of_outer_steps);
    summary.setParticleUpdates(
        summary.particle_count * summary.advection_steps,
        "particle_count*advection_steps");
    summary.extra_fields.push_back(
        {"observer_quantities", "Pressure;TotalMechanicalEnergy"});
    summary.extra_fields.push_back(
        {"observer_sampled", feature_samples > 0 ? "true" : "false"});
    summary.extra_fields.push_back(
        {"feature_samples", std::to_string(feature_samples)});
    summary.extra_fields.push_back(
        {"reference_status", reference_status});
    benchmark_recorder.writeSummary(summary);

    return 0;
}
