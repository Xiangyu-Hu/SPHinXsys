/**
 * @file twisting_column.cpp
 * @brief This is an example of solid with classic Neohookean model
 * to demonstrate the robustness of the formulation with Kirchhoff stress decomposition.
 * @author Chi Zhang and Xiangyu Hu
 * @ref DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include <algorithm>
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 6.0;                  /**< X-direction size. */
Real PH = 1.0;                  /**< Y-direction size. */
Real PW = 1.0;                  /**< Z-direction size. */
Real particle_spacing_ref = PH / 10.0;
StdVec<Vecd> observation_location = {Vecd(PL, 0.0, 0.0)};
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 1100.0; /**< Reference density. */
Real poisson = 0.45;  /**< Poisson ratio. */
Real Youngs_modulus = 1.7e7;
Real angular_0 = -300.0;
//------------------------------------------------------------------------------
// define a velocity profile for initial condition
//------------------------------------------------------------------------------
class VelocityProfile : public ReturnFunction<Vecd>
{
    Real angular_0_ = angular_0;
    Real PL_ = PL;

  public:
    VelocityProfile() {};

    Vecd operator()(const Vec3d &position)
    {
        Real x = position[0];
        Real y = position[1];
        Real z = position[2];
        Real angular_velocity = angular_0_ * sin((M_PI * x) / (2.0 * PL_));
        Real local_radius = sqrt(pow(y, 2) + pow(z, 2));
        Real angular = atan2(y, z);

        if (x > 0.0)
        {
            return Vecd(0.0, angular_velocity * local_radius * cos(angular),
                        -angular_velocity * local_radius * sin(angular));
        }
        else
        {
            return Vecd::Zero();
        }
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    const paper_bench::BenchmarkDefaults benchmark_defaults{
        PH / 10.0, 0.5, 0.5 / 250.0,
        {{"standard", PH / 10.0},
         {"coarse", PH / 5.0},
         {"fine", PH / 20.0}}};
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av, benchmark_defaults);
    TickCount benchmark_start = TickCount::now();
    particle_spacing_ref = Real(benchmark_config.dp);

    TimeInterval interval_output;
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "twisting_column_3d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    interval_output += TickCount::now() - environment_io_start;
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    // Build up an SPHSystem and IO environment.
    // Please the make sure the global domain bounds are correctly defined.
    //----------------------------------------------------------------------
    const Real BW = particle_spacing_ref * 0.0;
    const Real SL = particle_spacing_ref;
    Vecd halfsize_column(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
    Vecd translation_column(0.5 * (PL - SL), 0.0, 0.0);
    Vecd halfsize_holder(
        0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
    Vecd translation_holder(-0.5 * (SL + BW), 0.0, 0.0);
    Vec3d domain_lower_bound(
        -SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
    Vec3d domain_upper_bound(
        PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
    //----------------------------------------------------------------------
    //	Setup geometry first.
    //----------------------------------------------------------------------
    auto &column_shape = sph_system.addShape<ComplexShape>("Column");
    column_shape.add<GeometricShapeBox>(Transform(translation_column), halfsize_column);
    auto &holder_shape = sph_system.addShape<GeometricShapeBox>(Transform(translation_holder), halfsize_holder, "Holder");
    column_shape.add<GeometricShapeBox>(holder_shape);
    //----------------------------------------------------------------------
    // Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    auto &column = sph_system.addBody<SolidBody>(column_shape);
    column.defineMatterMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    column.generateParticles<BaseParticles, Lattice>();
    BodyRegionByParticle holder(column, holder_shape);

    auto &my_observer = sph_system.addBody<ObserverBody>("MyObserver");
    my_observer.generateParticles<ObserverParticles>(observation_location);
    const size_t initial_particle_count =
        column.getBaseParticles().TotalRealParticles() +
        my_observer.getBaseParticles().TotalRealParticles();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    auto &column_inner = sph_system.addInnerRelation(column, ConfigType::Lagrangian);
    auto &my_observer_contact = sph_system.addContactRelation(my_observer, column, ConfigType::Lagrangian);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
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
    auto &host_methods = sph_solver.getHostMethodContainer();
    auto &main_methods = sph_solver.getMainMethodContainer();
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(column));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(column_inner));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(my_observer_contact));

    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);
    auto &column_acoustic_step_1st_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, NoKernelCorrectionCK>(column_inner);
    auto &column_acoustic_step_2nd_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(column_inner);

    auto &column_acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(column, 0.5);
    auto &constraint_holder = main_methods.addStateDynamics<FixBodyPartConstraintCK>(holder);
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<VelocityProfile>>(column, "Velocity").exec();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_displacement = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>(my_observer_contact, "Position");
    auto &write_velocity = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>(my_observer_contact, "Velocity");
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    size_t acoustic_steps = 1;
    int screening_interval = 100;
    int observation_interval = screening_interval;
    const Real end_time = Real(benchmark_config.end_time);
    const Real output_interval = Real(benchmark_config.output_interval);
    const bool verification_mode = !benchmark_config.benchmark_mode;
    const bool write_heavy_output = benchmark_config.output_enabled;
    size_t feature_samples = 0;
    auto &state_recording =
        time_stepper.addTriggerByInterval(output_interval);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_acoustic_step;
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    lagrangian_configuration.exec();
    column_linear_correction_matrix.exec();
    const TimeInterval init_interval =
        TickCount::now() - initialization_start;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    {
        TickCount initial_io_start = TickCount::now();
        write_displacement.writeToFile(acoustic_steps);
        write_velocity.writeToFile(acoustic_steps);
        feature_samples++;
        if (write_heavy_output)
        {
            body_state_recorder.writeToFile();
        }
        interval_output += TickCount::now() - initial_io_start;
    }
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    size_t executed_acoustic_steps = 0;
    while (!time_stepper.isEndTime(end_time))
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(column_acoustic_time_step);
        column_acoustic_step_1st_half.exec(acoustic_dt);
        constraint_holder.exec();
        column_acoustic_step_2nd_half.exec(acoustic_dt);
        executed_acoustic_steps++;
        interval_acoustic_step += TickCount::now() - time_instance;

        /** Output body state during the simulation according output_interval. */
        time_instance = TickCount::now();
        /** screen output, write body observables and restart files  */
        if (acoustic_steps == 1 || acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (acoustic_steps % observation_interval == 0)
        {
            write_displacement.writeToFile(acoustic_steps);
            write_velocity.writeToFile(acoustic_steps);
            feature_samples++;
        }

        if (write_heavy_output && state_recording())
        {
            body_state_recorder.writeToFile();
        }
        interval_output += TickCount::now() - time_instance;

        acoustic_steps++;
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    std::string reference_status = "not_checked_benchmark";
    if (verification_mode && write_heavy_output)
    {
        TickCount regression_io_start = TickCount::now();
        if (sph_system.GenerateRegressionData())
        {
            write_displacement.generateDataBase(0.005);
            write_velocity.generateDataBase(0.005);
            reference_status = "generated";
        }
        else
        {
            write_displacement.testResult();
            write_velocity.testResult();
            reference_status = "tested";
        }
        interval_output += TickCount::now() - regression_io_start;
    }
    else if (verification_mode)
    {
        reference_status = "not_checked_output_disabled";
    }

    const TimeInterval total_wall = TickCount::now() - benchmark_start;
    const double compute_seconds = std::max(
        0.0, total_wall.seconds() - init_interval.seconds() -
                 interval_output.seconds());
    std::cout << "Total wall time: " << total_wall.seconds()
              << " seconds (compute: " << compute_seconds
              << ", initialization: " << init_interval.seconds()
              << ", I/O: " << interval_output.seconds() << ")." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";

    paper_bench::BenchmarkSummary summary;
    summary.dp = particle_spacing_ref;
    summary.initial_particle_count = initial_particle_count;
    summary.particle_count =
        column.getBaseParticles().TotalRealParticles() +
        my_observer.getBaseParticles().TotalRealParticles();
    summary.physical_time = time_stepper.getPhysicalTime();
    summary.outer_steps = executed_acoustic_steps;
    summary.acoustic_steps = executed_acoustic_steps;
    summary.wall_seconds = total_wall.seconds();
    summary.compute_seconds = compute_seconds;
    summary.io_seconds = interval_output.seconds();
    summary.init_seconds = init_interval.seconds();
    summary.time_per_outer_step =
        executed_acoustic_steps == 0
            ? 0.0
            : compute_seconds /
                  static_cast<double>(executed_acoustic_steps);
    summary.setParticleUpdates(
        summary.particle_count * summary.acoustic_steps,
        "particle_count*acoustic_steps");
    summary.acoustic_component_seconds =
        paper_bench::formatMetric(interval_acoustic_step.seconds());
    summary.extra_fields.push_back(
        {"observer_quantities", "Position;Velocity"});
    summary.extra_fields.push_back(
        {"observer_sampled", feature_samples > 0 ? "true" : "false"});
    summary.extra_fields.push_back(
        {"feature_samples", std::to_string(feature_samples)});
    summary.extra_fields.push_back(
        {"reference_status", reference_status});
    benchmark_recorder.writeSummary(summary);

    return 0;
}
