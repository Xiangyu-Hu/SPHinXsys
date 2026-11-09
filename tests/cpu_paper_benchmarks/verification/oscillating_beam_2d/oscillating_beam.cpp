/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D oscillation beam example-one body version           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases, also the first case for            *
 * understanding SPH method for solid simulation.                              *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include <algorithm>
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;                  // beam length
Real PH = 0.02;                 // for thick plate; =0.01 for thin plate
Real SL = 0.06;                 // depth of the insert
// reference particle spacing
Real global_resolution = PH / 10.0;
Real BW = global_resolution * 4; // boundary width, at least three particles
//----------------------------------------------------------------------
//	Material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         // reference density
Real Youngs_modulus = 2.0e6; // reference Youngs modulus
Real poisson = 0.3975;       // Poisson ratio
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class LinearProfile : public ReturnFunction<Vecd>
{
    Real vf_ = vf;
    Real c0_;
    Real kl_ = kl;
    Real M_ = M;
    Real N_ = N;
    Real Q_ = Q;
    Real PL_ = PL;

  public:
    LinearProfile(Real c0) : c0_(c0) {};

    Vecd operator()(const Vec2d &position)
    {
        Real x = position[0] / PL_;
        Vecd result = Vec2d::Zero();
        if (x > 0.0)
        {
            result[1] = vf_ * c0_ * (M_ * (cos(kl_ * x) - cosh(kl_ * x)) - N_ * (sin(kl_ * x) - sinh(kl_ * x))) / Q_;
        }
        return result;
    }
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    const paper_bench::BenchmarkDefaults benchmark_defaults{
        PH / 10.0, 1.0, 1.0 / 100.0,
        {{"standard", PH / 10.0},
         {"coarse", PH / 5.0},
         {"fine", PH / 20.0}}};
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av, benchmark_defaults);
    TickCount benchmark_start = TickCount::now();
    global_resolution = Real(benchmark_config.dp);
    BW = global_resolution * 4.0;

    TimeInterval interval_output;
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "oscillating_beam_2d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    interval_output += TickCount::now() - environment_io_start;
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vec2d(-SL - BW, -PL / 2.0),
        Vec2d(PL + 3.0 * BW, PL / 2.0));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox beam_base_shape(
        BoundingBoxd(Vec2d(-SL - BW, -PH / 2 - BW),
                     Vec2d(0.0, PH / 2 + BW)),
        "BeamBase");
    GeometricShapeBox beam_column(
        BoundingBoxd(Vec2d(-SL, -PH / 2),
                     Vec2d(PL, PH / 2)),
        "BeamColumn");
    auto &beam_shape = sph_system.addShape<ComplexShape>("BeamBody");
    beam_shape.add(&beam_base_shape);
    beam_shape.add(&beam_column);
    auto &beam_body = sph_system.addBody<SolidBody>(beam_shape);
    auto &beam_material = beam_body.defineMatterMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam_body.generateParticles<BaseParticles, Lattice>();
    auto &beam_constrain_shape = sph_system.addShape<ComplexShape>("BeamConstrain");
    beam_constrain_shape.add(&beam_base_shape);
    beam_constrain_shape.subtract(&beam_column);
    auto &beam_base = beam_body.addBodyPart<BodyRegionByParticle>(beam_constrain_shape);

    auto &beam_observer = sph_system.addBody<ObserverBody>("BeamObserver");
    beam_observer.defineAdaptationRatios(1.15, 2.0);
    beam_observer.generateParticles<ObserverParticles>(observation_location);
    const size_t initial_particle_count =
        beam_body.getBaseParticles().TotalRealParticles() +
        beam_observer.getBaseParticles().TotalRealParticles();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    auto &beam_body_inner = sph_system.addInnerRelation(beam_body, ConfigType::Lagrangian);
    auto &beam_observer_contact = sph_system.addContactRelation(beam_observer, beam_body, ConfigType::Lagrangian);
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
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    auto &main_methods = sph_solver.getMainMethodContainer();

    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&main_methods.addCellLinkedListDynamics(beam_body));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(beam_body_inner));
    lagrangian_configuration.add(&main_methods.addRelationDynamics(beam_observer_contact));
    lagrangian_configuration.add(&main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(beam_body_inner));

    ParticleDynamicsGroup beam_acoustic_integration;
    beam_acoustic_integration.add(&main_methods.addInteractionDynamicsWithUpdate<
                                   solid_dynamics::StructureNumericalDamping, SaintVenantKirchhoffSolid>(beam_body_inner));
    beam_acoustic_integration.add(&main_methods.addInteractionDynamicsOneLevel<
                                   solid_dynamics::StructureIntegration1stHalfPK2, SaintVenantKirchhoffSolid>(beam_body_inner));
    beam_acoustic_integration.add(&main_methods.addStateDynamics<FixBodyPartConstraintCK>(beam_base));
    beam_acoustic_integration.add(&main_methods.addInteractionDynamicsOneLevel<
                                   solid_dynamics::StructureIntegration2ndHalf>(beam_body_inner));

    auto &acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(beam_body);
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<LinearProfile>>(
                    beam_body, "Velocity", beam_material.ReferenceSoundSpeed())
        .exec();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_displacement = main_methods.addObserveRegression<
        RegressionTestEnsembleAverage, Vecd>(beam_observer_contact, "Position");
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    size_t time_steps = 1;
    int screening_interval = 500;
    int observation_interval = screening_interval;
    const Real end_time = Real(benchmark_config.end_time);
    const Real output_interval = Real(benchmark_config.output_interval);
    const bool verification_mode = !benchmark_config.benchmark_mode;
    const bool write_heavy_output = benchmark_config.output_enabled;
    size_t feature_samples = 0;
    auto &state_recording =
        time_stepper.addTriggerByInterval(output_interval);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    lagrangian_configuration.exec();
    const TimeInterval init_interval =
        TickCount::now() - initialization_start;
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_acoustic_step;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    {
        TickCount initial_io_start = TickCount::now();
        write_displacement.writeToFile(time_steps);
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
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        beam_acoustic_integration.exec(acoustic_dt);
        executed_acoustic_steps++;
        interval_acoustic_step += TickCount::now() - time_instance;

        /** Output body state during the simulation according output_interval. */
        time_instance = TickCount::now();
        /** screen output, write body observables and restart files  */
        if (time_steps == 1 || time_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << time_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (time_steps % observation_interval == 0)
        {
            write_displacement.writeToFile(time_steps);
            feature_samples++;
        }

        if (write_heavy_output && state_recording())
        {
            body_state_recorder.writeToFile();
        }
        interval_output += TickCount::now() - time_instance;

        time_steps++;
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
            write_displacement.generateDataBase(
                Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
            reference_status = "generated";
        }
        else
        {
            write_displacement.testResult();
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
    summary.dp = global_resolution;
    summary.initial_particle_count = initial_particle_count;
    summary.particle_count =
        beam_body.getBaseParticles().TotalRealParticles() +
        beam_observer.getBaseParticles().TotalRealParticles();
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
    summary.extra_fields.push_back({"observer_quantities", "Position"});
    summary.extra_fields.push_back(
        {"observer_sampled", feature_samples > 0 ? "true" : "false"});
    summary.extra_fields.push_back(
        {"feature_samples", std::to_string(feature_samples)});
    summary.extra_fields.push_back(
        {"reference_status", reference_status});
    benchmark_recorder.writeSummary(summary);

    return 0;
}
