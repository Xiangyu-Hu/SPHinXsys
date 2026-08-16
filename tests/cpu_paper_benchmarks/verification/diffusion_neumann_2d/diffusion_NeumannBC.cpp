/**
 * @file 	diffusion_NeumannBC.cpp
 * @brief 	2D test of diffusion problem with Neumann boundary condition.
 * @details This is the first case to validate multiple boundary conditions.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include <algorithm>
#include <cmath>
#include <memory>

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real global_resolution = H / 100.0;
Real BW = global_resolution * 2.0;
//----------------------------------------------------------------------
//	Basic parameters for diffusion properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
std::string species_name = "Phi";
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 100.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real heat_flux = 900.0; // from the Neumann boundary
//----------------------------------------------------------------------
//	Generate 2D geometrics used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> thermal_domain_edge_points{
    Vecd(0.0, 0.0), Vecd(0.0, H),
    Vecd(L, H), Vecd(L, 0.0), Vecd(0.0, 0.0)};

std::vector<Vecd> left_region_edge_points{
    Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW), Vecd(0.4 * L, H + BW),
    Vecd(0.4 * L, H), Vecd(0.3 * L, H)};

std::vector<Vecd> right_region_edge_points{
    Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW), Vecd(0.7 * L, H + BW),
    Vecd(0.7 * L, H), Vecd(0.6 * L, H)};

std::vector<Vecd> heat_flux_region_edge_points{
    Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0), Vecd(0.55 * L, 0),
    Vecd(0.55 * L, -BW), Vecd(0.45 * L, -BW)};
//----------------------------------------------------------------------
//	Define Shapes used for SPH bodies and body-parts.
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(thermal_domain_edge_points, GeometricOps::add);
    }
};

class DirichletWallBoundary : public MultiPolygonShape
{
  public:
    explicit DirichletWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(left_region_edge_points, GeometricOps::add);
        multi_polygon_.addPolygon(right_region_edge_points, GeometricOps::add);
    }
};

class NeumannWallBoundary : public MultiPolygonShape
{
  public:
    explicit NeumannWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(heat_flux_region_edge_points, GeometricOps::add);
    }
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the middle line. */
    size_t number_of_observation_points = 5;
    Real range_of_measure = L;
    Real start_of_measure = 0;

    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(0.5 * L, range_of_measure * Real(i) /
                                                Real(number_of_observation_points - 1) +
                                            start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    const paper_bench::BenchmarkDefaults benchmark_defaults{
        H / 100.0, 1.0, 0.1,
        {{"standard", H / 100.0},
         {"coarse", H / 50.0},
         {"medium", H / 100.0},
         {"fine", H / 200.0}}};
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av, benchmark_defaults);
    TickCount benchmark_start = TickCount::now();
    global_resolution = benchmark_config.dp;
    BW = global_resolution * 2.0;
    left_region_edge_points = {
        Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW),
        Vecd(0.4 * L, H + BW), Vecd(0.4 * L, H), Vecd(0.3 * L, H)};
    right_region_edge_points = {
        Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW),
        Vecd(0.7 * L, H + BW), Vecd(0.7 * L, H), Vecd(0.6 * L, H)};
    heat_flux_region_edge_points = {
        Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0),
        Vecd(0.55 * L, 0), Vecd(0.55 * L, -BW),
        Vecd(0.45 * L, -BW)};
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "diffusion_neumann_2d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    TimeInterval io_interval = TickCount::now() - environment_io_start;
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    diffusion_body.defineMatterMaterial<Solid>();
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_Dirichlet(sph_system, makeShared<DirichletWallBoundary>("DirichletWallBoundary"));
    wall_Dirichlet.defineMatterMaterial<Solid>();
    wall_Dirichlet.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_Neumann(sph_system, makeShared<NeumannWallBoundary>("NeumannWallBoundary"));
    wall_Neumann.defineMatterMaterial<Solid>();
    wall_Neumann.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    const size_t initial_particle_count =
        diffusion_body.getBaseParticles().TotalRealParticles() +
        wall_Dirichlet.getBaseParticles().TotalRealParticles() +
        wall_Neumann.getBaseParticles().TotalRealParticles() +
        temperature_observer.getBaseParticles().TotalRealParticles();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    MultiPolygonShape left_region_shape(MultiPolygon(left_region_edge_points), "LeftRegion");
    BodyRegionByParticle wall_Dirichlet_left_region(wall_Dirichlet, left_region_shape);

    MultiPolygonShape right_region_shape(MultiPolygon(right_region_edge_points), "RightRegion");
    BodyRegionByParticle wall_Dirichlet_right_region(wall_Dirichlet, right_region_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> diffusion_body_contact(diffusion_body, {&wall_Dirichlet, &wall_Neumann});
    Contact<> temperature_observer_contact(temperature_observer, {&diffusion_body});
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
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_Neumann);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> diffusion_body_cell_linked_list(diffusion_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_Dirichlet_cell_linked_list(wall_Dirichlet);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_Neumann_cell_linked_list(wall_Neumann);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>>
        diffusion_body_update_complex_relation(diffusion_body_inner, diffusion_body_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> observer_contact_relation(temperature_observer_contact);

    IsotropicDiffusion isotropic_diffusion(species_name, diffusion_coeff);
    auto contact_dirichlet_view = makeRelationView(diffusion_body_contact, wall_Dirichlet);
    auto contact_neumann_view = makeRelationView(diffusion_body_contact, wall_Neumann);

    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, &isotropic_diffusion);
    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Neumann<IsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Neumann<IsotropicDiffusion>, NoKernelCorrectionCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(diffusion_body_inner, &isotropic_diffusion),
                                 DynamicsArgs(contact_dirichlet_view, &isotropic_diffusion),
                                 DynamicsArgs(contact_neumann_view, &isotropic_diffusion));
    //----------------------------------------------------------------------
    //	Specify initial condition if necessary.
    //----------------------------------------------------------------------
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, ConstantValue<Real>>>
        diffusion_initial_condition(diffusion_body, species_name, initial_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<BodyRegionByParticle, ConstantValue<Real>>>
        left_initial_condition(wall_Dirichlet_left_region, species_name, left_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<BodyRegionByParticle, ConstantValue<Real>>>
        right_initial_condition(wall_Dirichlet_right_region, species_name, right_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<SPHBody, ConstantValue<Real>>>
        wall_Neumann_initial_condition(wall_Neumann, species_name + "Flux", heat_flux);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_states(sph_system);
    using TemperatureObservation =
        ObservedQuantityRecording<MainExecutionPolicy, Real,
                                  RestoringCorrection>;
    using TemperatureRegression =
        RegressionTestEnsembleAverage<TemperatureObservation>;
    std::unique_ptr<TemperatureObservation> benchmark_temperature_observer;
    std::unique_ptr<TemperatureRegression> verification_temperature_observer;
    if (benchmark_config.benchmark_mode)
    {
        benchmark_temperature_observer =
            std::make_unique<TemperatureObservation>(
                temperature_observer_contact, species_name);
    }
    else
    {
        verification_temperature_observer =
            std::make_unique<TemperatureRegression>(
                temperature_observer_contact, species_name);
    }
    auto write_temperature = [&](size_t iteration)
    {
        if (benchmark_temperature_observer)
        {
            benchmark_temperature_observer->writeToFile(iteration);
        }
        else
        {
            verification_temperature_observer->writeToFile(iteration);
        }
    };
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();

    diffusion_body_cell_linked_list.exec();
    wall_Dirichlet_cell_linked_list.exec();
    wall_Neumann_cell_linked_list.exec();

    diffusion_body_update_complex_relation.exec();
    observer_contact_relation.exec();

    diffusion_initial_condition.exec();
    left_initial_condition.exec();
    right_initial_condition.exec();
    wall_Neumann_initial_condition.exec();
    const TimeInterval init_interval =
        TickCount::now() - initialization_start;
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = Real(benchmark_config.end_time);
    Real Observe_time = 0.01 * end_time;
    Real Output_Time = Real(benchmark_config.output_interval);
    Real dt = 0.0;
    size_t number_of_outer_steps = 0;
    const bool verification_mode = !benchmark_config.benchmark_mode;
    const bool write_heavy_output = benchmark_config.output_enabled;
    size_t feature_samples = 0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    {
        TickCount initial_io_start = TickCount::now();
        write_temperature(ite);
        feature_samples++;
        if (write_heavy_output)
        {
            write_states.writeToFile();
        }
        io_interval += TickCount::now() - initial_io_start;
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time &&
               sv_physical_time->getValue() < end_time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time &&
                   integration_time < Output_Time &&
                   sv_physical_time->getValue() < end_time)
            {
                dt = get_time_step_size.exec();
                if (!std::isfinite(dt) || dt <= 0.0)
                {
                    throw std::runtime_error(
                        "Diffusion time-step estimator returned a non-positive or non-finite value");
                }
                dt = std::min(
                    {dt,
                     Observe_time - relaxation_time,
                     Output_Time - integration_time,
                     end_time - sv_physical_time->getValue()});
                if (ite % 500 == 0)
                {
                    TickCount console_io_start = TickCount::now();
                    std::cout << "N=" << ite << " Time: "
                              << sv_physical_time->getValue() << "	dt: "
                              << dt << "\n";
                    io_interval += TickCount::now() - console_io_start;
                }

                diffusion_relaxation_rk2.exec(dt);

                ite++;
                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
            }
        }

        {
            TickCount periodic_io_start = TickCount::now();
            write_temperature(ite);
            feature_samples++;
            if (write_heavy_output)
            {
                write_states.writeToFile();
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
            verification_temperature_observer->generateDataBase(1.0e-3, 1.0e-3);
            reference_status = "generated";
        }
        else if (sph_system.RestartStep() == 0)
        {
            verification_temperature_observer->testResult();
            reference_status = "tested";
        }
        else
        {
            reference_status = "not_checked_restart";
        }
        io_interval += TickCount::now() - regression_io_start;
    }
    else if (verification_mode)
    {
        reference_status = "not_checked_output_disabled";
    }

    TickCount final_console_io_start = TickCount::now();
    std::cout << "Total physical time for computation: "
              << sv_physical_time->getValue() << " seconds." << std::endl;
    io_interval += TickCount::now() - final_console_io_start;

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
        diffusion_body.getBaseParticles().TotalRealParticles() +
        wall_Dirichlet.getBaseParticles().TotalRealParticles() +
        wall_Neumann.getBaseParticles().TotalRealParticles() +
        temperature_observer.getBaseParticles().TotalRealParticles();
    summary.physical_time = sv_physical_time->getValue();
    summary.outer_steps = number_of_outer_steps;
    summary.wall_seconds = total_wall.seconds();
    summary.compute_seconds = compute_seconds;
    summary.io_seconds = io_interval.seconds();
    summary.init_seconds = init_interval.seconds();
    summary.time_per_outer_step =
        number_of_outer_steps == 0
            ? 0.0
            : compute_seconds / static_cast<double>(number_of_outer_steps);
    summary.setParticleUpdates(
        summary.particle_count * static_cast<std::size_t>(ite),
        "particle_count*diffusion_steps");
    summary.diffusion_seconds = paper_bench::formatMetric(compute_seconds);
    summary.extra_fields.push_back({"diffusion_steps", std::to_string(ite)});
    summary.extra_fields.push_back(
        {"diffusion_time_per_step",
         ite == 0 ? "0"
                  : std::to_string(compute_seconds /
                                   static_cast<double>(ite))});
    summary.extra_fields.push_back({"observer_available", "true"});
    summary.extra_fields.push_back(
        {"observer_sampled", feature_samples > 0 ? "true" : "false"});
    summary.extra_fields.push_back(
        {"feature_samples", std::to_string(feature_samples)});
    summary.extra_fields.push_back({"reference_status", reference_status});
    benchmark_recorder.writeSummary(summary);
    return 0;
}