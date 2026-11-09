/**
 * @file  oscillating_beam.cpp
 * @brief GPU paper A2: 2D oscillating beam (extra_source_and_tests copy).
 * @details Duration V2 T=2 s; tip probes + kinetic energy; VTP on bound field schedule.
 * @note    GPU paper experiment (docs/gpu_paper_discussion).
 */
#include "gpu_paper_recording.h"
#include "sphinxsys.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace SPH;
using namespace sphinxsys_gpu_paper;

//----------------------------------------------------------------------
// Geometry / material (ENOG-matched source values; algorithm may differ).
//----------------------------------------------------------------------
Real PL = 0.2;
Real PH = 0.02;
Real SL = 0.06;
Real ph_over_dp = 10.0; // L1=10, L2=20, L3=40
Real global_resolution = PH / ph_over_dp;
Real BW = global_resolution * 4;
BoundingBoxd system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0));

Real rho0_s = 1.0e3;
Real Youngs_modulus = 2.0e6;
Real poisson = 0.3975;

Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;

StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};

void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        ph_over_dp = Real(std::atof(v));
    if (ph_over_dp < Real(5.0))
        ph_over_dp = Real(10.0);
    global_resolution = PH / ph_over_dp;
    BW = global_resolution * 4;
    system_domain_bounds = BoundingBoxd(Vec2d(-SL - BW, -PL / 2.0), Vec2d(PL + 3.0 * BW, PL / 2.0));
}

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

int main(int ac, char *av[])
{
    applyResolutionFromEnv();

    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    GeometricShapeBox beam_base_shape(
        BoundingBoxd(Vec2d(-SL - BW, -PH / 2 - BW), Vec2d(0.0, PH / 2 + BW)), "BeamBase");
    GeometricShapeBox beam_column(
        BoundingBoxd(Vec2d(-SL, -PH / 2), Vec2d(PL, PH / 2)), "BeamColumn");

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

    auto &beam_body_inner = sph_system.addInnerRelation(beam_body, ConfigType::Lagrangian);
    auto &beam_observer_contact = sph_system.addContactRelation(beam_observer, beam_body, ConfigType::Lagrangian);

    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.getHostMethodContainer();
    auto &main_methods = sph_solver.getMainMethodContainer();

    auto &beam_cell_list = main_methods.addCellLinkedListDynamics(beam_body);
    auto &beam_inner_rel = main_methods.addRelationDynamics(beam_body_inner);
    auto &beam_obs_rel = main_methods.addRelationDynamics(beam_observer_contact);
    auto &beam_correction = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(beam_body_inner);

    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&beam_cell_list);
    lagrangian_configuration.add(&beam_inner_rel);
    lagrangian_configuration.add(&beam_obs_rel);
    lagrangian_configuration.add(&beam_correction);

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

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &observe_tip_position = main_methods.addObserveRecorder<Vecd>(beam_observer_contact, "Position");
    auto &observe_tip_velocity = main_methods.addObserveRecorder<Vecd>(beam_observer_contact, "Velocity");
    auto &kinetic_energy = main_methods.addReduceDynamics<TotalKineticEnergyCK>(beam_body);

    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    const char *end_env = std::getenv("GPU_PAPER_END_TIME");
    Real end_time = end_env ? Real(std::atof(end_env)) : Real(2.0);

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "A2_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "validate";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/A2_oscillating_beam/schedules/"
              "A2_validate_T2_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "A2", mode);
    recorder.setMeta("case_name", "2d_gpu_paper_a2_oscillating_beam");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("global_resolution", std::to_string(global_resolution));
    recorder.setMeta("PH_over_dp", std::to_string(ph_over_dp));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("end_time", std::to_string(end_time));
    recorder.setMeta("PL", std::to_string(PL));
    recorder.setMeta("PH", std::to_string(PH));
    recorder.setMeta("Youngs_modulus", std::to_string(Youngs_modulus));
    recorder.setMeta("poisson", std::to_string(poisson));
    recorder.setMeta("rho0_s", std::to_string(rho0_s));
    recorder.setMeta("vf", std::to_string(vf));
    recorder.setMeta("schedule_json", schedule_path);
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "no_explicit_transfer_API_in_this_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected",
                     "detail deferred to profile short window policy");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected",
                                         "only_on_B2_B3_profile_window");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "needs_valid_telemetry_coverage");
    recorder.registerMissing({"strain_energy", "not_collected",
                              "no dedicated CK elastic-energy reduce wired yet; kinetic recorded",
                              "case_a2", "add elastic energy reduce if required for paper table"});
    recorder.registerMissing({"ENOG_reference_error", "not_collected",
                              "ENOG period/amplitude digitization not bound into this run",
                              "case_a2", "post-process tip probes vs ENOG table after formal runs"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        for (int j = 0; j < 100; ++j)
            probe_times.push_back(double(j) * double(end_time) / 99.0);
        recorder.registerMissing({"probe_schedule_json_unreadable", "not_collected",
                                  "failed to parse schedule; fell back to uniform100 of end_time",
                                  schedule_path, "fix GPU_PAPER_SCHEDULE_JSON path"});
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    auto sample_probes = [&](double t, double requested)
    {
        // Total-Lagrangian: do NOT rebuild LinearCorrectionMatrix here — only refresh observer contact.
        beam_obs_rel.exec();
        observe_tip_position.writeToFile(0);
        observe_tip_velocity.writeToFile(0);
        auto tip_pos = beam_observer.getBaseParticles().getVariableDataByName<Vecd>("Position");
        auto tip_vel = beam_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
        Real x = tip_pos[0][0];
        Real y = tip_pos[0][1];
        Real vx = tip_vel[0][0];
        Real vy = tip_vel[0][1];
        Real ke = kinetic_energy.exec();

        auto emit = [&](const std::string &qid, Real value, const std::string &unit)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = "tip";
            ps.body = "BeamBody";
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = unit;
            ps.x = double(PL);
            ps.y = 0.0;
            recorder.addProbe(ps);
        };
        emit("tip_x", x, "m");
        emit("tip_y", y, "m");
        emit("tip_vx", vx, "m/s");
        emit("tip_vy", vy, "m/s");
        emit("kinetic_energy", ke, "J");

        recorder.addDiagnostic(t, "BeamBody", "N_active",
                               double(beam_body.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "BeamBody", "tip_y", double(y), "m", "measured");
        recorder.addDiagnostic(t, "BeamBody", "kinetic_energy", double(ke), "J", "measured");
    };

    size_t time_steps = 1;
    int screening_interval = 500;
    auto wall0 = std::chrono::steady_clock::now();

    // Reference configuration once (matches upstream oscillating_beam_sycl).
    lagrangian_configuration.exec();

    TimeInterval interval_output;
    TimeInterval interval_acoustic_step;

    if (enable_vtp)
        body_state_recorder.writeToFile();
    {
        double req0 = probe_times.empty() ? 0.0 : probe_times[0];
        sample_probes(time_stepper.getPhysicalTime(), req0);
        if (!probe_times.empty())
            probe_cursor = 1;
        if (enable_vtp && !field_times.empty() && field_times[0] <= 0.0)
            field_cursor = 1;
    }

    TickCount t0 = TickCount::now();
    int dt_fail = 0;
    while (!time_stepper.isEndTime(end_time))
    {
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        if (!(acoustic_dt > Real(0.0)) || !std::isfinite(double(acoustic_dt)))
        {
            std::cout << "ERROR: non-positive/non-finite acoustic_dt=" << acoustic_dt
                      << " at N=" << time_steps
                      << " Time=" << time_stepper.getPhysicalTime() << std::endl;
            dt_fail = 1;
            break;
        }
        beam_acoustic_integration.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        if (time_steps == 1 || time_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << time_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        double t_now = time_stepper.getPhysicalTime();
        double requested = 0.0;
        if (dueForSample(t_now, probe_times, probe_cursor, &requested))
            sample_probes(t_now, requested);
        if (enable_vtp && dueForSample(t_now, field_times, field_cursor, &requested))
            body_state_recorder.writeToFile();

        interval_output += TickCount::now() - time_instance;
        time_steps++;
    }

    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";

    recorder.stages().add("acoustic_step", interval_acoustic_step.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("acoustic", int(time_steps), 0, 0, 0, 0, 0, 0, "solid_ck", "measured");
    recorder.addBodyStats(time_stepper.getPhysicalTime(), "BeamBody",
                          int(beam_body.getBaseParticles().TotalRealParticles()),
                          int(beam_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    const char *status = dt_fail ? "timestep_collapse" : "ok";
    recorder.finalize(status, 0.0, time_stepper.getPhysicalTime(), wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval_output.seconds(),
                      0.0, interval_output.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return dt_fail;
}
