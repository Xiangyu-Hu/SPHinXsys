/**
 * @file passive_cantilever.cpp
 * @brief GPU paper B4: 3D passive NeoHookean cantilever (CK / SYCL).
 * @details Port of tests/3d_examples/test_3d_passive_cantilever_neohookean.
 *          Gravity retained via GravityForceCK; classic physical damping has no CK
 *          counterpart yet — omitted and declared in meta/missing.
 *          Duration V2: End_Time default 1.0; PH/dp L1=6 (source), L2=12, ... L5=96.
 */
#include "gpu_paper_recording.h"
#include "sphinxsys.h"
#include "gpu_paper_sycl_group.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace SPH;
using namespace sphinxsys_gpu_paper;

Real total_physical_time = 1.0;
Real PL = 0.1;
Real PH = 0.04;
Real PW = 0.04;
Real SL = 0.02;
Real ph_over_dp = 6.0; // source default PH/dp
Real particle_spacing_ref = PH / ph_over_dp;
Real BW = particle_spacing_ref * 4.0;

Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_cantilever(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
Vec3d domain_lower_bound(-SL - BW, -BW, -BW);
Vec3d domain_upper_bound(PL + BW, PH + BW, PW + BW);
StdVec<Vecd> observation_location = {Vecd(PL, PH, PW)};

Real rho0_s = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5.0e4;
Real physical_viscosity_classic = 50.0; // classic only; not wired in CK
Real gravity_g = 9.8;

void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        ph_over_dp = Real(std::atof(v));
    if (ph_over_dp < Real(2.0))
        ph_over_dp = Real(6.0);

    particle_spacing_ref = PH / ph_over_dp;
    BW = particle_spacing_ref * 4.0;
    // keep geometric SL fixed (source); only spacing/BW follow resolution
    halfsize_cantilever = Vecd(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
    translation_cantilever = Vecd(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
    halfsize_holder = Vecd(0.5 * SL, 0.5 * PH, 0.5 * PW);
    translation_holder = Vecd(-0.5 * SL, 0.5 * PH, 0.5 * PW);
    domain_lower_bound = Vec3d(-SL - BW, -BW, -BW);
    domain_upper_bound = Vec3d(PL + BW, PH + BW, PW + BW);
    observation_location = {Vecd(PL, PH, PW)};
}

int main(int ac, char *av[])
{
    applyResolutionFromEnv();
    applyGpuPaperSyclGroupFromEnv();

    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    auto &cantilever_shape = sph_system.addShape<ComplexShape>("Cantilever");
    cantilever_shape.add<GeometricShapeBox>(Transform(translation_cantilever), halfsize_cantilever);
    auto &holder_shape =
        sph_system.addShape<GeometricShapeBox>(Transform(translation_holder), halfsize_holder, "Holder");
    cantilever_shape.add<GeometricShapeBox>(holder_shape);

    auto &cantilever = sph_system.addBody<SolidBody>(cantilever_shape);
    cantilever.defineMatterMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    cantilever.generateParticles<BaseParticles, Lattice>();
    BodyRegionByParticle holder(cantilever, holder_shape);

    auto &tip_observer = sph_system.addBody<ObserverBody>("TipObserver");
    tip_observer.generateParticles<ObserverParticles>(observation_location);

    auto &cantilever_inner = sph_system.addInnerRelation(cantilever, ConfigType::Lagrangian);
    auto &tip_observer_contact = sph_system.addContactRelation(tip_observer, cantilever, ConfigType::Lagrangian);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();

    auto &cantilever_cell_list = main_methods.addCellLinkedListDynamics(cantilever);
    auto &cantilever_inner_rel = main_methods.addRelationDynamics(cantilever_inner);
    auto &obs_rel = main_methods.addRelationDynamics(tip_observer_contact);
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&cantilever_cell_list);
    lagrangian_configuration.add(&cantilever_inner_rel);
    lagrangian_configuration.add(&obs_rel);

    auto &linear_correction =
        main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(cantilever_inner);
    auto &acoustic_1st = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, NoKernelCorrectionCK>(cantilever_inner);
    auto &acoustic_2nd = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(cantilever_inner);
    auto &acoustic_dt = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(cantilever, 0.5);
    auto &constraint_holder = main_methods.addStateDynamics<FixBodyPartConstraintCK>(holder);

    // Classic uses IncreaseToFullGravity with t_full=0 → constant g from t=0.
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    auto &apply_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(cantilever, gravity);

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &observe_position = main_methods.addObserveRecorder<Vecd>(tip_observer_contact, "Position");
    auto &observe_velocity = main_methods.addObserveRecorder<Vecd>(tip_observer_contact, "Velocity");
    auto &kinetic_energy = main_methods.addReduceDynamics<TotalKineticEnergyCK>(cantilever);

    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    if (const char *end_env = std::getenv("GPU_PAPER_END_TIME"))
        total_physical_time = Real(std::atof(end_env));

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "B4_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "benchmark";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/B4_passive_cantilever/schedules/"
              "B4_benchmark_T1_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    const int n_particles = int(cantilever.getBaseParticles().TotalRealParticles());

    RunRecorder recorder(run_dir, run_id, "B4", mode);
    recorder.setMeta("case_name", "3d_gpu_paper_b4_passive_cantilever");
    recorder.setMeta("source_case", "test_3d_passive_cantilever_neohookean");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("particle_spacing_ref", std::to_string(particle_spacing_ref));
    recorder.setMeta("PH_over_dp", std::to_string(ph_over_dp));
    recorder.setMeta("sycl_work_group_size", gpuPaperGroupMetaString());
    recorder.setMeta("N_particles_lattice", std::to_string(n_particles));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_ck_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("end_time", std::to_string(total_physical_time));
    recorder.setMeta("PL", std::to_string(PL));
    recorder.setMeta("PH", std::to_string(PH));
    recorder.setMeta("PW", std::to_string(PW));
    recorder.setMeta("SL", std::to_string(SL));
    recorder.setMeta("poisson", std::to_string(poisson));
    recorder.setMeta("Youngs_modulus", std::to_string(Youngs_modulus));
    recorder.setMeta("rho0_s", std::to_string(rho0_s));
    recorder.setMeta("gravity_g", std::to_string(gravity_g));
    recorder.setMeta("physical_viscosity_classic", std::to_string(physical_viscosity_classic));
    recorder.setMeta("physical_damping", "omitted_no_CK");
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");
    recorder.setMeta("api", "SPHSolver_NeoHookean_GravityForceCK_NoPhysicalDamping");
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "ck_device_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected", "");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected", "");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "");
    recorder.registerMissing({"physical_damping_CK", "not_collected",
                              "classic DampingProjectionInner not ported to CK; numerical damping only",
                              "case_b4", "optional add CK damping or document difference vs classic"});
    recorder.registerMissing({"strain_energy", "not_collected",
                              "elastic energy reduce not wired; kinetic recorded",
                              "case_b4", "optional"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        for (int j = 0; j <= 50; ++j)
            probe_times.push_back(double(j) * double(total_physical_time) / 50.0);
        recorder.registerMissing({"probe_schedule_json_unreadable", "not_collected",
                                  "failed to parse schedule; fell back to uniform grid",
                                  schedule_path, "fix GPU_PAPER_SCHEDULE_JSON"});
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    auto sample_probes = [&](double t, double requested)
    {
        obs_rel.exec();
        observe_position.writeToFile(0);
        observe_velocity.writeToFile(0);
        auto tip_pos = tip_observer.getBaseParticles().getVariableDataByName<Vecd>("Position");
        auto tip_vel = tip_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
        Real ke = kinetic_energy.exec();

        auto emit = [&](const std::string &qid, Real value, const std::string &unit)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = "tip";
            ps.body = "Cantilever";
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = unit;
            ps.x = double(PL);
            ps.y = double(PH);
            ps.z = double(PW);
            recorder.addProbe(ps);
        };
        emit("pos_x", tip_pos[0][0], "m");
        emit("pos_y", tip_pos[0][1], "m");
        emit("pos_z", tip_pos[0][2], "m");
        emit("vel_x", tip_vel[0][0], "m/s");
        emit("vel_y", tip_vel[0][1], "m/s");
        emit("vel_z", tip_vel[0][2], "m/s");
        emit("deflection_y", tip_pos[0][1] - PH, "m");

        ProbeSample pke;
        pke.physical_time_s = t;
        pke.requested_time_s = requested;
        pke.probe_id = "body";
        pke.body = "Cantilever";
        pke.quantity = "kinetic_energy";
        pke.value = double(ke);
        pke.unit = "J";
        recorder.addProbe(pke);

        recorder.addDiagnostic(t, "Cantilever", "N_active",
                               double(cantilever.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "Cantilever", "tip_y", double(tip_pos[0][1]), "m", "measured");
        recorder.addDiagnostic(t, "Cantilever", "deflection_y", double(tip_pos[0][1] - PH), "m", "measured");
        recorder.addDiagnostic(t, "Cantilever", "kinetic_energy", double(ke), "J", "measured");
    };

    size_t acoustic_steps = 1;
    int screening_interval = 100;
    auto wall0 = std::chrono::steady_clock::now();
    TimeInterval interval_output;
    TimeInterval interval_acoustic_step;

    std::cout << std::fixed << std::setprecision(6)
              << "[B4] PH/dp=" << ph_over_dp << " dp=" << particle_spacing_ref
              << " N=" << n_particles << " End_Time=" << total_physical_time
              << " damping=omitted_no_CK\n";

    lagrangian_configuration.exec();
    linear_correction.exec();

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
    const Real dt_min = Real(1.0e-12);

    while (!time_stepper.isEndTime(total_physical_time))
    {
        TickCount time_instance = TickCount::now();
        Real step_dt = time_stepper.incrementPhysicalTime(acoustic_dt);
        if (!(step_dt > dt_min) || !std::isfinite(double(step_dt)))
        {
            std::cout << "ERROR: bad acoustic_dt=" << step_dt
                      << " at N=" << acoustic_steps
                      << " Time=" << time_stepper.getPhysicalTime() << std::endl;
            dt_fail = 1;
            break;
        }

        apply_gravity.exec();
        acoustic_1st.exec(step_dt);
        constraint_holder.exec();
        acoustic_2nd.exec(step_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        if (acoustic_steps == 1 || acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
                      << "	Time = " << time_stepper.getPhysicalTime()
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        double t_now = time_stepper.getPhysicalTime();
        double requested = 0.0;
        if (dueForSample(t_now, probe_times, probe_cursor, &requested))
            sample_probes(t_now, requested);
        if (enable_vtp && dueForSample(t_now, field_times, field_cursor, &requested))
            body_state_recorder.writeToFile();
        interval_output += TickCount::now() - time_instance;

        acoustic_steps++;
    }

    if (!dt_fail)
        sample_probes(time_stepper.getPhysicalTime(), double(total_physical_time));

    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9)
              << "interval_acoustic_step = " << interval_acoustic_step.seconds() << "\n";

    recorder.stages().add("acoustic_step", interval_acoustic_step.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("acoustic", int(acoustic_steps), 0, 0, 0, 0, 0, 0, "solid_ck", "measured");
    recorder.addBodyStats(time_stepper.getPhysicalTime(), "Cantilever",
                          int(cantilever.getBaseParticles().TotalRealParticles()),
                          int(cantilever.getBaseParticles().ParticlesBound()),
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
