/**
 * @file  twisting_column.cpp
 * @brief GPU paper B3: 3D twisting column (extra_source copy of SYCL CK NeoHookean case).
 * @details Duration V2 End_Time=0.5; PH/dp L1..L5; tip + off-axis probes; gated VTP.
 * @note    Source nu=0.45 / angular=-300; ENOG paper nu differs — declare when comparing.
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
Real total_physical_time = 0.5;
Real PL = 6.0;
Real PH = 1.0;
Real PW = 1.0;
Real ph_over_dp = 10.0; // L2 source default: PH/dp=10 → dp=0.1
Real particle_spacing_ref = PH / ph_over_dp;
Real BW = particle_spacing_ref * 0.0;
Real SL = particle_spacing_ref * 1.0;
Vecd halfsize_column(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_column(0.5 * (PL - SL), 0.0, 0.0);
Vecd halfsize_holder(0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
Vecd translation_holder(-0.5 * (SL + BW), 0.0, 0.0);
Vec3d domain_lower_bound(-SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
Vec3d domain_upper_bound(PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
// tip center + off-axis material-space probe (Duration V2)
StdVec<Vecd> observation_location = {Vecd(PL, 0.0, 0.0), Vecd(PL, 0.35, 0.0)};
//----------------------------------------------------------------------
Real rho0_s = 1100.0;
Real poisson = 0.45;
Real Youngs_modulus = 1.7e7;
Real angular_0 = -300.0;

void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        ph_over_dp = Real(std::atof(v));
    if (ph_over_dp < Real(3.0))
        ph_over_dp = Real(10.0);
    particle_spacing_ref = PH / ph_over_dp;
    BW = particle_spacing_ref * 0.0;
    SL = particle_spacing_ref * 1.0;
    halfsize_column = Vecd(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
    translation_column = Vecd(0.5 * (PL - SL), 0.0, 0.0);
    halfsize_holder = Vecd(0.5 * (SL + BW), 0.5 * (PH + BW), 0.5 * (PW + BW));
    translation_holder = Vecd(-0.5 * (SL + BW), 0.0, 0.0);
    domain_lower_bound = Vec3d(-SL - BW, -0.5 * (PH + BW), -0.5 * (PW + BW));
    domain_upper_bound = Vec3d(PL, 0.5 * (PH + BW), 0.5 * (PW + BW));
    observation_location = {Vecd(PL, 0.0, 0.0), Vecd(PL, 0.35, 0.0)};
}

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
        return Vecd::Zero();
    };
};

int main(int ac, char *av[])
{
    applyResolutionFromEnv();

    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    auto &column_shape = sph_system.addShape<ComplexShape>("Column");
    column_shape.add<GeometricShapeBox>(Transform(translation_column), halfsize_column);
    auto &holder_shape = sph_system.addShape<GeometricShapeBox>(Transform(translation_holder), halfsize_holder, "Holder");
    column_shape.add<GeometricShapeBox>(holder_shape);

    auto &column = sph_system.addBody<SolidBody>(column_shape);
    column.defineMatterMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    column.generateParticles<BaseParticles, Lattice>();
    BodyRegionByParticle holder(column, holder_shape);

    auto &my_observer = sph_system.addBody<ObserverBody>("MyObserver");
    my_observer.generateParticles<ObserverParticles>(observation_location);

    auto &column_inner = sph_system.addInnerRelation(column, ConfigType::Lagrangian);
    auto &my_observer_contact = sph_system.addContactRelation(my_observer, column, ConfigType::Lagrangian);

    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.getHostMethodContainer();
    auto &main_methods = sph_solver.getMainMethodContainer();

    auto &column_cell_list = main_methods.addCellLinkedListDynamics(column);
    auto &column_inner_rel = main_methods.addRelationDynamics(column_inner);
    auto &obs_rel = main_methods.addRelationDynamics(my_observer_contact);
    ParticleDynamicsGroup lagrangian_configuration;
    lagrangian_configuration.add(&column_cell_list);
    lagrangian_configuration.add(&column_inner_rel);
    lagrangian_configuration.add(&obs_rel);

    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);
    auto &column_acoustic_step_1st_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, NeoHookeanSolid, NoKernelCorrectionCK>(column_inner);
    auto &column_acoustic_step_2nd_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(column_inner);

    auto &column_acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(column, 0.5);
    auto &constraint_holder = main_methods.addStateDynamics<FixBodyPartConstraintCK>(holder);
    host_methods.addStateDynamics<VariableAssignment, SpatialDistribution<VelocityProfile>>(column, "Velocity").exec();

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &observe_position = main_methods.addObserveRecorder<Vecd>(my_observer_contact, "Position");
    auto &observe_velocity = main_methods.addObserveRecorder<Vecd>(my_observer_contact, "Velocity");
    auto &kinetic_energy = main_methods.addReduceDynamics<TotalKineticEnergyCK>(column);

    TimeStepper &time_stepper = sph_solver.getTimeStepper();

    if (const char *end_env = std::getenv("GPU_PAPER_END_TIME"))
        total_physical_time = Real(std::atof(end_env));

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "B3_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "benchmark";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/B3_twisting_column/schedules/"
              "B3_benchmark_T0p5_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "B3", mode);
    recorder.setMeta("case_name", "3d_gpu_paper_b3_twisting_column");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("particle_spacing_ref", std::to_string(particle_spacing_ref));
    recorder.setMeta("PH_over_dp", std::to_string(ph_over_dp));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("end_time", std::to_string(total_physical_time));
    recorder.setMeta("PL", std::to_string(PL));
    recorder.setMeta("PH", std::to_string(PH));
    recorder.setMeta("PW", std::to_string(PW));
    recorder.setMeta("poisson", std::to_string(poisson));
    recorder.setMeta("Youngs_modulus", std::to_string(Youngs_modulus));
    recorder.setMeta("rho0_s", std::to_string(rho0_s));
    recorder.setMeta("angular_0", std::to_string(angular_0));
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");
    recorder.setMeta("material_note", "NeoHookean_nu0.45_not_ENOG_nu0.49");
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "no_explicit_transfer_API_in_this_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected",
                     "detail deferred to profile short window policy");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected",
                                         "only_on_B2_B3_profile_window");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "needs_valid_telemetry_coverage");
    recorder.registerMissing({"strain_energy", "not_collected",
                              "no dedicated NeoHookean elastic-energy reduce wired yet; kinetic recorded",
                              "case_b3", "add elastic energy reduce if required"});
    recorder.registerMissing({"section_twist_angle_fit", "not_collected",
                              "end-face particle ID fit not wired; off-axis observer recorded instead",
                              "case_b3", "optional post-process persistent-ID section angle"});
    recorder.registerMissing({"ENOG_reference_error", "not_collected",
                              "source nu=0.45 vs paper 0.49/0.499; not same-config quantitative",
                              "case_b3", "declare mismatch or freeze matched config separately"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        if (mode == "benchmark")
        {
            for (int j = 0; j <= 50; ++j)
                probe_times.push_back(0.01 * double(j));
        }
        else
        {
            for (int j = 0; j < 100; ++j)
                probe_times.push_back(double(j) * double(total_physical_time) / 99.0);
        }
        recorder.registerMissing({"probe_schedule_json_unreadable", "not_collected",
                                  "failed to parse schedule; fell back to mode default grid",
                                  schedule_path, "fix GPU_PAPER_SCHEDULE_JSON path"});
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    auto sample_probes = [&](double t, double requested)
    {
        obs_rel.exec();
        observe_position.writeToFile(0);
        observe_velocity.writeToFile(0);
        auto tip_pos = my_observer.getBaseParticles().getVariableDataByName<Vecd>("Position");
        auto tip_vel = my_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
        Real ke = kinetic_energy.exec();

        auto emit = [&](const std::string &pid, const std::string &qid, Real value,
                        const std::string &unit, Real x, Real y, Real z)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = pid;
            ps.body = "Column";
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = unit;
            ps.x = double(x);
            ps.y = double(y);
            ps.z = double(z);
            recorder.addProbe(ps);
        };
        // tip center
        emit("tip_center", "pos_x", tip_pos[0][0], "m", PL, 0, 0);
        emit("tip_center", "pos_y", tip_pos[0][1], "m", PL, 0, 0);
        emit("tip_center", "pos_z", tip_pos[0][2], "m", PL, 0, 0);
        emit("tip_center", "vel_x", tip_vel[0][0], "m/s", PL, 0, 0);
        emit("tip_center", "vel_y", tip_vel[0][1], "m/s", PL, 0, 0);
        emit("tip_center", "vel_z", tip_vel[0][2], "m/s", PL, 0, 0);
        // off-axis
        if (my_observer.getBaseParticles().TotalRealParticles() > 1)
        {
            emit("tip_offaxis", "pos_x", tip_pos[1][0], "m", PL, 0.35, 0);
            emit("tip_offaxis", "pos_y", tip_pos[1][1], "m", PL, 0.35, 0);
            emit("tip_offaxis", "pos_z", tip_pos[1][2], "m", PL, 0.35, 0);
            emit("tip_offaxis", "vel_y", tip_vel[1][1], "m/s", PL, 0.35, 0);
            emit("tip_offaxis", "vel_z", tip_vel[1][2], "m/s", PL, 0.35, 0);
            Real r = std::sqrt(double(tip_pos[1][1] * tip_pos[1][1] + tip_pos[1][2] * tip_pos[1][2]));
            Real ang = std::atan2(double(tip_pos[1][1]), double(tip_pos[1][2]));
            emit("tip_offaxis", "radius", Real(r), "m", PL, 0.35, 0);
            emit("tip_offaxis", "atan2_yz", Real(ang), "rad", PL, 0.35, 0);
        }
        emit("column", "kinetic_energy", ke, "J", 0, 0, 0);

        recorder.addDiagnostic(t, "Column", "N_active",
                               double(column.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "Column", "kinetic_energy", double(ke), "J", "measured");
        recorder.addDiagnostic(t, "Column", "tip_offaxis_y", double(tip_pos[1][1]), "m", "measured");
    };

    size_t acoustic_steps = 1;
    int screening_interval = 100;
    auto wall0 = std::chrono::steady_clock::now();
    TimeInterval interval_output;
    TimeInterval interval_acoustic_step;

    lagrangian_configuration.exec();
    column_linear_correction_matrix.exec();

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
        Real acoustic_dt = time_stepper.incrementPhysicalTime(column_acoustic_time_step);
        if (!(acoustic_dt > dt_min) || !std::isfinite(double(acoustic_dt)))
        {
            std::cout << "ERROR: bad acoustic_dt=" << acoustic_dt
                      << " at N=" << acoustic_steps
                      << " Time=" << time_stepper.getPhysicalTime() << std::endl;
            dt_fail = 1;
            break;
        }
        column_acoustic_step_1st_half.exec(acoustic_dt);
        constraint_holder.exec();
        column_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        if (acoustic_steps == 1 || acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
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

        acoustic_steps++;
    }

    if (!dt_fail)
        sample_probes(time_stepper.getPhysicalTime(), double(total_physical_time));

    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";

    recorder.stages().add("acoustic_step", interval_acoustic_step.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("acoustic", int(acoustic_steps), 0, 0, 0, 0, 0, 0, "solid_ck", "measured");
    recorder.addBodyStats(time_stepper.getPhysicalTime(), "Column",
                          int(column.getBaseParticles().TotalRealParticles()),
                          int(column.getBaseParticles().ParticlesBound()),
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
