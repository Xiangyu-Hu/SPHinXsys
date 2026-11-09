/**
 * @file  dambreak.cpp
 * @brief GPU paper B2: 3D dambreak (extra_source copy of tests_sycl CK path).
 * @details Duration V2 End_Time=10; LH/dp L1..L5; benchmark window51 / validate denser probes.
 * @note    Geometry matches CPC appendix (DL=5.366, DH=2, DW=0.5, probe y=0.2 mid-width).
 */
#include "gpu_paper_recording.h"
#include "sphinxsys.h"
#include "gpu_paper_sycl_group.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace SPH;
using namespace sphinxsys_gpu_paper;
//----------------------------------------------------------------------
Real global_resolution = 0.05;
Real BW = global_resolution * 4;
Real DL = 5.366;
Real DH = 2.0;
Real DW = 0.5;
Real LL = 2.0;
Real LH = 1.0;
Real LW = 0.5;
Real lh_over_dp = 20.0; // L2 default = source 0.05
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f; // =20 with H=g=1 → CPC sound speed
//----------------------------------------------------------------------
void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        lh_over_dp = Real(std::atof(v));
    if (lh_over_dp < Real(5.0))
        lh_over_dp = Real(20.0);
    global_resolution = LH / lh_over_dp;
    BW = global_resolution * 4;
}

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
    observation_points.push_back(Vecd(DL, 0.2, 0.5 * DW)); // CPC primary
    observation_points.push_back(Vecd(DL, 0.24, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.252, 0.5 * DW));
    observation_points.push_back(Vecd(DL, 0.266, 0.5 * DW));
    return observation_points;
}

int main(int ac, char *av[])
{
    applyResolutionFromEnv();
    applyGpuPaperSyclGroupFromEnv();

    BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    WaterBlock initial_water_block("WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(createObservationPoints());

    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, wall_boundary);
    Contact<> fluid_observer_contact(fluid_observer, water_block);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_block_update_complex_relation(water_block_inner, water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(fluid_observer_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);

    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    StateDynamics<MainExecutionPolicy, GravityForceCK<Gravity>> constant_gravity(water_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
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

    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "PositionDivergence");

    ReducedQuantityRecording<MainExecutionPolicy, TotalMechanicalEnergyCK>
        record_water_mechanical_energy(water_block, gravity);
    ObservedQuantityRecording<MainExecutionPolicy, Real>
        fluid_observer_pressure(fluid_observer_contact, "Pressure");

    const char *end_env = std::getenv("GPU_PAPER_END_TIME");
    Real end_time = end_env ? Real(std::atof(end_env)) : Real(10.0);

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "B2_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "benchmark";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/B2_dambreak/schedules/"
              "B2_benchmark_T10_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "B2", mode);
    recorder.setMeta("case_name", "3d_gpu_paper_b2_dambreak");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("global_resolution", std::to_string(global_resolution));
    recorder.setMeta("LH_over_dp", std::to_string(lh_over_dp));
    recorder.setMeta("sycl_work_group_size", gpuPaperGroupMetaString());
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("end_time", std::to_string(end_time));
    recorder.setMeta("DL", std::to_string(DL));
    recorder.setMeta("DH", std::to_string(DH));
    recorder.setMeta("DW", std::to_string(DW));
    recorder.setMeta("c_f", std::to_string(c_f));
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");
    recorder.setMeta("geometry_note", "CPC_appendix_matched_DL5.366_DH2_DW0.5_probe_y0.2");
    recorder.setMeta("water_wall_bounding", "none_single_phase");
    recorder.setMeta("dt_guard", "min_dt_and_max_acoustic_substeps_and_nan");
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "no_explicit_transfer_API_in_this_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected",
                     "detail deferred to profile short window policy");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected",
                                         "only_on_B2_B3_profile_window");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "needs_valid_telemetry_coverage");
    recorder.registerMissing({"CPC_digitized_reference_error", "not_collected",
                              "CPC fig10.1 curves not yet bound into this run",
                              "case_b2", "post-process probes vs digitized CPC after formal runs"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        if (mode == "benchmark")
        {
            for (int j = 0; j <= 50; ++j)
                probe_times.push_back(0.2 * double(j));
        }
        else
        {
            for (int j = 0; j < 100; ++j)
                probe_times.push_back(double(j) * double(end_time) / 99.0);
        }
        recorder.registerMissing({"probe_schedule_json_unreadable", "not_collected",
                                  "failed to parse schedule; fell back to mode default grid",
                                  schedule_path, "fix GPU_PAPER_SCHEDULE_JSON path"});
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;

    const StdVec<Real> probe_y = {0.01, 0.1, 0.2, 0.24, 0.252, 0.266};
    int state_fail = 0; // 1=dt, 2=nan, 3=acoustic_substep_cap

    auto sample_probes = [&](double t, double requested)
    {
        fluid_observer_contact_relation.exec();
        fluid_observer_pressure.writeToFile(0);
        record_water_mechanical_energy.writeToFile(0);
        auto p_obs = fluid_observer.getBaseParticles().getVariableDataByName<Real>("Pressure");
        Real mech_e = 0.0;
        {
            auto &bp = water_block.getBaseParticles();
            size_t n = bp.TotalRealParticles();
            auto pos = bp.getVariableDataByName<Vecd>("Position");
            auto vel = bp.getVariableDataByName<Vecd>("Velocity");
            auto mass = bp.getVariableDataByName<Real>("Mass");
            auto rho = bp.getVariableDataByName<Real>("Density");
            Real mass_sum = 0, x_front = -1e30, vmax = 0, rho_dev = 0, ke = 0, pe = 0;
            for (size_t i = 0; i < n; ++i)
            {
                mass_sum += mass[i];
                x_front = std::max(x_front, Real(pos[i][0]));
                vmax = std::max(vmax, Real(vel[i].norm()));
                rho_dev = std::max(rho_dev, std::abs(Real(rho[i] / rho0_f) - Real(1.0)));
                ke += Real(0.5) * mass[i] * Real(vel[i].squaredNorm());
                pe += mass[i] * gravity_g * Real(pos[i][1]);
            }
            if (n == 0)
                x_front = 0;
            mech_e = ke + pe;
            Real p_cpc = p_obs[2];

            auto emit = [&](const std::string &pid, const std::string &qid, Real value,
                            const std::string &unit, Real x, Real y, Real z)
            {
                ProbeSample ps;
                ps.physical_time_s = t;
                ps.requested_time_s = requested;
                ps.probe_id = pid;
                ps.body = "WaterBody";
                ps.quantity = qid;
                ps.value = double(value);
                ps.unit = unit;
                ps.x = double(x);
                ps.y = double(y);
                ps.z = double(z);
                recorder.addProbe(ps);
            };
            for (size_t k = 0; k < probe_y.size() && k < size_t(fluid_observer.getBaseParticles().TotalRealParticles()); ++k)
            {
                emit("p_y" + std::to_string(probe_y[k]), "pressure", p_obs[k], "Pa_code", DL, probe_y[k], 0.5 * DW);
            }
            emit("water", "mechanical_energy", mech_e, "J_code", 0, 0, 0);
            emit("water", "mass", mass_sum, "kg_code", 0, 0, 0);
            emit("water", "x_front", x_front, "m", x_front, 0, 0);
            emit("water", "vmax", vmax, "m/s", 0, 0, 0);
            emit("water", "rho_dev_max", rho_dev, "1", 0, 0, 0);

            recorder.addDiagnostic(t, "WaterBody", "N_active", double(n), "1", "measured");
            recorder.addDiagnostic(t, "WaterBody", "pressure_y0p2", double(p_cpc), "Pa_code", "measured");
            recorder.addDiagnostic(t, "WaterBody", "x_front", double(x_front), "m", "measured");
            recorder.addDiagnostic(t, "WaterBody", "mechanical_energy", double(mech_e), "J_code", "measured");
            recorder.addDiagnostic(t, "WaterBody", "vmax", double(vmax), "m/s", "measured");
            recorder.addDiagnostic(t, "WaterBody", "rho_dev_max", double(rho_dev), "1", "measured");

            if (!std::isfinite(double(mech_e)) || !std::isfinite(double(vmax)) || !std::isfinite(double(rho_dev)) ||
                !std::isfinite(double(p_cpc)) || !std::isfinite(double(mass_sum)))
            {
                std::cout << "ERROR: non-finite probe state at Time=" << t
                          << " mech_e=" << mech_e << " vmax=" << vmax
                          << " rho_dev=" << rho_dev << " p=" << p_cpc << std::endl;
                state_fail = 2;
            }
        }
    };

    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();

    size_t number_of_iterations = 0;
    int screen_output_interval = 100;

    auto wall0 = std::chrono::steady_clock::now();
    TimeInterval interval_output;
    TimeInterval interval_acoustic;
    TimeInterval interval_advection;
    TimeInterval interval_config;

    if (enable_vtp)
        body_states_recording.writeToFile();
    {
        double req0 = probe_times.empty() ? 0.0 : probe_times[0];
        sample_probes(sv_physical_time->getValue(), req0);
        if (!probe_times.empty())
            probe_cursor = 1;
        if (enable_vtp && !field_times.empty() && field_times[0] <= 0.0)
            field_cursor = 1;
    }

    TickCount t1 = TickCount::now();
    const Real dt_min = Real(1.0e-10);
    const int max_acoustic_substeps = 20000;

    while (sv_physical_time->getValue() < end_time && state_fail == 0)
    {
        TickCount ta = TickCount::now();
        fluid_density_summation.exec();
        fluid_density_regularization.exec();
        water_advection_step_setup.exec();
        Real advection_dt = fluid_advection_time_step.exec();
        if (!(advection_dt > dt_min) || !std::isfinite(double(advection_dt)))
        {
            std::cout << "ERROR: bad advection_dt=" << advection_dt
                      << " at N=" << number_of_iterations
                      << " Time=" << sv_physical_time->getValue() << std::endl;
            state_fail = 1;
            break;
        }
        fluid_boundary_indicator.exec();
        fluid_linear_correction_matrix.exec();
        interval_advection += TickCount::now() - ta;

        Real relaxation_time = 0.0;
        Real acoustic_dt = 0.0;
        int acoustic_substeps = 0;
        while (relaxation_time < advection_dt)
        {
            ta = TickCount::now();
            acoustic_dt = fluid_acoustic_time_step.exec();
            if (!(acoustic_dt > dt_min) || !std::isfinite(double(acoustic_dt)))
            {
                std::cout << "ERROR: bad acoustic_dt=" << acoustic_dt
                          << " at N=" << number_of_iterations
                          << " Time=" << sv_physical_time->getValue()
                          << " substeps=" << acoustic_substeps << std::endl;
                state_fail = 1;
                break;
            }
            if (++acoustic_substeps > max_acoustic_substeps)
            {
                std::cout << "ERROR: acoustic substep cap " << max_acoustic_substeps
                          << " exceeded at N=" << number_of_iterations
                          << " Time=" << sv_physical_time->getValue()
                          << " advection_dt=" << advection_dt
                          << " acoustic_dt=" << acoustic_dt << std::endl;
                state_fail = 3;
                break;
            }
            fluid_acoustic_step_1st_half.exec(acoustic_dt);
            fluid_acoustic_step_2nd_half.exec(acoustic_dt);
            relaxation_time += acoustic_dt;
            sv_physical_time->incrementValue(acoustic_dt);
            interval_acoustic += TickCount::now() - ta;

            ta = TickCount::now();
            double t_now = sv_physical_time->getValue();
            double requested = 0.0;
            if (dueForSample(t_now, probe_times, probe_cursor, &requested))
            {
                sample_probes(t_now, requested);
                if (state_fail)
                    break;
            }
            if (enable_vtp && dueForSample(t_now, field_times, field_cursor, &requested))
                body_states_recording.writeToFile();
            interval_output += TickCount::now() - ta;
        }
        if (state_fail)
            break;

        water_update_particle_position.exec();

        if (number_of_iterations % screen_output_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                      << sv_physical_time->getValue()
                      << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt
                      << "	acoustic_substeps = " << acoustic_substeps << "\n";
        }
        number_of_iterations++;

        ta = TickCount::now();
        if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            particle_sort.exec();
        water_cell_linked_list.exec();
        water_block_update_complex_relation.exec();
        interval_config += TickCount::now() - ta;
    }

    if (state_fail == 0)
        sample_probes(sv_physical_time->getValue(), double(end_time));

    TimeInterval tt = TickCount::now() - t1 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_advection = "
              << interval_advection.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic = "
              << interval_acoustic.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_config = "
              << interval_config.seconds() << "\n";

    recorder.stages().add("advection_physics", interval_advection.seconds());
    recorder.stages().add("acoustic_step", interval_acoustic.seconds());
    recorder.stages().add("configuration", interval_config.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("advection", int(number_of_iterations), 0, 0, 0, 0, 0, 0, "fluid_ck", "measured");
    recorder.addBodyStats(sv_physical_time->getValue(), "WaterBody",
                          int(water_block.getBaseParticles().TotalRealParticles()),
                          int(water_block.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    const char *status = "ok";
    if (state_fail == 1)
        status = "timestep_collapse";
    else if (state_fail == 2)
        status = "non_finite_state";
    else if (state_fail == 3)
        status = "acoustic_substep_cap";
    recorder.finalize(status, 0.0, sv_physical_time->getValue(), wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval_output.seconds(),
                      0.0, interval_output.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return state_fail == 0 ? 0 : 1;
}
