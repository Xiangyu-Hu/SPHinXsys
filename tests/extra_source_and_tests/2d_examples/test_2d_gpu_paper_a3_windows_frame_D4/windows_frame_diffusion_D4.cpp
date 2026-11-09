/**
 * @file  windows_frame_diffusion_D4.cpp
 * @brief GPU paper A3: ISO 10077-2 D4 windows-frame diffusion — SYCL/CK path.
 * @details Duration V2 End_Time=0.02 (pseudo-time). Uses DiffusionRelaxationCK + Robin CK.
 */
#include "windows_frame_diffusion_D4.h"
#include "sphinxsys.h"
#include "gpu_paper_local_isotropic_diffusion.h"
#include "gpu_paper_diffusion_robin.h"
#include "gpu_paper_recording.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>

using namespace SPH;
using namespace sphinxsys_gpu_paper;

void applyResolutionFromEnv()
{
    Real scale = 1.0;
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        scale = Real(std::atof(v));
    if (scale < Real(1.0))
        scale = Real(1.0);
    global_resolution = Real(0.001) / scale;
    BW = global_resolution * Real(2.0);
    system_domain_bounds = BoundingBoxd(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
}

int main(int ac, char *av[])
{
    applyResolutionFromEnv();

    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    SolidBody diffusion_body(sph_system, makeShared<MultiPolygonShape>(createOverallStructureBody(), "DiffusionBody"));
    diffusion_body.defineMatterMaterial<Solid>();
    diffusion_body.addMaterialProperty<GpuPaperLocalIsotropicDiffusion>(species_name, wood_cond, epdm_cond);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody boundary_Robin_in(sph_system, makeShared<MultiPolygonShape>(createInternalAirBody(), "InternalConvectionBoundary"));
    boundary_Robin_in.defineMatterMaterial<Solid>();
    boundary_Robin_in.generateParticles<BaseParticles, Lattice>();

    SolidBody boundary_Robin_ex(sph_system, makeShared<MultiPolygonShape>(createExternalAirBody(), "ExternalConvectionBoundary"));
    boundary_Robin_ex.defineMatterMaterial<Solid>();
    boundary_Robin_ex.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    StdVec<Vecd> obs_pts = createObservationPoints();
    temperature_observer.generateParticles<ObserverParticles>(obs_pts);

    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> contact_robin_in(diffusion_body, boundary_Robin_in);
    Contact<> contact_robin_ex(diffusion_body, boundary_Robin_ex);
    Contact<> temperature_observer_contact(temperature_observer, diffusion_body);

    BodyRegionByParticle epdm_body(diffusion_body, makeShared<MultiPolygonShape>(createEPDMBody()));
    BodyRegionByParticle panel_body(diffusion_body, makeShared<MultiPolygonShape>(createPanelBody()));
    BodyRegionByParticle ac_body1(diffusion_body, makeShared<MultiPolygonShape>(createACBody1()));
    BodyRegionByParticle ac_body2(diffusion_body, makeShared<MultiPolygonShape>(createACBody2()));
    BodyRegionByParticle ac_open_body1(diffusion_body, makeShared<MultiPolygonShape>(createACOpenBody1()));
    BodyRegionByParticle decreased_convection_body(
        boundary_Robin_in, makeShared<MultiPolygonShape>(createDecreasedInternalConvectionBody()));

    // Host one-shot geometry / material assignment (classic helpers).
    SimpleDynamics<NormalDirectionFromBodyShape> diffusion_body_normal_direction(diffusion_body);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_in(boundary_Robin_in);
    SimpleDynamics<NormalDirectionFromBodyShape> Robin_normal_direction_ex(boundary_Robin_ex);
    SimpleDynamics<LocalDiffusivityDefinition> epdm_diffusivity(epdm_body, epdm_cond);
    SimpleDynamics<LocalDiffusivityDefinition> panel_diffusivity(panel_body, pane_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_diffusivity(ac_body1, ac1_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac2_diffusivity(ac_body2, ac2_cond);
    SimpleDynamics<LocalDiffusivityDefinition> ac1_open_diffusivity(ac_open_body1, ac1_open_cond);
    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> diffusion_cell_linked_list(diffusion_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> robin_in_cell_linked_list(boundary_Robin_in);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> robin_ex_cell_linked_list(boundary_Robin_ex);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>, Contact<>>
        diffusion_update_relation(diffusion_body_inner, contact_robin_in, contact_robin_ex);
    UpdateRelation<MainExecutionPolicy, Contact<>> observer_contact_relation(temperature_observer_contact);

    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, GpuPaperLocalIsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Robin<GpuPaperLocalIsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Robin<GpuPaperLocalIsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, GpuPaperLocalIsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Robin<GpuPaperLocalIsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Robin<GpuPaperLocalIsotropicDiffusion>, NoKernelCorrectionCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(diffusion_body_inner, &diffusion_body.getMaterialProperty<GpuPaperLocalIsotropicDiffusion>()),
                                 DynamicsArgs(contact_robin_in, &diffusion_body.getMaterialProperty<GpuPaperLocalIsotropicDiffusion>()),
                                 DynamicsArgs(contact_robin_ex, &diffusion_body.getMaterialProperty<GpuPaperLocalIsotropicDiffusion>()));

    // Robin BC values on contact bodies (registered by Robin<> ctor above).
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_in(boundary_Robin_in);
    SimpleDynamics<RobinBoundaryDefinition> robin_boundary_condition_ex(boundary_Robin_ex);
    SimpleDynamics<LocalConvectionDefinition> decreased_convection_initial_condition(
        decreased_convection_body, convection_i_decreased);

    GetDiffusionTimeStepSize get_time_step_size(diffusion_body);

    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_states(sph_system);
    ObservedQuantityRecording<MainExecutionPolicy, Real, RestoringCorrection>
        write_solid_temperature(temperature_observer_contact, species_name);
    // Host reduce after device→host sync via ObservedQuantityRecording / variable access.
    ReduceDynamics<Average<QuantitySummation<Real, SPHBody>>> heat_in_sum(
        diffusion_body, species_name + "TransferFromInternalConvectionBoundary");
    ReduceDynamics<Average<QuantitySummation<Real, SPHBody>>> heat_ex_sum(
        diffusion_body, species_name + "TransferFromExternalConvectionBoundary");

    // Host init before device kernels.
    setup_diffusion_initial_condition.exec();
    epdm_diffusivity.exec();
    panel_diffusivity.exec();
    ac1_diffusivity.exec();
    ac2_diffusivity.exec();
    ac1_open_diffusivity.exec();
    diffusion_body_normal_direction.exec();
    Robin_normal_direction_in.exec();
    Robin_normal_direction_ex.exec();

    diffusion_cell_linked_list.exec();
    robin_in_cell_linked_list.exec();
    robin_ex_cell_linked_list.exec();
    diffusion_update_relation.exec();
    observer_contact_relation.exec();

    robin_boundary_condition_in.exec();
    robin_boundary_condition_ex.exec();
    decreased_convection_initial_condition.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    const char *end_env = std::getenv("GPU_PAPER_END_TIME");
    Real End_Time = end_env ? Real(std::atof(end_env)) : Real(0.02);
    Real Observe_time = Real(0.01) * End_Time;
    Real Output_Time = Real(0.1) * End_Time;
    Real dt = 0.0;

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "A3_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "validate";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/A3_windows_frame_D4/schedules/"
              "A3_validate_T0p02_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    Real res_scale = Real(0.001) / global_resolution;

    RunRecorder recorder(run_dir, run_id, "A3", mode);
    recorder.setMeta("case_name", "2d_gpu_paper_a3_windows_frame_D4");
    recorder.setMeta("variant", "D4");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("global_resolution", std::to_string(global_resolution));
    recorder.setMeta("resolution_scale", std::to_string(res_scale));
    recorder.setMeta("coordinate", "pseudo_time");
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_ck_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("api", "DiffusionRelaxationCK_Robin_GpuPaperLocalIsotropicDiffusion");
    recorder.setMeta("end_time", std::to_string(End_Time));
    recorder.setMeta("schedule_json", schedule_path);
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "ck_device_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected", "");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected", "");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "");
    recorder.registerMissing({"COMSOL_reference_L2D_Uf", "not_collected",
                              "reference table not bound into this run",
                              "case_a3", "post-process vs COMSOL/ISO after formal runs"});

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        for (int j = 0; j < 100; ++j)
            probe_times.push_back(double(j) * double(End_Time) / 99.0);
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;
    std::vector<Real> prev_T(obs_pts.size(), initial_temperature);

    auto sample_probes = [&](double t, double requested)
    {
        observer_contact_relation.exec();
        write_solid_temperature.writeToFile(ite);
        auto *Tobs = temperature_observer.getBaseParticles().getVariableDataByName<Real>(species_name);
        Real q_in = heat_in_sum.exec();
        Real q_ex = heat_ex_sum.exec();
        Real res_l2 = 0.0;
        Real res_linf = 0.0;
        for (size_t i = 0; i < obs_pts.size(); ++i)
        {
            Real Ti = Tobs[i];
            Real dT = Ti - prev_T[i];
            res_l2 += dT * dT;
            res_linf = SMAX(res_linf, std::abs(dT));
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = "T_obs_" + std::to_string(i);
            ps.body = "DiffusionBody";
            ps.quantity = "Phi";
            ps.value = double(Ti);
            ps.unit = "C";
            ps.x = double(obs_pts[i][0]);
            ps.y = double(obs_pts[i][1]);
            recorder.addProbe(ps);
            prev_T[i] = Ti;
        }
        res_l2 = std::sqrt(res_l2 / Real(obs_pts.size()));

        auto emit = [&](const std::string &qid, Real value, const std::string &unit)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = "domain";
            ps.body = "DiffusionBody";
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = unit;
            recorder.addProbe(ps);
        };
        emit("heat_transfer_internal", q_in, "W_per_m_proxy");
        emit("heat_transfer_external", q_ex, "W_per_m_proxy");
        emit("heat_balance_residual", q_in + q_ex, "W_per_m_proxy");
        emit("temp_change_L2", res_l2, "C");
        emit("temp_change_Linf", res_linf, "C");
        recorder.addDiagnostic(t, "DiffusionBody", "N_active",
                               double(diffusion_body.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "DiffusionBody", "heat_balance_residual", double(q_in + q_ex), "W_per_m_proxy",
                               "measured");
        recorder.addDiagnostic(t, "DiffusionBody", "temp_change_L2", double(res_l2), "C", "measured");
    };

    auto wall0 = std::chrono::steady_clock::now();
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;

    dt = get_time_step_size.exec();
    if (enable_vtp)
        write_states.writeToFile();
    {
        double req0 = probe_times.empty() ? 0.0 : probe_times[0];
        sample_probes(physical_time, req0);
        if (!probe_times.empty())
            probe_cursor = 1;
        if (enable_vtp && !field_times.empty() && field_times[0] <= 0.0)
            field_cursor = 1;
    }

    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
                if (ite % 500 == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << ite << " Time: "
                              << physical_time << "	dt: " << dt << "\n";
                }

                diffusion_relaxation_rk2.exec(dt);
                ite++;
                dt = get_time_step_size.exec();
                if (!(dt > Real(0.0)) || !std::isfinite(double(dt)))
                {
                    std::cout << "ERROR: bad dt=" << dt << " at N=" << ite << " Time=" << physical_time << std::endl;
                    recorder.finalize("timestep_collapse", 0.0, physical_time,
                                      std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count(),
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                    return 1;
                }
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                double requested = 0.0;
                if (dueForSample(physical_time, probe_times, probe_cursor, &requested))
                    sample_probes(physical_time, requested);
                if (enable_vtp && dueForSample(physical_time, field_times, field_cursor, &requested))
                    write_states.writeToFile();
            }
        }

        TickCount t2 = TickCount::now();
        write_solid_temperature.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TickCount::interval_t tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << physical_time << " seconds." << std::endl;

    sample_probes(physical_time, double(End_Time));

    recorder.stages().add("diffusion_compute", tt.seconds());
    recorder.stages().add("vtp_or_io", interval.seconds());
    recorder.addStepBlock("diffusion", ite, 0, 0, 0, 0, 0, 0, "sycl_ck", "measured");
    recorder.addBodyStats(physical_time, "DiffusionBody",
                          int(diffusion_body.getBaseParticles().TotalRealParticles()),
                          int(diffusion_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    recorder.finalize("ok", 0.0, physical_time, wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval.seconds(),
                      0.0, interval.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return 0;
}
