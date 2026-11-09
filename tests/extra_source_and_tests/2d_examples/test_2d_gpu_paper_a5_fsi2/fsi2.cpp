/**
 * @file fsi2.cpp
 * @brief GPU paper A5: cylinder + flexible beam FSI2 — SYCL/CK path.
 * @details Duration V2 End_Time=200. BC adapted from classic periodic+sponge to
 *          emitter/disposer (SYCL fish/cylinder pattern). Solid: SVK + FixBodyPartConstraintCK.
 *          Insert body: relax → reload/Reload.xml (classic FSI2 workflow), then physics.
 *          No StructureNumericalDamping (matches classic FSI2).
 */
#include "fsi2.h"
#include "gpu_paper_recording.h"
#include "sphinxsys.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>

using namespace SPH;
using namespace sphinxsys_gpu_paper;

int main(int ac, char *av[])
{
    applyGpuPaperResolutionFromEnv();

    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    // Classic FSI2 defaults: physics starts from relaxed InsertedBody reload.
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    if (const char *v = std::getenv("GPU_PAPER_RELAX"))
        sph_system.setRunParticleRelaxation(std::string(v) == "1");
    if (const char *v = std::getenv("GPU_PAPER_RELOAD"))
        sph_system.setReloadParticles(std::string(v) == "1");
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av); // --relax / --reload override env
#endif

    //----------------------------------------------------------------------
    // Particle relaxation for body-fitted InsertedBody (cylinder + beam).
    // Writes ./reload/Reload.xml then exits (two-step workflow like classic FSI2).
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        RelaxationSystem relaxation_system(system_domain_bounds, global_resolution);
        auto &insert_shape = relaxation_system.addShape<Insert>("InsertedBody");
        auto &insert_relax = relaxation_system.addBody<SolidBody>(insert_shape);
        insert_relax.defineAdaptationRatios(1.15, 2.0);
        LevelSetShape &insert_level_set = insert_relax.defineBodyLevelSetShape().writeLevelSet();
        insert_relax.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_insert_surface(insert_relax);
        Inner<> insert_relax_inner(insert_relax);

        SPHSolver relax_solver(relaxation_system);
        auto &main_methods = relax_solver.getMainMethodContainer();
        auto &host_methods = relax_solver.getHostMethodContainer();

        auto &insert_cell_linked_list = main_methods.addCellLinkedListDynamics(insert_relax);
        auto &insert_update_inner = main_methods.addRelationDynamics(insert_relax_inner);
        auto &random_insert_particles =
            host_methods.addStateDynamics<RandomizeParticlePositionCK>(insert_relax, Real(0.25));
        auto &relaxation_residual =
            main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(insert_relax_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(insert_relax, insert_level_set);
        auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(insert_relax);
        auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(insert_relax);
        auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_insert_surface);
        auto &body_state_recorder =
            main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(relaxation_system);
        auto &write_particle_reload_files =
            main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&insert_relax});

        random_insert_particles.exec();
        level_set_bounding.exec();
        body_state_recorder.writeToFile(0);

        int ite_p = 0;
        const int relax_steps = 1000;
        while (ite_p < relax_steps)
        {
            insert_cell_linked_list.exec();
            insert_update_inner.exec();
            relaxation_residual.exec();
            Real relaxation_step = relaxation_scaling.exec();
            update_particle_position.exec(relaxation_step);
            level_set_bounding.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        write_particle_reload_files.writeToFile();
        std::cout << "Wrote reload file under ./reload/Reload.xml (cwd="
                  << "bin working directory)" << std::endl;
        return 0;
    }

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody insert_body(sph_system, makeShared<Insert>("InsertedBody"));
    insert_body.defineAdaptationRatios(1.15, 2.0);
    insert_body.defineBodyLevelSetShape().writeLevelSet();
    insert_body.defineMatterMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    if (sph_system.ReloadParticles())
        insert_body.generateParticles<BaseParticles, Reload>(insert_body.Name());
    else
        insert_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody beam_observer(sph_system, "BeamObserver");
    StdVec<Vecd> beam_observation_location = {0.5 * (BRT + BRB)};
    beam_observer.defineAdaptationRatios(1.15, 2.0);
    beam_observer.generateParticles<ObserverParticles>(beam_observation_location);

    BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createBeamBaseShape()));

    Inner<> water_block_inner(water_block);
    Inner<> insert_body_inner(insert_body, ConfigType::Lagrangian);
    Contact<> fluid_wall_contact(water_block, wall_boundary);
    Contact<> fluid_structure_contact(water_block, insert_body);
    Contact<> insert_body_contact(insert_body, water_block);
    Contact<> beam_observer_contact(beam_observer, insert_body, ConfigType::Lagrangian);

    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.getHostMethodContainer();
    auto &main_methods = sph_solver.getMainMethodContainer();

    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary).exec();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(insert_body).exec();

    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &insert_cell_linked_list = main_methods.addCellLinkedListDynamics(insert_body);
    auto &water_update_complex =
        main_methods.addRelationDynamics(water_block_inner, fluid_wall_contact, fluid_structure_contact);
    auto &insert_inner_build = main_methods.addRelationDynamics(insert_body_inner);
    auto &insert_contact_relation = main_methods.addRelationDynamics(insert_body_contact);
    auto &observer_contact_relation = main_methods.addRelationDynamics(beam_observer_contact);

    auto &insert_corrected_configuration =
        main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(insert_body_inner);

    auto &initialize_displacement =
        main_methods.addStateDynamics<FSI::InitializeDisplacementCK>(insert_body);
    auto &update_average_velocity =
        main_methods.addStateDynamics<FSI::UpdateAverageVelocityAndAccelerationCK>(insert_body);

    auto &insert_stress_1st =
        main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration1stHalfPK2, SaintVenantKirchhoffSolid>(insert_body_inner);
    auto &constraint_beam_base =
        main_methods.addStateDynamics<FixBodyPartConstraintCK>(beam_base);
    auto &insert_stress_2nd =
        main_methods.addInteractionDynamicsOneLevel<
            solid_dynamics::StructureIntegration2ndHalf>(insert_body_inner);
    auto &insert_acoustic_dt =
        main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(insert_body);
    auto &insert_update_normal =
        main_methods.addStateDynamics<solid_dynamics::UpdateElasticNormalDirectionCK>(insert_body);

    auto &water_advection_step_setup =
        main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position =
        main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);

    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_wall_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_structure_contact);

    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
                        fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_wall_contact)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(fluid_structure_contact)
            .addGeneralPostInteraction<FSI::PressureForceFromFluid, WithUpdate, AcousticRiemannSolverCK, NoKernelCorrectionCK>(
                insert_body_contact);

    // FreeStream density needs Indicator from FreeSurfaceIndicationCK.
    auto &fluid_boundary_indicator =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::FreeSurfaceIndicationCK>(water_block_inner)
            .addPostContactInteraction(fluid_wall_contact)
            .addPostContactInteraction(fluid_structure_contact);

    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(water_block_inner)
            .addPostContactInteraction(fluid_wall_contact)
            .addPostContactInteraction(fluid_structure_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeStream>(water_block);

    OrientedBoxByCell emitter_buffer_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));

    auto &transport_correction =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(fluid_wall_contact)
            .addPostContactInteraction<Boundary, NoKernelCorrectionCK>(fluid_structure_contact)
            .addPostStateDynamics<fluid_dynamics::TransportVelocityCorrectionCK, NoLimiter>(
                water_block, Real(0.25))
            .addPostStateDynamics<ConstantConstraintCK, Vecd>(emitter_buffer_part, "Displacement", Vecd::Zero());

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<
                        fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(fluid_wall_contact)
            .addPostContactInteraction<Wall, Viscosity, NoKernelCorrectionCK>(fluid_structure_contact)
            .addGeneralPostInteraction<FSI::ViscousForceFromFluid, WithUpdate, Viscosity, NoKernelCorrectionCK>(
                insert_body_contact);

    auto &fluid_advection_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionViscousTimeStepCK>(water_block, U_f);
    auto &fluid_acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(water_block);

    OrientedBoxByParticle emitter_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);
    auto &emitter_injection =
        main_methods.addStateDynamics<fluid_dynamics::EmitterInflowInjectionCK>(emitter_part);
    auto &inflow_condition =
        main_methods.addStateDynamics<fluid_dynamics::EmitterInflowConditionCK, ParabolicInflowSpeed>(
            emitter_buffer_part, ParabolicInflowSpeed(U_f, Real(2.0), Real(0.5) * DH));

    Vec2d correct_disposer_translation = Vec2d(DL, -0.25 * DH) + disposer_halfsize;
    OrientedBoxByCell disposer_part(
        water_block,
        OrientedBox(xAxis, Transform(Vec2d(correct_disposer_translation)), disposer_halfsize));
    auto &disposer_indication =
        main_methods.addStateDynamics<fluid_dynamics::WithinDisposerIndication>(disposer_part);
    auto &particle_deletion =
        main_methods.addStateDynamics<fluid_dynamics::OutflowParticleDeletion>(water_block);

    auto &write_real_body_states =
        main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    ObservedQuantityRecording<MainExecutionPolicy, Vecd>
        write_beam_tip_displacement(beam_observer_contact, "Position");

    //----------------------------------------------------------------------
    wall_cell_linked_list.exec();
    insert_cell_linked_list.exec();
    water_cell_linked_list.exec();
    insert_inner_build.exec();
    water_update_complex.exec();
    insert_contact_relation.exec();
    observer_contact_relation.exec();
    insert_corrected_configuration.exec();

    fluid_boundary_indicator.exec();
    fluid_density_regularization.exec();
    water_block.getBaseParticles().addEvolvingVariable<Vecd>("Force");
    water_advection_step_setup.exec();
    transport_correction.exec();
    fluid_viscous_force.exec();

    const char *end_env = std::getenv("GPU_PAPER_END_TIME");
    Real end_time = end_env ? Real(std::atof(end_env)) : Real(200.0);
    int screen_output_interval = 100;

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "A5_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "validate";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/A5_fsi2/schedules/A5_validate_T200_bound.json";
    bool enable_vtp = false;
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "A5", mode);
    recorder.setMeta("case_name", "2d_gpu_paper_a5_fsi2");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("global_resolution", std::to_string(global_resolution));
    recorder.setMeta("D_over_dp", std::to_string(d_over_dp));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_ck_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("api", "SPHSolver_FSI_CK_emitter_disposer");
    recorder.setMeta("bc_note", "classic_periodic_sponge_replaced_by_emitter_disposer");
    recorder.setMeta("insert_reload", sph_system.ReloadParticles() ? "1" : "0");
    recorder.setMeta("end_time", std::to_string(end_time));
    recorder.setMeta("schedule_json", schedule_path);
    recorder.setMeta("enable_vtp", enable_vtp ? "1" : "0");

    // Alignment knobs (CPU/GPU tip drift is dominated by emitter/disposer + FP):
    // - disable particle sort to avoid extra neighbor-order divergence
    // - optional fixed acoustic dt so both ends share the same time grid
    bool disable_particle_sort = false;
    if (const char *v = std::getenv("GPU_PAPER_DISABLE_PARTICLE_SORT"))
        disable_particle_sort = (std::string(v) == "1");
    Real fixed_acoustic_dt = Real(-1.0);
    if (const char *v = std::getenv("GPU_PAPER_FIXED_ACOUSTIC_DT"))
        fixed_acoustic_dt = Real(std::atof(v));
    recorder.setMeta("disable_particle_sort", disable_particle_sort ? "1" : "0");
    recorder.setMeta("fixed_acoustic_dt", fixed_acoustic_dt > Real(0.0) ? std::to_string(fixed_acoustic_dt) : "adaptive");
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected", "ck_device_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected", "");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected", "");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "");
    recorder.registerMissing({"period_amplitude_validated", "not_collected",
                              "175-200 window period must be checked after full T=200",
                              "case_a5", "post-process tip y(t)"});
    recorder.registerMissing({"classic_periodic_bc", "adapted",
                              "SYCL port uses emitter/disposer instead of periodic+sponge",
                              "case_a5", "document in paper methods"});
    if (!sph_system.ReloadParticles())
    {
        recorder.registerMissing({"insert_body_reload", "disabled",
                                  "running Lattice InsertedBody; prefer relax then reload",
                                  "case_a5", "run with GPU_PAPER_RELAX=1 then GPU_PAPER_RELOAD=1"});
    }

    std::vector<double> probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    std::vector<double> field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    if (probe_times.empty())
    {
        for (int j = 0; j < 51; ++j)
            probe_times.push_back(double(j) * double(end_time) / 50.0);
    }
    size_t probe_cursor = 0;
    size_t field_cursor = 0;
    Vecd tip0 = 0.5 * (BRT + BRB);
    // Keep last good tip: ObservingQuantityCK("Position") overwrites observer Position;
    // a zero sample is almost always empty-neighbor interpolation, not physical motion.
    Vecd tip_last = tip0;
    int tip_reject_count = 0;

    auto sample_probes = [&](double t, double requested)
    {
        // Total-Lagrangian observer contact must stay on the initial neighbor list.
        // Rebuilding after insert_cell_linked_list (current Position CLL) makes
        // InitialPosition queries miss neighbors → interpolated tip becomes (0,0).
        write_beam_tip_displacement.writeToFile(0);
        const Vecd tip_raw = write_beam_tip_displacement.getObservedQuantity()[0];
        const Vecd tip_disp = tip_raw - tip0;
        const Real tip_disp2 = tip_disp.squaredNorm();
        // Physical tip motion is O(0.1); a failed interpolate yields ~tip0 (e.g. (0,0) → disp²≈36).
        const bool tip_ok = std::isfinite(double(tip_raw[0])) && std::isfinite(double(tip_raw[1])) &&
                            tip_disp2 < Real(1.0);
        Vecd tip = tip_ok ? tip_raw : tip_last;
        if (!tip_ok)
            tip_reject_count++;
        else
            tip_last = tip_raw;
        // Keep probe station at tip0 so Position is not left at a failed (0,0) sample.
        beam_observer.getBaseParticles().getVariableDataByName<Vecd>("Position")[0] = tip0;
        Real x = tip[0], y = tip[1];
        auto emit = [&](const std::string &qid, Real value)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = "beam_tip";
            ps.body = "InsertedBody";
            ps.quantity = qid;
            ps.value = double(value);
            ps.unit = "m";
            ps.x = double(tip0[0]);
            ps.y = double(tip0[1]);
            recorder.addProbe(ps);
        };
        emit("tip_x", x);
        emit("tip_y", y);
        emit("tip_disp_x", x - tip0[0]);
        emit("tip_disp_y", y - tip0[1]);
        recorder.addDiagnostic(t, "InsertedBody", "tip_y", double(y), "m", tip_ok ? "measured" : "held_last_good");
        recorder.addDiagnostic(t, "WaterBody", "N_active",
                               double(water_block.getBaseParticles().TotalRealParticles()), "1", "measured");
    };

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_trigger = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t number_of_iterations = 0;
    int dt_fail = 0;

    auto wall0 = std::chrono::steady_clock::now();
    TickCount t1 = TickCount::now();
    TimeInterval interval_output;

    if (enable_vtp)
        write_real_body_states.writeToFile();
    {
        double req0 = probe_times.empty() ? 0.0 : probe_times[0];
        sample_probes(time_stepper.getPhysicalTime(), req0);
        if (!probe_times.empty())
            probe_cursor = 1;
    }

    while (!time_stepper.isEndTime(end_time) && dt_fail == 0)
    {
        Real acoustic_dt = fixed_acoustic_dt > Real(0.0)
                               ? time_stepper.incrementPhysicalTime(fixed_acoustic_dt)
                               : time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        if (!(acoustic_dt > Real(1.0e-12)) || !std::isfinite(double(acoustic_dt)))
        {
            std::cout << "ERROR: bad acoustic_dt=" << acoustic_dt
                      << " at N=" << number_of_iterations << std::endl;
            dt_fail = 1;
            break;
        }

        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        inflow_condition.exec();
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);

        initialize_displacement.exec();
        Real dt_s_sum = 0.0;
        int guard = 0;
        int solid_sub_count = 0;
        while (dt_s_sum < acoustic_dt && guard++ < 10000)
        {
            Real dt_s_candidate = insert_acoustic_dt.exec();
            if (!(dt_s_candidate > Real(0.0)) || !std::isfinite(double(dt_s_candidate)))
            {
                std::cout << "ERROR: bad solid dt_s=" << dt_s_candidate << std::endl;
                dt_fail = 1;
                break;
            }
            Real dt_s = SMIN(dt_s_candidate, acoustic_dt - dt_s_sum);
            insert_stress_1st.exec(dt_s);
            constraint_beam_base.exec();
            insert_stress_2nd.exec(dt_s);
            dt_s_sum += dt_s;
            solid_sub_count++;
        }
        if (dt_fail)
            break;
        update_average_velocity.exec(acoustic_dt);

        double physical_time = time_stepper.getPhysicalTime();
        double requested = 0.0;
        TickCount to = TickCount::now();
        if (dueForSample(physical_time, probe_times, probe_cursor, &requested))
            sample_probes(physical_time, requested);
        if (enable_vtp && dueForSample(physical_time, field_times, field_cursor, &requested))
            write_real_body_states.writeToFile();
        interval_output += TickCount::now() - to;

        if (advection_trigger(fluid_advection_time_step))
        {
            number_of_iterations++;
            Real Dt = advection_trigger.getInterval();
            water_update_particle_position.exec();
            insert_update_normal.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "	Time = " << physical_time
                          << "	Dt = " << Dt
                          << "	dt_s_sub = " << solid_sub_count << "\n";
            }

            emitter_injection.exec();
            disposer_indication.exec();
            particle_deletion.exec();
            if (!disable_particle_sort && number_of_iterations % 100 == 0)
                particle_sort.exec();

            insert_cell_linked_list.exec();
            water_cell_linked_list.exec();
            water_update_complex.exec();
            insert_contact_relation.exec();

            fluid_boundary_indicator.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
        }
    }

    double physical_time = time_stepper.getPhysicalTime();
    if (!dt_fail)
        sample_probes(physical_time, double(end_time));

    TimeInterval tt = TickCount::now() - t1 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    recorder.stages().add("solver_exclusive", tt.seconds());
    recorder.stages().add("vtp_or_io", interval_output.seconds());
    recorder.addStepBlock("advection", int(number_of_iterations), 0, 0, 0, 0, 0, 0, "fsi_sycl_ck", "measured");
    recorder.addBodyStats(physical_time, "WaterBody",
                          int(water_block.getBaseParticles().TotalRealParticles()),
                          int(water_block.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");
    recorder.addBodyStats(physical_time, "InsertedBody",
                          int(insert_body.getBaseParticles().TotalRealParticles()),
                          int(insert_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    recorder.setMeta("tip_reject_count", std::to_string(tip_reject_count));
    if (tip_reject_count > 0)
    {
        recorder.registerMissing({"tip_interpolate_reject", "held_last_good",
                                  "some tip samples rejected (near-zero/NaN interpolate); held last good tip",
                                  "case_a5", "should be 0 after Lagrangian observer fix"});
        std::cout << "[gpu_paper] tip_reject_count=" << tip_reject_count << std::endl;
    }

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    const char *status = dt_fail ? "timestep_collapse" : "ok";
    recorder.finalize(status, 0.0, physical_time, wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval_output.seconds(),
                      0.0, interval_output.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return dt_fail;
}
