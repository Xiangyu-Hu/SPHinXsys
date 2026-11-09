/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	GPU paper A1: 2D mixed Poiseuille (extra_source_and_tests copy).
 * @details Instrumented with V4 recording (probes/stages/summary). Do not use original tests/ path.
 * @note    GPU paper experiment (docs/gpu_paper_discussion).
 */

#include "gpu_paper_recording.h"
#include "sphinxsys.h" // SPHinXsys Library.
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
using namespace SPH;
using namespace sphinxsys_gpu_paper;

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.004;                    /**< Channel length. */
Real DH = 0.001;                    /**< Channel height. */
Real dh_over_dp = 20.0;             /**< L1=20, L2=40, L3=80 (Duration/V4). */
Real global_resolution = DH / dh_over_dp;
Real BW = global_resolution * 4;
StdVec<Vecd> observer_location;
BoundingBoxd system_domain_bounds(
    Vec2d(-2.0 * BW, -2.0 * BW),
    Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.1;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
Real c_f = 10.0 * U_f;

//----------------------------------------------------------------------
//  Geometric shapes for the channel and boundaries.
//----------------------------------------------------------------------
Real bidirectional_buffer_length = 5.0 * global_resolution;
Vec2d bidirectional_buffer_halfsize(
    0.5 * bidirectional_buffer_length, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d left_indicator_translation(0.0, 0.5 * DH);
Vec2d right_bidirectional_translation(
    DL - 0.5 * bidirectional_buffer_length, 0.5 * DH);
Vec2d right_indicator_translation(DL, 0.5 * DH);
Vec2d right_disposer_translation(
    DL - 0.5 * bidirectional_buffer_length, 0.5 * DH);
Vec2d normal(1.0, 0.0);

/** Apply GPU_PAPER_DH_OVER_DP before constructing bodies. */
void applyResolutionFromEnv()
{
    if (const char *v = std::getenv("GPU_PAPER_DH_OVER_DP"))
        dh_over_dp = Real(std::atof(v));
    if (dh_over_dp < Real(5.0))
        dh_over_dp = Real(20.0);
    global_resolution = DH / dh_over_dp;
    BW = global_resolution * 4;
    system_domain_bounds = BoundingBoxd(Vec2d(-2.0 * BW, -2.0 * BW), Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));
    bidirectional_buffer_length = 5.0 * global_resolution;
    bidirectional_buffer_halfsize = Vec2d(0.5 * bidirectional_buffer_length, 0.5 * DH);
    left_bidirectional_translation = bidirectional_buffer_halfsize;
    left_indicator_translation = Vec2d(0.0, 0.5 * DH);
    right_bidirectional_translation = Vec2d(DL - 0.5 * bidirectional_buffer_length, 0.5 * DH);
    right_indicator_translation = Vec2d(DL, 0.5 * DH);
    right_disposer_translation = Vec2d(DL - 0.5 * bidirectional_buffer_length, 0.5 * DH);
}
//----------------------------------------------------------------------
//  Inlet velocity profile for the left boundary (Poiseuille-like).
//----------------------------------------------------------------------
class InflowVelocityPrescribed : public VelocityPrescribed<WeaklyCompressibleFluid>
{
  public:
    InflowVelocityPrescribed(Real DH, Real U_f, Real mu_f)
        : VelocityPrescribed<WeaklyCompressibleFluid>(),
          DH_(DH), U_f_(U_f), tau_(0.1) {};

    Real getAxisVelocity(const Vecd &input_position, const Real &input_axis_velocity, Real time)
    {
        Real y_centered = input_position[1];
        Real u_steady = U_f_ * (1.0 - math::pow((Real(2) * y_centered / DH_), Real(2)));
        Real transient_factor = 1.0 - math::exp(-time / tau_);
        return u_steady * transient_factor;
    };

    Real DH_, U_f_, tau_;
};
//----------------------------------------------------------------------
//  Helper function for the analytical solution.
//----------------------------------------------------------------------
Real poiseuille_2d_u_steady(Real y)
{
    // Shift y so that y_centered = 0 at the channel center.
    Real y_centered = y - 0.5 * DH;
    return U_f * (1.0 - std::pow((2.0 * y_centered / DH), 2));
}
//----------------------------------------------------------------------
//  Outlet pressure condition classes for right side.
//----------------------------------------------------------------------
class InletInflowPressureConditionRight : public BaseStateCondition
{
  public:
    InletInflowPressureConditionRight(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser) {}

        Real operator()(size_t /*index_i*/, Real /*time*/)
        {
            return Outlet_pressure;
        }
    };
};
//----------------------------------------------------------------------
//  Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name)
        : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_body_shape;
        water_body_shape.emplace_back(0.0, 0.0);
        water_body_shape.emplace_back(0.0, DH);
        water_body_shape.emplace_back(DL, DH);
        water_body_shape.emplace_back(DL, 0.0);
        water_body_shape.emplace_back(0.0, 0.0);

        multi_polygon_.addPolygon(water_body_shape, GeometricOps::add);
    }
};
//----------------------------------------------------------------------
//  Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name)
        : MultiPolygonShape(shape_name)
    {
        // Outer boundary
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.emplace_back(0.0, -BW);
        outer_wall_shape.emplace_back(0.0, DH + BW);
        outer_wall_shape.emplace_back(DL, DH + BW);
        outer_wall_shape.emplace_back(DL, -BW);
        outer_wall_shape.emplace_back(0.0, -BW);

        // Inner boundary
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.emplace_back(-BW, 0.0);
        inner_wall_shape.emplace_back(-BW, DH);
        inner_wall_shape.emplace_back(DL + BW, DH);
        inner_wall_shape.emplace_back(DL + BW, 0.0);
        inner_wall_shape.emplace_back(-BW, 0.0);

        multi_polygon_.addPolygon(outer_wall_shape, GeometricOps::add);
        multi_polygon_.addPolygon(inner_wall_shape, GeometricOps::sub);
    }
};
//----------------------------------------------------------------------
//  Validate velocity from observer with analytical solution
//----------------------------------------------------------------------
int velocity_validation(
    const std::vector<Vecd> &observer_location,
    const std::vector<Vecd> &observer_vel,
    Real (*analytical_solution)(Real),
    Real tolerance_factor,
    Real U_f)
{
    size_t total_passed = 0;
    size_t total_failed = 0;
    std::vector<std::string> messages;

    // Loop over each observer point and compare the x-component of the velocity.
    for (size_t index = 0; index < observer_location.size(); ++index)
    {
        Real y = observer_location[index][1];
        Real vel_x_analytical = analytical_solution(y);
        Real vel_x_simulation = observer_vel[index][0];

        Real error = std::abs((vel_x_simulation - vel_x_analytical) / U_f);
        std::ostringstream msg;
        msg << "Measure at observer index " << index
            << " | Analytical: " << vel_x_analytical
            << " | Simulation: " << vel_x_simulation
            << " | Error: " << error;
        messages.push_back(msg.str());

        if (error <= tolerance_factor)
        {
            total_passed++;
        }
        else
        {
            total_failed++;
        }
    }

    // Print summary
    std::cout << "Detailed error measures:\n";
    for (const auto &msg : messages)
    {
        std::cout << msg << "\n";
    }
    std::cout << "[TEST SUMMARY] Velocity Validation:\n"
              << "Total Observations: " << observer_location.size() << "\n"
              << "Passed: " << total_passed << "\n"
              << "Failed: " << total_failed << "\n";

    // Final assertion for unit testing
    if (total_failed != 0)
    {
        std::cout << "Test failed with " << total_failed << " mismatches. Check log for details.";
        return 1;
    }
    return 0;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    applyResolutionFromEnv();
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_body.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> particle_buffer(0.5);
    water_body.generateParticlesWithReserve<BaseParticles, Lattice>(particle_buffer);

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMatterMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();
    // Add observer
    {
        int num_points = 15;
        // Avoid deploy observer too close to wall
        Real y_start = 2.0 * global_resolution;
        Real y_end = DH - 2.0 * global_resolution;
        Real total_range = y_end - y_start;
        Real dy = total_range / (num_points - 1);

        for (int i = 0; i < num_points; ++i)
        {
            Real y_i = y_start + i * dy;
            observer_location.push_back(Vecd(0.5 * DL, y_i));
        }
    }

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ObserverParticles>(observer_location);
    // //----------------------------------------------------------------------
    // //	Creating body parts.
    // //----------------------------------------------------------------------
    OrientedBoxByCell left_emitter_by_cell(water_body, OrientedBox(xAxis, Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    left_emitter_by_cell.writeOrientedBoxToVtp();
    OrientedBoxByCell right_emitter_by_cell(water_body, OrientedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_disposer_translation)), bidirectional_buffer_halfsize));
    right_emitter_by_cell.writeOrientedBoxToVtp();

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    // ----------------------------------------------------------------------
    Inner<> water_body_inner(water_body);
    Contact<> water_wall_contact(water_body, wall);
    Contact<> velocity_observer_contact(velocity_observer, water_body);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_body_update_complex_relation(water_body_inner, water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(velocity_observer_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_body);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal_direction(wall); // run on CPU
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_body);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_update_particle_position(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_body_inner, 0.5), water_wall_contact);   
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionForOpenBoundaryFlowCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::CompressionSummation<Inner<>, Contact<>>>
        fluid_density_summation(water_body_inner, water_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, WeaklyCompressibleFluid, Internal, ExcludeBufferParticles>>
        fluid_density_regularization(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, KernelGradientIntegralCorrectedForOpenBoundaryFlowComplex> kernel_gradient_integral(water_body_inner, water_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionCK<SPHBody, TruncatedLinear, BulkParticles>> transport_correction(water_body);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionViscousTimeStepCK> fluid_advection_time_step(water_body, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>> fluid_acoustic_time_step(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_body_inner, water_wall_contact);
    fluid_dynamics::BidirectionalBoundaryCK<MainExecutionPolicy, LinearCorrectionCK, InflowVelocityPrescribed>
        bidirectional_velocity_condition_left(left_emitter_by_cell, DH, U_f, mu_f);
    fluid_dynamics::BidirectionalBoundaryCK<MainExecutionPolicy, LinearCorrectionCK, PressurePrescribed<WeaklyCompressibleFluid>>
        bidirectional_pressure_condition_right(right_emitter_by_cell, Outlet_pressure);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::OutflowParticleDeletion> out_flow_particle_deletion(water_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_body, "Pressure");
    body_states_recording.addToWrite<int>(water_body, "BufferIndicator");
    ObservedQuantityRecording<MainExecutionPolicy, Vecd, RestoringCorrection> write_centerline_velocity(velocity_observer_contact, "Velocity");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_normal_direction.exec();
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_body_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();
    fluid_boundary_indicator.exec();
    bidirectional_velocity_condition_left.tagBufferParticles();
    bidirectional_pressure_condition_right.tagBufferParticles();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control + GPU-paper recording
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    size_t screen_output_interval = 100;
    // Duration V2: formal A1 end_time = 4 s (source default was 2 s).
    const char *end_env = std::getenv("GPU_PAPER_END_TIME");
    Real end_time = end_env ? Real(std::atof(end_env)) : Real(4.0);
    Real output_interval = end_time; // default: one coarse block; VTP gated below
    const char *out_iv = std::getenv("GPU_PAPER_OUTPUT_INTERVAL");
    if (out_iv)
        output_interval = Real(std::atof(out_iv));

    std::string run_id = std::getenv("GPU_PAPER_RUN_ID") ? std::getenv("GPU_PAPER_RUN_ID") : "A1_smoke";
    std::string mode = std::getenv("GPU_PAPER_MODE") ? std::getenv("GPU_PAPER_MODE") : "validate";
    std::string run_dir = std::getenv("GPU_PAPER_RUN_DIR")
                              ? std::getenv("GPU_PAPER_RUN_DIR")
                              : (std::string("gpu_paper_runs/") + run_id);
    // Default: original-case physical time T=2 + V4 uniform100 bind (override via env).
    std::string schedule_path =
        std::getenv("GPU_PAPER_SCHEDULE_JSON")
            ? std::getenv("GPU_PAPER_SCHEDULE_JSON")
            : "/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/results/A1_poiseuille/schedules/"
              "A1_validate_T4_bound.json";
    bool enable_vtp = (mode != "benchmark");
    if (std::getenv("GPU_PAPER_ENABLE_VTP"))
        enable_vtp = std::string(std::getenv("GPU_PAPER_ENABLE_VTP")) == "1";

    RunRecorder recorder(run_dir, run_id, "A1", mode);
    recorder.setMeta("case_name", "2d_gpu_paper_a1_poiseuille");
    recorder.setMeta("sizeof_Real", std::to_string(sizeof(Real)));
    recorder.setMeta("global_resolution", std::to_string(global_resolution));
    recorder.setMeta("DH_over_dp", std::to_string(dh_over_dp));
    recorder.setMeta("backend",
#if SPHINXSYS_USE_SYCL
                     "sycl_gpu"
#else
                     "cpu_ck"
#endif
    );
    recorder.setMeta("Re", std::to_string(Re));
    recorder.setMeta("end_time", std::to_string(end_time));
    recorder.setMeta("Inlet_pressure", std::to_string(Inlet_pressure));
    recorder.setMeta("Outlet_pressure", std::to_string(Outlet_pressure));
    recorder.setMeta("U_f", std::to_string(U_f));
    recorder.setMeta("mu_f", std::to_string(mu_f));
    recorder.writeHeaders();
    recorder.addTransfer("solver", "unresolved", "na", 0, 0, 0, 0, "not_collected",
                         "no_explicit_transfer_API_in_this_path", "seeded_at_start");
    recorder.addSync("solver", "wait_paths_unresolved", 0, 0, true, "not_collected",
                     "detail deferred to profile short window policy");
    recorder.addKernelProfilePlaceholder("na", "suite_default", "not_collected",
                                         "only_on_B2_B3_profile_window");
    recorder.addEnergyPlaceholder("gpu_window", "not_collected", "needs_valid_telemetry_coverage");
    recorder.registerMissing({"inlet_outlet_cumulative_mass_flux", "not_collected",
                              "Bidirectional buffer has inject/delete; dedicated flux integrator not wired yet",
                              "case_a1", "add buffer mass accounting if needed for paper C checks"});

    std::vector<double> probe_times;
    std::vector<double> field_times;
    probe_times = loadProbeTimesFromBoundSchedule(schedule_path);
    field_times = loadFieldTimesFromBoundSchedule(schedule_path);
    recorder.setMeta("schedule_json", schedule_path);
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
        write_centerline_velocity.writeToFile(number_of_iterations);
        fluid_observer_contact_relation.exec();
        auto observer_vel = velocity_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");

        Real u_center = 0.0;
        if (!observer_location.empty())
            u_center = observer_vel[observer_location.size() / 2][0];

        // Trapezoidal Q ~ integral u dy across observer stations (2D channel unit depth).
        Real Q_est = 0.0;
        for (size_t i = 0; i + 1 < observer_location.size(); ++i)
        {
            Real y0 = observer_location[i][1];
            Real y1 = observer_location[i + 1][1];
            Real u0 = observer_vel[i][0];
            Real u1 = observer_vel[i + 1][0];
            Q_est += Real(0.5) * (u0 + u1) * (y1 - y0);
        }
        Real Q_analytical = (Real(2.0) / Real(3.0)) * U_f * DH;

        auto emit = [&](const std::string &probe_id, const std::string &quantity, Real value,
                        const std::string &unit, Real x, Real y)
        {
            ProbeSample ps;
            ps.physical_time_s = t;
            ps.requested_time_s = requested;
            ps.probe_id = probe_id;
            ps.body = "WaterBody";
            ps.quantity = quantity;
            ps.value = double(value);
            ps.unit = unit;
            ps.x = double(x);
            ps.y = double(y);
            recorder.addProbe(ps);
        };

        emit("centerline_mid", "u", u_center, "m/s", 0.5 * DL, 0.5 * DH);
        emit("section_mid", "Q_observer_trap", Q_est, "m^2/s", 0.5 * DL, 0.5 * DH);
        emit("section_mid", "Q_analytical_steady", Q_analytical, "m^2/s", 0.5 * DL, 0.5 * DH);
        emit("bc", "pressure_drop_prescribed", Inlet_pressure - Outlet_pressure, "Pa", 0.0, 0.5 * DH);

        for (size_t i = 0; i < observer_location.size(); ++i)
        {
            emit("profile_y" + std::to_string(i), "u", observer_vel[i][0], "m/s",
                 observer_location[i][0], observer_location[i][1]);
            emit("profile_y" + std::to_string(i), "u_analytical_steady",
                 poiseuille_2d_u_steady(observer_location[i][1]), "m/s",
                 observer_location[i][0], observer_location[i][1]);
        }

        recorder.addDiagnostic(t, "WaterBody", "N_active",
                               double(water_body.getBaseParticles().TotalRealParticles()), "1", "measured");
        recorder.addDiagnostic(t, "WaterBody", "u_center", double(u_center), "m/s", "measured");
        recorder.addDiagnostic(t, "WaterBody", "Q_observer_trap", double(Q_est), "m^2/s", "measured");
    };

    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    auto wall0 = std::chrono::steady_clock::now();
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    if (enable_vtp)
        body_states_recording.writeToFile();
    double req0 = probe_times.empty() ? 0.0 : probe_times[0];
    sample_probes(sv_physical_time->getValue(), req0);
    if (!probe_times.empty())
        probe_cursor = 1;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            fluid_density_summation.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            kernel_gradient_integral.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
            time_instance = TickCount::now();
            Real advection_dt = fluid_advection_time_step.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            time_instance = TickCount::now();
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = SMIN(fluid_acoustic_time_step.exec(), advection_dt);
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                bidirectional_velocity_condition_left.applyBoundaryCondition(acoustic_dt);
                bidirectional_pressure_condition_right.applyBoundaryCondition(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            water_update_particle_position.exec();
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	Dt = " << advection_dt << "	dt = " << acoustic_dt << "\n";
            }
            number_of_iterations++;
            /** inflow emitter injection*/
            bidirectional_velocity_condition_left.injectParticles();
            bidirectional_pressure_condition_right.injectParticles();
            bidirectional_velocity_condition_left.indicateOutFlowParticles();
            bidirectional_pressure_condition_right.indicateOutFlowParticles();
            out_flow_particle_deletion.exec();
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
            fluid_boundary_indicator.exec();
            bidirectional_velocity_condition_left.tagBufferParticles();
            bidirectional_pressure_condition_right.tagBufferParticles();

            double t_now = sv_physical_time->getValue();
            double requested = 0.0;
            if (dueForSample(t_now, probe_times, probe_cursor, &requested))
                sample_probes(t_now, requested);
            if (enable_vtp && dueForSample(t_now, field_times, field_cursor, &requested))
                body_states_recording.writeToFile();
        }

        TickCount t2 = TickCount::now();
        // VTP only on bound field schedule (Duration V2), not every output_interval.
        fluid_observer_contact_relation.exec();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    recorder.stages().add("computing_time_step", interval_computing_time_step.seconds());
    recorder.stages().add("pressure_relaxation", interval_computing_pressure_relaxation.seconds());
    recorder.stages().add("updating_configuration", interval_updating_configuration.seconds());
    recorder.stages().add("vtp_or_io", interval.seconds());

    recorder.addStepBlock("advection_outer", int(number_of_iterations), 0, 0, 0, 0, 0, 0, "mixed", "measured");
    recorder.addBodyStats(sv_physical_time->getValue(), "WaterBody",
                          int(water_body.getBaseParticles().TotalRealParticles()),
                          int(water_body.getBaseParticles().ParticlesBound()),
                          "not_sampled", 0, 0, 0, 0, 0, 0, "end_only");

    //----------------------------------------------------------------------
    //	Validation against analytical steady profile (original-case criterion)
    //----------------------------------------------------------------------
    fluid_observer_contact_relation.exec();
    auto observer_vel = velocity_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
    std::vector<Vecd> observer_vel_vec(observer_vel, observer_vel + observer_location.size());
    Real error_tolerance = 5 * 0.01; // Less than 5 percent when resolution is DH/20
    int validate_rc =
        velocity_validation(observer_location, observer_vel_vec, poiseuille_2d_u_steady, error_tolerance, U_f);

    // Final profile L2/Linf vs analytical (normalized by U_f)
    Real l2_num = 0.0, l2_den = 0.0, linf = 0.0;
    for (size_t i = 0; i < observer_location.size(); ++i)
    {
        Real e = (observer_vel_vec[i][0] - poiseuille_2d_u_steady(observer_location[i][1])) / U_f;
        l2_num += e * e;
        l2_den += 1.0;
        linf = std::max(linf, std::abs(e));
    }
    Real l2 = (l2_den > 0.0) ? std::sqrt(l2_num / l2_den) : 0.0;
    {
        std::ofstream ofs(run_dir + "/accuracy_cost.csv", std::ios::app);
        ofs << run_id << ",analytical_poiseuille,u_profile,L2," << l2 << ",U_f," << tt.seconds() << ","
            << global_resolution << "," << water_body.getBaseParticles().TotalRealParticles() << ","
            << (validate_rc == 0 ? "measured" : "measured_failed_tolerance") << "\n";
        ofs << run_id << ",analytical_poiseuille,u_profile,Linf," << linf << ",U_f," << tt.seconds() << ","
            << global_resolution << "," << water_body.getBaseParticles().TotalRealParticles() << ","
            << (validate_rc == 0 ? "measured" : "measured_failed_tolerance") << "\n";
    }
    recorder.addDiagnostic(sv_physical_time->getValue(), "WaterBody", "profile_error_L2_over_Uf", double(l2),
                           "1", "measured");
    recorder.addDiagnostic(sv_physical_time->getValue(), "WaterBody", "profile_error_Linf_over_Uf", double(linf),
                           "1", "measured");

    auto wall1 = std::chrono::steady_clock::now();
    double wall_process = std::chrono::duration<double>(wall1 - wall0).count();
    const std::string final_status = (validate_rc == 0) ? "ok" : "physics_validation_failed";
    recorder.finalize(final_status, 0.0, sv_physical_time->getValue(), wall_process,
                      0.0, 0.0, 0.0, tt.seconds(), tt.seconds() + interval.seconds(),
                      0.0, interval.seconds(), 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0);
    return validate_rc;
}
