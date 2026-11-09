/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed Poiseuille flow example.
 * @details This is a basic test case for mixed pressure/velocity inlet and outlet boundary conditions.
 * @author 	YuVirtonomy, Xiangyu Hu
 */

#include "sphinxsys.h" // SPHinXsys Library.
#include "benchmark_config.h"
#include "benchmark_recorder.h"
#include <algorithm>
using namespace SPH;

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.004;                    /**< Channel length. */
Real DH = 0.001;                    /**< Channel height. */
Real global_resolution = DH / 20.0; /**< Reference particle spacing. */
Real BW = global_resolution * 4;    /**< Extending width for BCs. */
StdVec<Vecd> observer_location;
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
const Real Inlet_pressure = 0.5;
const Real Outlet_pressure = -0.5;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = std::sqrt(rho0_f * std::pow(0.5 * DH, 3.0) *
                      std::abs(Inlet_pressure - Outlet_pressure) / (Re * DL));

/**
 * Analytical solution for a laminar Poiseuille flow in a 2D channel:
 *
 *    U_f = Δp * DH^2 / (8 * μ * L)
 *
 *  where:
 *    Δp = |p_in - p_out|
 *    DH  = channel height
 *    μ   = dynamic viscosity
 *    L   = channel length
 */
Real U_f = (DH * DH * std::abs(Inlet_pressure - Outlet_pressure)) /
           (8.0 * mu_f * DL);

// Compute speed of sound (c0) based on the pressure difference between inlet and outlet boundaries.
// Ensures density variations are limited to ~1% (WCSPH criterion), multiplied by 4 as a safety factor.
Real c_f = std::max(10.0 * U_f, sqrt(4 * (Inlet_pressure - Outlet_pressure) / (rho0_f * 0.01))); //

//----------------------------------------------------------------------
//  Geometric shapes for the channel and boundaries.
//----------------------------------------------------------------------
Real bidirectional_buffer_length = 3.0 * global_resolution;
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
struct VelocityValidationResult
{
    size_t passed = 0;
    size_t failed = 0;
    Real maximum_normalized_error = 0.0;
    double console_seconds = 0.0;

    bool passedAll() const { return failed == 0; }
};

VelocityValidationResult velocity_validation(
    const std::vector<Vecd> &observer_location,
    const std::vector<Vecd> &observer_vel,
    Real (*analytical_solution)(Real),
    Real tolerance_factor,
    Real U_f)
{
    VelocityValidationResult result;
    std::vector<std::string> messages;

    // Loop over each observer point and compare the x-component of the velocity.
    for (size_t index = 0; index < observer_location.size(); ++index)
    {
        Real y = observer_location[index][1];
        Real vel_x_analytical = analytical_solution(y);
        Real vel_x_simulation = observer_vel[index][0];

        Real error = std::abs((vel_x_simulation - vel_x_analytical) / U_f);
        result.maximum_normalized_error =
            std::max(result.maximum_normalized_error, error);
        std::ostringstream msg;
        msg << "Measure at observer index " << index
            << " | Analytical: " << vel_x_analytical
            << " | Simulation: " << vel_x_simulation
            << " | Error: " << error;
        messages.push_back(msg.str());

        if (error <= tolerance_factor)
        {
            result.passed++;
        }
        else
        {
            result.failed++;
        }
    }

    // Print summary
    TickCount console_io_start = TickCount::now();
    std::cout << "Detailed error measures:\n";
    for (const auto &msg : messages)
    {
        std::cout << msg << "\n";
    }
    std::cout << "[TEST SUMMARY] Velocity Validation:\n"
              << "Total Observations: " << observer_location.size() << "\n"
              << "Passed: " << result.passed << "\n"
              << "Failed: " << result.failed << "\n";

    // Final assertion for unit testing
    if (!result.passedAll())
    {
        std::cout << "Test failed with " << result.failed << " mismatches. Check log for details.";
    }
    result.console_seconds =
        (TickCount::now() - console_io_start).seconds();
    return result;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    const paper_bench::BenchmarkDefaults benchmark_defaults{
        DH / 20.0, 2.0, 0.25,
        {{"standard", DH / 20.0},
         {"coarse", DH / 10.0},
         {"medium", DH / 20.0},
         {"fine", DH / 40.0}}};
    paper_bench::BenchmarkConfig benchmark_config =
        paper_bench::BenchmarkConfig::parse(ac, av, benchmark_defaults);
    TickCount benchmark_start = TickCount::now();
    global_resolution = benchmark_config.dp;
    BW = global_resolution * 4.0;
    bidirectional_buffer_length = 3.0 * global_resolution;
    bidirectional_buffer_halfsize =
        Vec2d(0.5 * bidirectional_buffer_length, 0.5 * DH);
    left_bidirectional_translation = bidirectional_buffer_halfsize;
    right_bidirectional_translation =
        Vec2d(DL - 0.5 * bidirectional_buffer_length, 0.5 * DH);
    right_disposer_translation = right_bidirectional_translation;
    observer_location.clear();
    TickCount environment_io_start = TickCount::now();
    paper_bench::BenchmarkRecorder benchmark_recorder(
        "poiseuille_2d", benchmark_config);
    benchmark_recorder.activateRunDirectory();
    TimeInterval io_interval = TickCount::now() - environment_io_start;
    const double environment_io_seconds = io_interval.seconds();
    TickCount initialization_start = TickCount::now();
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(
        Vec2d(-2.0 * BW, -2.0 * BW),
        Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    benchmark_recorder.stageInputAssets();
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
    const size_t initial_particle_count =
        water_body.getBaseParticles().TotalRealParticles() +
        wall.getBaseParticles().TotalRealParticles() +
        velocity_observer.getBaseParticles().TotalRealParticles();
    // //----------------------------------------------------------------------
    // //	Creating body parts.
    // //----------------------------------------------------------------------
    OrientedBoxByCell left_emitter_by_cell(water_body, OrientedBox(xAxis, Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    OrientedBoxByCell right_emitter_by_cell(water_body, OrientedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_disposer_translation)), bidirectional_buffer_halfsize));
    if (benchmark_config.output_enabled)
    {
        TickCount geometry_io_start = TickCount::now();
        left_emitter_by_cell.writeOrientedBoxToVtp();
        right_emitter_by_cell.writeOrientedBoxToVtp();
        io_interval += TickCount::now() - geometry_io_start;
    }

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    // ----------------------------------------------------------------------
    Inner<> water_body_inner(water_body);
    Contact<> water_wall_contact(water_body, {&wall});
    Contact<> velocity_observer_contact(velocity_observer, {&water_body});
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
    StateDynamics<MainExecutionPolicy, LinearCorrectionMatrixScope<SPHBody, BulkParticles>>
        fluid_linear_correction_scope(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::CompressionSummation<Inner<>, Contact<>>>
        fluid_density_summation(water_body_inner, water_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::DensityRegularization<SPHBody, WeaklyCompressibleFluid, Internal, ExcludeBufferParticles>>
        fluid_density_regularization(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, KernelGradientIntegralCorrectedComplex> kernel_gradient_integral(water_body_inner, water_wall_contact);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionCK<SPHBody, TruncatedLinear, BulkParticles>> transport_correction(water_body);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_f);
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
    const double init_seconds = std::max(
        0.0, (TickCount::now() - initialization_start).seconds() -
                 (io_interval.seconds() - environment_io_seconds));
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    size_t number_of_outer_steps = 0;
    size_t number_of_acoustic_steps = 0;
    size_t screen_output_interval = 100;
    size_t observation_sample_interval = screen_output_interval * 2;
    Real end_time = Real(benchmark_config.end_time);
    Real output_interval = Real(benchmark_config.output_interval);
    const bool write_periodic_output = benchmark_config.output_enabled;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    if (write_periodic_output)
    {
        TickCount initial_io_start = TickCount::now();
        body_states_recording.writeToFile();
        write_centerline_velocity.writeToFile(number_of_iterations);
        io_interval += TickCount::now() - initial_io_start;
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval &&
               sv_physical_time->getValue() < end_time)
        {
            fluid_density_summation.exec();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            fluid_linear_correction_scope.exec();
            kernel_gradient_integral.exec();
            transport_correction.exec();
            fluid_viscous_force.exec();
            Real advection_dt = fluid_advection_time_step.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt &&
                   integration_time < output_interval &&
                   sv_physical_time->getValue() < end_time)
            {
                acoustic_dt = std::min(
                    {fluid_acoustic_time_step.exec(),
                     advection_dt - relaxation_time,
                     output_interval - integration_time,
                     end_time - sv_physical_time->getValue()});
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                bidirectional_velocity_condition_left.applyBoundaryCondition(acoustic_dt);
                bidirectional_pressure_condition_right.applyBoundaryCondition(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
                number_of_acoustic_steps++;
            }
            water_update_particle_position.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                TickCount console_io_start = TickCount::now();
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	Dt = " << advection_dt << "	dt = " << acoustic_dt << "\n";
                io_interval += TickCount::now() - console_io_start;
                if (write_periodic_output &&
                    number_of_iterations % observation_sample_interval == 0 &&
                    number_of_iterations != sph_system.RestartStep())
                {
                    TickCount observer_io_start = TickCount::now();
                    write_centerline_velocity.writeToFile(number_of_iterations);
                    io_interval += TickCount::now() - observer_io_start;
                }
            }
            number_of_iterations++;
            /** inflow emitter injection*/
            bidirectional_velocity_condition_left.injectParticles();
            bidirectional_pressure_condition_right.injectParticles();
            bidirectional_velocity_condition_left.indicateOutFlowParticles();
            bidirectional_pressure_condition_right.indicateOutFlowParticles();
            out_flow_particle_deletion.exec();
            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            fluid_boundary_indicator.exec();
            bidirectional_velocity_condition_left.tagBufferParticles();
            bidirectional_pressure_condition_right.tagBufferParticles();
        }

        if (write_periodic_output)
        {
            TickCount periodic_io_start = TickCount::now();
            body_states_recording.writeToFile();
            io_interval += TickCount::now() - periodic_io_start;
        }
        number_of_outer_steps++;
    }
    // A final observer sample is required by the analytical validation. This
    // is the only observer file write performed by benchmark mode.
    TickCount validation_sample_start = TickCount::now();
    write_centerline_velocity.writeToFile(number_of_iterations);
    io_interval += TickCount::now() - validation_sample_start;
    //----------------------------------------------------------------------
    //	GTest-based validation against analytical solution
    //----------------------------------------------------------------------
    TickCount validation_start = TickCount::now();
    // Get the velocity data from the observer body particles
    auto observer_vel = velocity_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
    // Validate observer velocities against analytical Poiseuille profile
    // Convert the pointer to a std::vector using the number of observer particles.
    std::vector<Vecd> observer_vel_vec(observer_vel, observer_vel + observer_location.size());
    Real error_tolerance = 5 * 0.01; // Less than 5 percent when resolution is DH/20
    const VelocityValidationResult validation =
        velocity_validation(observer_location, observer_vel_vec,
                            poiseuille_2d_u_steady, error_tolerance, U_f);
    const double validation_total_seconds =
        (TickCount::now() - validation_start).seconds();
    const double verification_seconds =
        std::max(0.0, validation_total_seconds - validation.console_seconds);
    const double io_seconds =
        io_interval.seconds() + validation.console_seconds;

    const TimeInterval total_wall = TickCount::now() - benchmark_start;
    const double compute_seconds = std::max(
        0.0, total_wall.seconds() - init_seconds -
                 io_seconds - verification_seconds);
    std::cout << "Total wall time: " << total_wall.seconds()
              << " seconds (compute: " << compute_seconds
              << ", initialization: " << init_seconds
              << ", I/O: " << io_seconds << ")." << std::endl;

    paper_bench::BenchmarkSummary summary;
    summary.dp = global_resolution;
    summary.initial_particle_count = initial_particle_count;
    summary.particle_count =
        water_body.getBaseParticles().TotalRealParticles() +
        wall.getBaseParticles().TotalRealParticles() +
        velocity_observer.getBaseParticles().TotalRealParticles();
    summary.physical_time = sv_physical_time->getValue();
    summary.outer_steps = number_of_outer_steps;
    summary.advection_steps = number_of_iterations;
    summary.acoustic_steps = number_of_acoustic_steps;
    summary.wall_seconds = total_wall.seconds();
    summary.compute_seconds = compute_seconds;
    summary.io_seconds = io_seconds;
    summary.init_seconds = init_seconds;
    summary.verification_seconds = verification_seconds;
    summary.time_per_outer_step =
        number_of_outer_steps == 0
            ? 0.0
            : compute_seconds / static_cast<double>(number_of_outer_steps);
    summary.status = validation.passedAll() ? "completed" : "validation_failed";
    summary.setParticleUpdates(
        summary.particle_count * summary.advection_steps,
        "particle_count*advection_steps");
    summary.extra_fields.push_back(
        {"validation_passed", validation.passedAll() ? "true" : "false"});
    summary.extra_fields.push_back(
        {"validation_pass_count", std::to_string(validation.passed)});
    summary.extra_fields.push_back(
        {"validation_fail_count", std::to_string(validation.failed)});
    summary.extra_fields.push_back(
        {"validation_max_normalized_error",
         std::to_string(validation.maximum_normalized_error)});
    summary.extra_fields.push_back(
        {"validation_tolerance", std::to_string(error_tolerance)});
    benchmark_recorder.writeSummary(summary);
    return validation.passedAll() ? 0 : 1;
}
