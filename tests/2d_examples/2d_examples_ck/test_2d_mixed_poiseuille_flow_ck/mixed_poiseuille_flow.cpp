/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed Poiseuille flow example.
 * @details This is a basic test case for mixed pressure/velocity inlet and outlet boundary conditions.
 * @author 	YuVirtonomy, Xiangyu Hu
 */

#include "sphinxsys_ck.h" // SPHinXsys Library.
#include "gtest/gtest.h"

using namespace SPH;

//----------------------------------------------------------------------
//  Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.004;                 /**< Channel length. */
Real DH = 0.001;                 /**< Channel height. */
Real resolution_ref = DH / 20.0; /**< Reference particle spacing. */
Real BW = resolution_ref * 4;    /**< Extending width for BCs. */
StdVec<Vecd> observer_location;
BoundingBox system_domain_bounds(
    Vec2d(-2.0 * BW, -2.0 * BW),
    Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));
//----------------------------------------------------------------------
//  Material parameters.
//----------------------------------------------------------------------
const Real Inlet_pressure = 0.2;
const Real Outlet_pressure = 0.1;
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

/** Choose a wave speed for the weakly compressible model. */
Real c_f = 10.0 * U_f;

//----------------------------------------------------------------------
//  Geometric shapes for the channel and boundaries.
//----------------------------------------------------------------------
Real bidirectional_buffer_length = 3.0 * resolution_ref;
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
class InletInflowConditionLeft : public BaseStateCondition
{
  public:
    InletInflowConditionLeft(BaseParticles *particles)
        : BaseStateCondition(particles),
          DH_(DH), U_f_(U_f), mu_f_(mu_f) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser),
              DH_ck_(encloser.DH_),
              U_f_ck_(encloser.U_f_),
              mu_f_ck_(encloser.mu_f_),
              tau_((DH_ck_ * DH_ck_) / (M_PI * M_PI * mu_f_ck_))
        {
        }

        void operator()(AlignedBox *aligned_box, UnsignedInt index_i, Real time)
        {
            // Shift the y-coordinate so that y_centered = 0 at the channel center.
            Real y_centered = pos_[index_i][1] - 0.5 * DH_ck_;

            // Steady-state analytical solution for Poiseuille flow
            Real u_steady = U_f_ck_ *
                            (1.0 - std::pow((2.0 * y_centered / DH_ck_), 2));

            // Transient factor that approaches 1 as time grows.
            Real transient_factor = 1.0 - std::exp(-time / tau_);

            // Time-dependent velocity profile
            Real u_ave = u_steady * transient_factor;
            if (aligned_box->checkInBounds(pos_[index_i]))
                vel_[index_i] = Vec2d(u_ave, 0.0);
        }

      protected:
        Real DH_ck_;
        Real U_f_ck_;
        Real mu_f_ck_;
        Real tau_;
    };
    Real DH_;
    Real U_f_;
    Real mu_f_;
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

        multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
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

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//  Validate velocity from observer with analytical solution
//----------------------------------------------------------------------
void velocity_validation(
    const std::vector<Vecd> &observer_location,
    const std::vector<Vecd> &observer_vel,
    Real (*analytical_solution)(Real),
    Real tolerance_factor,
    Real U_f)
{
    size_t total_passed = 0;
    size_t total_failed = 0;
    std::vector<std::string> failure_messages;

    // Loop over each observer point and compare the x-component of the velocity.
    for (size_t index = 0; index < observer_location.size(); ++index)
    {
        Real y = observer_location[index][1];
        Real vel_x_analytical = analytical_solution(y);
        Real vel_x_simulation = observer_vel[index][0];

        // Check if within tolerance
        if (std::abs((vel_x_simulation - vel_x_analytical) / vel_x_analytical) <= tolerance_factor)
        {
            total_passed++;
        }
        else
        {
            total_failed++;
            std::ostringstream msg;
            msg << "Mismatch at observer index " << index
                << " | Analytical: " << vel_x_analytical
                << " | Simulation: " << vel_x_simulation
                << " | Error: " << std::abs((vel_x_simulation - vel_x_analytical) / vel_x_analytical);
            failure_messages.push_back(msg.str());
        }
    }

    // Print summary
    std::cout << "[TEST SUMMARY] Velocity Validation:\n"
              << "Total Observations: " << observer_location.size() << "\n"
              << "Passed: " << total_passed << "\n"
              << "Failed: " << total_failed << "\n";

    // Print detailed failure messages if any
    if (!failure_messages.empty())
    {
        std::cout << "Detailed Failures:\n";
        for (const auto &msg : failure_messages)
        {
            std::cout << msg << "\n";
        }
    }

    // Final assertion for unit testing
    EXPECT_EQ(total_failed, 0) << "Test failed with " << total_failed << " mismatches. Check log for details.";
}

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_body.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_buffer(0.5);
    water_body.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_buffer);

    SolidBody wall(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Lattice>();
    // Add observer
    {
        int num_points = 15;
        // Avoid deploy observer too close to wall
        Real y_start = 2.0 * resolution_ref;
        Real y_end = DH - 2.0 * resolution_ref;
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
    AlignedBoxPartByCell left_emitter_by_cell(water_body, AlignedBox(xAxis, Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell right_emitter_by_cell(water_body, AlignedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_disposer_translation)), bidirectional_buffer_halfsize));

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    // ----------------------------------------------------------------------
    Relation<Inner<>> water_body_inner(water_body);
    Relation<Contact<>> water_wall_contact(water_body, {&wall});
    Relation<Contact<>> velocity_observer_contact(velocity_observer, {&water_body});
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    using SequencedExecutionPolicy = execution::SequencedPolicy;
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> water_cell_linked_list(water_body);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_body_update_complex_relation(water_body_inner, water_wall_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(velocity_observer_contact);
    ParticleSortCK<MainExecutionPolicy, QuickSort> particle_sort(water_body);
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
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_body_inner, 0.5), water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexInternalPressureBoundary>
        fluid_density_regularization(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK>
        transport_correction_ck(water_body_inner, water_wall_contact);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>> fluid_acoustic_time_step(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCKWithoutUpdate>
        zero_gradient_ck(water_body_inner, water_wall_contact);
    fluid_dynamics::VelocityBidirectionalConditionCK<MainExecutionPolicy, SequencedExecutionPolicy, NoKernelCorrectionCK, InletInflowConditionLeft>
        bidirectional_velocity_condition_left(left_emitter_by_cell, inlet_buffer);
    fluid_dynamics::PressureBidirectionalConditionCK<MainExecutionPolicy, SequencedExecutionPolicy, NoKernelCorrectionCK, InletInflowPressureConditionRight>
        bidirectional_pressure_condition_right(right_emitter_by_cell, inlet_buffer);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_body, "Pressure");
    ObservedQuantityRecording<MainExecutionPolicy, Vecd> write_centerline_velocity("Velocity", velocity_observer_contact);
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
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    size_t screen_output_interval = 100;
    size_t observation_sample_interval = screen_output_interval * 2;
    Real end_time = 2.0;
    Real output_interval = 0.1;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(MainExecutionPolicy{});
    write_centerline_velocity.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_viscous_force.exec();
            transport_correction_ck.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            fluid_linear_correction_matrix.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = SMIN(fluid_acoustic_time_step.exec(), advection_dt);
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                zero_gradient_ck.exec();
                bidirectional_velocity_condition_left.applyPressureCondition(acoustic_dt);
                bidirectional_velocity_condition_left.applyVelocityCondition();
                bidirectional_pressure_condition_right.applyPressureCondition(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            water_advection_step_close.exec();
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	Dt = " << advection_dt << "	dt = " << acoustic_dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;
            /** inflow emitter injection*/
            bidirectional_velocity_condition_left.injectParticles();
            bidirectional_pressure_condition_right.injectParticles();
            bidirectional_velocity_condition_left.deleteParticles();
            bidirectional_pressure_condition_right.deleteParticles();
            /** Update cell linked list and configuration. */
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
        }

        TickCount t2 = TickCount::now();

        body_states_recording.writeToFile(MainExecutionPolicy{});
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
    //----------------------------------------------------------------------
    //	GTest-based validation against analytical solution
    //----------------------------------------------------------------------
    // Get the velocity data from the observer body particles
    auto observer_vel = velocity_observer.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
    // Validate observer velocities against analytical Poiseuille profile
    // Convert the pointer to a std::vector using the number of observer particles.
    std::vector<Vecd> observer_vel_vec(observer_vel, observer_vel + observer_location.size());
    Real error_tolerance = 3 * 0.01; // Less than 3 percent when resolution is DH/20
    velocity_validation(observer_location, observer_vel_vec, poiseuille_2d_u_steady, error_tolerance, U_f);
    return 0;
}
