/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys_ck.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.004;                                             /**< Channel length. */
Real DH = 0.001;                                             /**< Channel height. */
Real resolution_ref = DH / 20.0;                             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
BoundingBox system_domain_bounds(Vec2d(-BW * 2, -BW * 2), Vec2d(DL + BW * 2, DH + BW * 2));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = -0.2;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
// Real U_f = 0.05;
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Real bidrectional_buffer_length = 3.0 * resolution_ref;
Vec2d bidirectional_buffer_halfsize = Vec2d(bidrectional_buffer_length * 0.5, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d left_indicator_translation = Vec2d(0.0, 0.5 * DH);

Vec2d right_bidirectional_translation = Vec2d(DL - 0.5 * bidrectional_buffer_length, 0.5 * DH);
Vec2d right_indicator_translation = Vec2d(DL, 0.5 * DH);

Vec2d right_disposer_translation = Vec2d(DL - 0.5 * bidrectional_buffer_length, 0.5 * DH);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public BaseStateCondition
{
  public:
    InletInflowCondition(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser){};

        void operator()(AlignedBox *aligned_box, UnsignedInt index_i)
        {
            vel_[index_i] = Vec2d(U_f, 0.0);
        };
    };
};
//----------------------------------------------------------------------
//	InletInflowpPressureCondition
//----------------------------------------------------------------------
class InletInflowpPressureCondition : public BaseStateCondition
{
  public:
    InletInflowpPressureCondition(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser){};

        Real operator()(size_t index_i, Real time)
        {
            Real rise_time = 0.0; // Time duration for full sine increase
            Real p = 0.0;
            Real target_pressure_ = Inlet_pressure;

            if (time < rise_time)
            {
                p = target_pressure_ * 0.5 * (1.0 - cos(M_PI * time / rise_time));
            }
            else
            {
                p = target_pressure_;
            }

            // p_[index_i] = p;
            return p;
        };
    };
};

class InletInflowpPressureConditionRight : public BaseStateCondition
{
  public:
    InletInflowpPressureConditionRight(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser){};

        Real operator()(size_t index_i, Real time)
        {
            Real rise_time = 0.0; // Time duration for full sine increase
            Real p = 0.0;
            Real target_pressure_ = Outlet_pressure;

            if (time < rise_time)
            {
                p = target_pressure_ * 0.5 * (1.0 - cos(M_PI * time / rise_time));
            }
            else
            {
                p = target_pressure_;
            }

            // p_[index_i] = p;
            return p;
        };
    };
};

//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_body_shape;
        water_body_shape.push_back(Vecd(0.0, 0.0));
        water_body_shape.push_back(Vecd(0.0, DH));
        water_body_shape.push_back(Vecd(DL, DH));
        water_body_shape.push_back(Vecd(DL, 0.0));
        water_body_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(0.0, -BW));
        outer_wall_shape.push_back(Vecd(0.0, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, -BW));
        outer_wall_shape.push_back(Vecd(0.0, -BW));
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(-BW, 0.0));
        inner_wall_shape.push_back(Vecd(-BW, DH));
        inner_wall_shape.push_back(Vecd(DL + BW, DH));
        inner_wall_shape.push_back(Vecd(DL + BW, 0.0));
        inner_wall_shape.push_back(Vecd(-BW, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

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

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ObserverParticles>(observer_location);
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    AlignedBoxPartByParticle left_emitter_by_particle(water_body, AlignedBox(xAxis, Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell left_emitter_by_cell(water_body, AlignedBox(xAxis, Transform(left_bidirectional_translation), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell right_emitter_by_cell(water_body, AlignedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_disposer_translation)), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell left_indicator_by_cell(water_body, AlignedBox(xAxis, Transform(left_indicator_translation), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell right_indicator_by_cell(water_body, AlignedBox(xAxis, Transform(Vec2d(right_indicator_translation)), bidirectional_buffer_halfsize));
    AlignedBoxPartByCell right_disposer(water_body, AlignedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_disposer_translation)), bidirectional_buffer_halfsize));

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
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
    StateDynamics<MainExecutionPolicy, fluid_dynamics::SurfaceIndicationByAlignedBoxCK> label_left_indicator(left_indicator_by_cell);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::SurfaceIndicationByAlignedBoxCK> label_right_indicator(right_indicator_by_cell);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK>
        transport_correction_ck(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCKWithoutUpdate>
        zero_gradient_ck(water_body_inner, water_wall_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_body, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK> fluid_acoustic_time_step(water_body);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_body_inner, water_wall_contact);

    StateDynamics<MainExecutionPolicy, fluid_dynamics::TagBufferParticlesCK> left_tag_buffer_particle_(left_emitter_by_cell);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::TagBufferParticlesCK> right_tag_buffer_particle_(right_emitter_by_cell);

    // template class BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, InletInflowpPressureCondition>;
    // template class BufferEmitterInflowInjectionCK<AlignedBoxPartByParticle, InletInflowpPressureCondition>;
    StateDynamics<MainExecutionPolicy, fluid_dynamics::InflowConditionCK<AlignedBoxPartByCell, InletInflowCondition>> inflow_condition(left_emitter_by_cell);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::PressureConditionCK<AlignedBoxPartByCell, NoKernelCorrectionCK, InletInflowpPressureCondition>> pressure_condition(left_emitter_by_cell);
    StateDynamics<execution::SequencedPolicy, fluid_dynamics::BufferEmitterInflowInjectionCK<AlignedBoxPartByCell, InletInflowpPressureCondition>> emitter_injection(left_emitter_by_cell, inlet_buffer);
    // StateDynamics<execution::SequencedPolicy, fluid_dynamics::BufferEmitterInflowInjectionCK<AlignedBoxPartByCell>> emitter_injection(left_emitter_by_cell, inlet_buffer);
    // StateDynamics<execution::SequencedPolicy, fluid_dynamics::EmitterInflowInjectionCK<AlignedBoxPartByParticle>> emitter_injection(left_emitter_by_particle, inlet_buffer);
    StateDynamics<SequencedExecutionPolicy, fluid_dynamics::DisposerOutflowDeletionCK> right_remove_particles(right_disposer);

    fluid_dynamics::BidirectionalBufferCK<MainExecutionPolicy, NoKernelCorrectionCK, InletInflowpPressureCondition>
        bidirectional_buffer_left(left_emitter_by_cell, inlet_buffer);

    fluid_dynamics::BidirectionalBufferCK<MainExecutionPolicy, NoKernelCorrectionCK, InletInflowpPressureConditionRight>
        bidirectional_buffer_right(right_emitter_by_cell, inlet_buffer);

    // InnerRelation water_block_inner(water_body);
    // ContactRelation water_block_contact(water_body, {&wall});
    // InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_body, "Pressure");
    body_states_recording.addToWrite<int>(water_body, "Indicator");
    body_states_recording.addToWrite<Real>(water_body, "Density");
    body_states_recording.addToWrite<Vecd>(water_body, "ZeroGradientResidue");
    body_states_recording.addToWrite<int>(water_body, "BufferParticleIndicator");
    body_states_recording.addToWrite<int>(water_body, "PreviousSurfaceIndicator");
    body_states_recording.addToWrite<int>(water_body, "WithScopeVerify");
    body_states_recording.addToWrite<int>(water_body, "DensitySummationVerify");
    body_states_recording.addToWrite<Real>(water_body, "Mass");
    // body_states_recording.addToWrite<Vecd>(water_block, "KernelSummation");

    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd>> write_centerline_velocity("Velocity", velocity_observer_contact);
    auto vel_ = water_body.getBaseParticles().getVariableDataByName<Vecd>("Velocity");
    auto buffer_particle_indicator_ = water_body.getBaseParticles().getVariableDataByName<int>("BufferParticleIndicator");

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
    bidirectional_buffer_left.tagBufferParticles();
    bidirectional_buffer_right.tagBufferParticles();
    // right_tag_buffer_particle_.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 1000;
    int observation_sample_interval = 1000;
    Real end_time = 10.0;
    Real output_interval = end_time / 100;
    Real total_time = 0.0;
    Real relax_time = 1.0;
    /** statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
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

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = fluid_acoustic_time_step.exec();
                acoustic_dt = std::min(acoustic_dt, advection_dt - relaxation_time);
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                zero_gradient_ck.exec();
                bidirectional_buffer_left.applyPressureCondition(acoustic_dt);
                bidirectional_buffer_right.applyPressureCondition(acoustic_dt);
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);

                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            water_advection_step_close.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	Dt = " << advection_dt << "	dt = " << acoustic_dt << "\n";
            }
            number_of_iterations++;

            /** inflow emitter injection*/
            bidirectional_buffer_left.injectParticles();
            // bidirectional_buffer_right.injectParticles();
            // bidirectional_buffer_left.deleteParticles();
            std::cout << "before remove: " << number_of_iterations << " \n";

            bidirectional_buffer_right.deleteParticles();
            std::cout << "after remove: " << number_of_iterations << " \n";

            // emitter_injection.exec();
            // right_remove_particles.exec();
            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                std::cout << "particle_sort.exec(); \n";
                // particle_sort.exec();
            }

            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            fluid_boundary_indicator.exec();
            bidirectional_buffer_left.tagBufferParticles();
            bidirectional_buffer_right.tagBufferParticles();
            // right_tag_buffer_particle_.exec();
            // label_left_indicator.exec();
            // label_right_indicator.exec();
            body_states_recording.writeToFile(MainExecutionPolicy{});
        }

        TickCount t2 = TickCount::now();

        body_states_recording.writeToFile(MainExecutionPolicy{});
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    if (sph_system.GenerateRegressionData())
    {
        write_centerline_velocity.generateDataBase(1.0e-3);
    }
    else
    {
        write_centerline_velocity.testResult();
    }

    return 0;
}
