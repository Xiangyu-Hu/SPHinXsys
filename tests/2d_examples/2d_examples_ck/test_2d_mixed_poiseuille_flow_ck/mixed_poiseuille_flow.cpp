/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
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
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.1;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL - 2.5 * resolution_ref, 0.5 * DH);
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
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
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
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ObserverParticles>(observer_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Relation<Inner<>> water_block_inner(water_block);
    Relation<Contact<>> water_block_contact(water_block, {&wall_boundary});
    Relation<Contact<>> velocity_observer_contact(velocity_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the main execution policy for this case.
    //----------------------------------------------------------------------
    using MainExecutionPolicy = execution::ParallelPolicy;
    using SequencedExecutionPolicy = execution::SequencedPolicy;
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_body_update_complex_relation(water_block_inner, water_block_contact);
    UpdateRelation<MainExecutionPolicy, Contact<>> fluid_observer_contact_relation(velocity_observer_contact);
    ParticleSortCK<MainExecutionPolicy, QuickSort> particle_sort(water_block);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_block);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(InteractArgs(water_block_inner, 0.5), water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        fluid_density_regularization(water_block_inner, water_block_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_block_inner, water_block_contact);

    // InnerRelation water_block_inner_ori(water_block);
    // ContactRelation water_wall_contact_ori(water_block, {&wall_boundary});
    // ComplexRelation water_block_complex_ori(water_block_inner_ori, water_wall_contact_ori);

    // InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner_ori, water_wall_contact_ori);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK>
        transport_correction_ck(water_block_inner, water_block_contact);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step(water_block, U_f);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK> fluid_acoustic_time_step(water_block);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_block_inner, water_block_contact);

    AlignedBoxPartByParticle left_emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    StateDynamics<MainExecutionPolicy, fluid_dynamics::InflowConditionCK<AlignedBoxPartByParticle, InletInflowCondition>> inflow_condition(left_emitter);
    StateDynamics<SequencedExecutionPolicy, fluid_dynamics::EmitterInflowInjectionCK> emitter_injection(left_emitter, in_outlet_particle_buffer);

    AlignedBoxPartByCell right_disposer(water_block, AlignedBox(xAxis, Transform(Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    StateDynamics<SequencedExecutionPolicy, fluid_dynamics::DisposerOutflowDeletionCK> right_remove_particles(right_disposer, in_outlet_particle_buffer);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Vecd>(water_block, "ZeroGradientResidue");
    // body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<MainExecutionPolicy, Vecd>> write_centerline_velocity("Velocity", velocity_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_body_update_complex_relation.exec();
    fluid_observer_contact_relation.exec();
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
        while (integration_time < output_interval)
        {
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            fluid_viscous_force.exec();
            fluid_linear_correction_matrix.exec();
            transport_correction_ck.exec();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_acoustic_step_1st_half.exec(acoustic_dt);
                inflow_condition.exec();
                fluid_acoustic_step_2nd_half.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            water_advection_step_close.exec();
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << "	Total Time = " << total_time
                          << "	Physical Time = " << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;
            // std::cout << "298\n";
            /** inflow emitter injection*/
            emitter_injection.exec();
            // right_remove_particles.exec();

            // first do injection for all buffers
            // left_bidirection_buffer.injection.exec();
            // right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            // left_bidirection_buffer.deletion.exec();
            // right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                std::cout << "particle_sort\n";
                particle_sort.exec();
            }

            water_cell_linked_list.exec();
            water_body_update_complex_relation.exec();
            fluid_observer_contact_relation.exec();
            // water_block.updateCellLinkedList();
            // water_block_complex_ori.updateConfiguration();

            // std::cout << "319\n";
            body_states_recording.writeToFile(MainExecutionPolicy{});
            fluid_boundary_indicator.exec();

            // std::cout << "321\n";
            // left_bidirection_buffer.tag_buffer_particles.exec();
            // right_bidirection_buffer.tag_buffer_particles.exec();
        }
        // std::cout << "326\n";
        body_states_recording.writeToFile(MainExecutionPolicy{});
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
