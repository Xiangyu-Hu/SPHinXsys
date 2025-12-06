/*
 * @file modified_T_shaped_pipe.cpp
 * @brief This is the benchmark test of multi -inlet and multi - outlet.
 * @details We consider a flow with one inlet and two outlets in a T - shaped pipe in 2D.
 * @author Xiangyu Hu,Shuoguo Zhang
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.2;                                               /**< Reference length. */
Real DH = 0.1;                                               /**< Reference and the height of main channel. */
Real DL1 = 0.75 * DL;                                        /**< The length of the main channel. */
Real resolution_ref = 0.005;                                 /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;                        /**< Reference size of the emitter buffer to impose inflow condition. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */

//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real Outlet_pressure = 0;
Real rho0_f = 1000.0;                                                 /**< Reference density of fluid. */
Real Re = 100.0;                                                      /**< Reynolds number. */
Real U_f = 1.0;                                                       /**< Characteristic velocity. */
Real mu_f = rho0_f * U_f * DH / Re;                                   /**< Dynamics viscosity. */
Real c_f = 10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * (DL - DL1))); /** Reference sound speed needs to consider the flow speed in the narrow channels. */
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        return p;
    }
};

struct UpOutflowPressure
{
    template <class BoundaryConditionType>
    UpOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

struct DownOutflowPressure
{
    template <class BoundaryConditionType>
    DownOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real curent_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBox &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the water block in T shape polygon. */
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH),
    Vecd(DL, 2.0 * DH), Vecd(DL, -DH), Vecd(DL1, -DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0)};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
    Vecd(-DL_sponge, -BW), Vecd(-DL_sponge, DH + BW), Vecd(DL1 - BW, DH + BW), Vecd(DL1 - BW, 2.0 * DH),
    Vecd(DL + BW, 2.0 * DH), Vecd(DL + BW, -DH), Vecd(DL1 - BW, -DH), Vecd(DL1 - BW, -BW), Vecd(-DL_sponge, -BW)};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
    Vecd(-DL_sponge - BW, 0.0), Vecd(-DL_sponge - BW, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH + BW),
    Vecd(DL, 2.0 * DH + BW), Vecd(DL, -DH - BW), Vecd(DL1, -DH - BW), Vecd(DL1, 0.0), Vecd(-DL_sponge - BW, 0.0)};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
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
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -DH - BW), Vec2d(DL + BW, 2.0 * DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setGenerateRegressionData(false);
    sph_system.handleCommandlineOptions(ac, av);
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
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    //----------------------------------------------------------------------
    // Left buffer
    //----------------------------------------------------------------------
    Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize;
    AlignedBoxByCell left_emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    //----------------------------------------------------------------------
    // Up buffer
    //----------------------------------------------------------------------
    Vec2d up_buffer_halfsize = Vec2d(0.5 * BW, 0.75);
    Vec2d up_buffer_translation = Vec2d(0.5 * (DL + DL1), 2.0 * DH - 0.5 * BW);
    AlignedBoxByCell up_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(-0.5 * Pi), Vec2d(up_buffer_translation)), up_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<UpOutflowPressure> up_bidirection_buffer(up_emitter, in_outlet_particle_buffer);
    //----------------------------------------------------------------------
    // Down buffer
    //----------------------------------------------------------------------
    Vec2d down_buffer_halfsize = Vec2d(0.5 * BW, 0.75);
    Vec2d down_buffer_translation = Vec2d(0.5 * (DL + DL1), -DH + 0.5 * BW);
    AlignedBoxByCell down_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(0.5 * Pi), Vec2d(down_buffer_translation)), down_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<DownOutflowPressure> down_bidirection_buffer(down_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<UpOutflowPressure>> up_inflow_pressure_condition(up_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<DownOutflowPressure>> down_inflow_pressure_condition(down_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");

    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_centerline_velocity("Velocity", velocity_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    up_bidirection_buffer.tag_buffer_particles.exec();
    down_bidirection_buffer.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 30.0;                /**< End time. */
    Real Output_Time = end_time / 300.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                       /**< Default acoustic time step sizes. */
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
    body_states_recording.writeToFile();
    write_centerline_velocity.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                up_inflow_pressure_condition.exec(dt);
                down_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_bidirection_buffer.injection.exec();
            up_bidirection_buffer.injection.exec();
            down_bidirection_buffer.injection.exec();

            left_bidirection_buffer.deletion.exec();
            up_bidirection_buffer.deletion.exec();
            down_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            up_bidirection_buffer.tag_buffer_particles.exec();
            down_bidirection_buffer.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        velocity_observer_contact.updateConfiguration();
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
