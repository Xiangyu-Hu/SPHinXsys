/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 0.004;                                             /**< Channel length. */
Real DH = 0.001;                                             /**< Channel height. */
Real resolution_ref = DH / 20.0;                             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.1;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
Real c_f = 10.0 * U_f;
/**
 * @brief buffer parameters.
 */
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL - 2.5 * resolution_ref, 0.5 * DH);
Vec2d normal = Vec2d(1.0, 0.0);

/**
 * @brief 	Pressure boundary definition.
 */
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

/**
 * @brief 	inflow velocity definition.
 */
struct InflowVelocity
{
    Real u_ave;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ave(0.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;

        u_ave = fabs(Inlet_pressure - Outlet_pressure) * (position[1] + 0.5 * DH) * (position[1] + 0.5 * DH - DH) / (2.0 * mu_f * DL) +
                (4.0 * fabs(Inlet_pressure - Outlet_pressure) * DH * DH) /
                    (mu_f * DL * Pi * Pi * Pi) * sin(Pi * (position[1] + 0.5 * DH) / DH) * exp(-(Pi * Pi * mu_f * run_time) / (DH * DH));

        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;

        return target_velocity;
    }
};

/**
 * @brief 	Fluid body definition.
 */
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
/**
 * @brief 	Wall boundary body definition.
 */
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

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setGenerateRegressionData(false);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<Lattice>(in_outlet_particle_buffer);
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<Lattice>();

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<Observer>(observer_location);
    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** delete outflow particles */
    BodyAlignedBoxByCell left_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        boundary_indicator(water_block_inner, water_block_contact);
    /** bidirectional buffer */
    BodyAlignedBoxByCell left_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer, xAxis);
    BodyAlignedBoxByCell right_emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_emitter_inflow_injection(right_emitter, in_outlet_particle_buffer, xAxis);
    /** output parameters */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** momentum equation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** mass equation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** pressure boundary condition. */
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_disposer);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact);

    /**
     * @brief Output.
     */
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_centerline_velocity("Velocity", velocity_observer_contact);
    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();

    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 10.0;   /**< End time. */
    Real Output_Time = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    write_centerline_velocity.writeToFile(number_of_iterations);
    /**
     * @brief 	Main loop starts here.
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
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
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerline_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();
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
