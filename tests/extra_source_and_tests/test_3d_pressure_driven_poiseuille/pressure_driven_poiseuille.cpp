/**
 * @file 	pulsatile_poiseuille_flow.cpp
 * @brief 	2D pulsatile poiseuille flow example
 * @details This is the one of the basic test cases for pressure boundary condition and bidirectional buffer.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "base_particles.h"
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

/** Name space. */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real CR = 0.0005;               /**< Channel radius. */
Real CL = CR * 2 * 15;          /**< Channel length. */
Real resolution_ref = CR / 5.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0; /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(CL + BW, 2.0 * CR + BW, 2.0 * CR + BW));
/**
 * @brief Material properties of the fluid.
 */
Real Inlet_pressure = 10.0;
Real Outlet_pressure = 5.0;
Real rho0_f = 1000.0;
Real Re = 200.0;
Real mu_f = sqrt(2 * rho0_f * pow(CR, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (4.0 * Re * CL)); // using max velocity
Real U_f = pow(CR, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (4.0 * mu_f * CL);
Real c_f = 10.0 * U_f;
Vecd normal = Vecd(1.0, 0.0, 0.0);

int water_resolution(10);

Vecd bidirectional_buffer_halfsize = Vecd(2.5 * resolution_ref, 1.5 * CR, 1.5 * CR);
Vecd left_bidirectional_translation = bidirectional_buffer_halfsize;
Vecd right_bidirectional_translation = Vecd(CL - 2.5 * resolution_ref, 1.5 * CR, 1.5 * CR);

class LeftInflowPressure : public fluid_dynamics::FlowPressureBuffer
{
  public:
    LeftInflowPressure(BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector) {}
    Real getTargetPressure(size_t index_i, Real dt) override
    {
        //        Real pressure = Inlet_pressure * cos(GlobalStaticVariables::physical_time_);
        Real pressure = Inlet_pressure;
        return pressure;
    }
    void setupDynamics(Real dt = 0.0) override {}
};

class RightInflowPressure : public fluid_dynamics::FlowPressureBuffer
{
  public:
    RightInflowPressure(BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector) {}
    Real getTargetPressure(size_t index_i, Real dt) override
    {
        Real pressure = Outlet_pressure;
        return pressure;
    }
    void setupDynamics(Real dt = 0.0) override {}
};

class LeftBidirectionalBufferCondition : public fluid_dynamics::BidirectionalBuffer
{
  public:
    LeftBidirectionalBufferCondition(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                     size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {

        /*constant pressure*/
        Real pressure = Inlet_pressure;
        return pressure;
    }
};

class RightBidirectionalBufferCondition : public fluid_dynamics::BidirectionalBuffer
{
  public:
    RightBidirectionalBufferCondition(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                      size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        Vecd translation_cylinder(0.5 * CL, CR, CR);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0, 0), CR,
                                       0.5 * CL, water_resolution, translation_cylinder);
    }
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        Vecd translation_cylinder(0.5 * CL, CR, CR);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0, 0), CR + BW,
                                       0.5 * CL + BW, water_resolution, translation_cylinder);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0, 0), CR,
                                            0.5 * CL + 2.0 * BW, water_resolution, translation_cylinder);
    }
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    std::cout << " Re =" << rho0_f * U_f * CR * 2 / mu_f << ", Re prescribed= " << Re << std::endl;
    std::cout << " U_f =" << U_f << std::endl;

    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.defineBodyLevelSetShape();
    // water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.defineBodyLevelSetShape();
    // wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});

    ComplexRelation water_block_complex(water_block_inner, water_block_contact);

    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    /**
     * @brief 	Methods used for time stepping.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    BodyAlignedBoxByCell left_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vecd(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);

    LeftBidirectionalBufferCondition left_emitter_inflow_injection(
        water_block, makeShared<AlignedBoxShape>(Transform(Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize), 10, xAxis);
    RightBidirectionalBufferCondition right_emitter_inflow_injection(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(right_bidirectional_translation)), bidirectional_buffer_halfsize), 10, xAxis);

    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("PositionDivergence");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);

    BodyRegionByCell left_pressure_region(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<LeftInflowPressure> left_pressure_condition(left_pressure_region, normal);
    BodyRegionByCell right_pressure_region(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(Vecd(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<RightInflowPressure> right_pressure_condition(right_pressure_region, normal);

    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact); /**
                                                                                * @brief Output.
                                                                                */
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.sph_bodies_);
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
    Real end_time = 5.0;    /**< End time. */
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
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                pressure_relaxation.exec(dt);
                kernel_summation.exec();
                left_pressure_condition.exec(dt);
                right_pressure_condition.exec(dt);
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
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
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

    return 0;
}
