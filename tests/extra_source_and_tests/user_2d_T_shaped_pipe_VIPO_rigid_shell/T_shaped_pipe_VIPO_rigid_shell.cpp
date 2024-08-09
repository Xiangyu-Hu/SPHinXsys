/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.2;                                               /**< Reference length. */
Real DH = 0.1;                                               /**< Reference and the height of main channel. */
Real DL1 = 0.75 * DL;                                        /**< The length of the main channel. */
Real resolution_ref = 0.002;                                 /**< Initial reference particle spacing. */
Real BW = resolution_ref * 1.0;                                
Real buffer_width = resolution_ref * 4.0;                                /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;                        /**< Reference size of the emitter buffer to impose inflow condition. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * BW);
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real Outlet_pressure = 0.1;
Real rho0_f = 1060;                                                 /**< Reference density of fluid. */
Real Re = 100.0;                                                      /**< Reynolds number. */
Real U_f = 0.5;                                                      /**< Characteristic velocity. */
Real mu_f = 0.00355; /**< Dynamics viscosity. */
Real c_f =10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * (DL - DL1)));
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

class ShellShape : public MultiPolygonShape
{
  public:
    explicit ShellShape(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

/** Particle generator and constraint boundary for shell baffle. */
class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real DL_sponge_;
    Real resolution_ref_;
    Real wall_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real resolution_ref, Real wall_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          DL_sponge_(20 * resolution_ref),
          resolution_ref_(resolution_ref), wall_thickness_(wall_thickness){};
    void prepareGeometricData() override
    {
        auto particle_number_mid_surface_01 = int((DL1 + DL_sponge_) / resolution_ref_);
        //std::cout << " particle_number_mid_surface_01 = " << particle_number_mid_surface_01 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_01 - 1; i++)
        {
            Real x = -DL_sponge_ + (Real(i) + 0.5) * resolution_ref_;
            // upper wall
            Real y1 = DH + 0.5 * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_ref_);
            Vec2d normal_direction_1 = Vec2d(0, 1.0);
            addSurfaceProperties(normal_direction_1, wall_thickness_);
            // lower wall
            Real y2 = - 0.5 * resolution_ref_; // lower wall
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_ref_);
            Vec2d normal_direction_2 = Vec2d(0, -1.0);
            addSurfaceProperties(normal_direction_2, wall_thickness_);
        }

        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_ref_, DH + 0.5 * resolution_ref_), resolution_ref_);
        addSurfaceProperties(Vec2d(-1.0, 1.0).normalized(), wall_thickness_);
        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_ref_, - 0.5 * resolution_ref_), resolution_ref_);
        addSurfaceProperties(Vec2d(-1.0, -1.0).normalized(), wall_thickness_);

        auto particle_number_mid_surface_02 = int(DH / resolution_ref_);
        //std::cout << " particle_number_mid_surface_02 = " << particle_number_mid_surface_02 << std::endl;
        for (int i = 1; i < particle_number_mid_surface_02; i++)
        {
            // upper wall
            Real y1 = DH + (Real(i) + 0.5) * resolution_ref_;
            Real x = DL1 - 0.5 * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_ref_);
            Vec2d normal_direction = Vec2d(-1.0, 0);
            addSurfaceProperties(normal_direction, wall_thickness_);
            // lower wall
            Real y2 =  -(Real(i) + 0.5) * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_ref_);
            addSurfaceProperties(normal_direction, wall_thickness_);
        }

        auto particle_number_mid_surface_03 = int(1.5 * DH / resolution_ref_);
        //std::cout << " particle_number_mid_surface_03 = " << particle_number_mid_surface_03 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_03; i++)
        {
            Real y1 = 0.5 * DH + (Real(i) + 0.5) * resolution_ref_;
            // upper wall
            Real x = DL + 0.5 * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_ref_);
            Vec2d normal_direction = Vec2d(1.0, 0);
            addSurfaceProperties(normal_direction, wall_thickness_);
            // lower wall
            Real y2 = 0.5 * DH - (Real(i) + 0.5) * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_ref_);
            addSurfaceProperties(normal_direction, wall_thickness_);
        }
    }
};

//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
            aligned_box_(boundary_condition.getAlignedBox()){}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;

        return target_velocity;
    }
};

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

struct UpOutflowPressure
{
    template <class BoundaryConditionType>
    UpOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
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

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -DH - BW), Vec2d(DL + BW, 2.0 * DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment(); // handle command line arguments
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody shell_body(sph_system, makeShared<ShellShape>("ShellBody"));
    shell_body.defineAdaptation<SPHAdaptation>(1.15, 2.0);
    shell_body.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(sph_system);
    shell_body.defineMaterial<Solid>();
    shell_body.generateParticles<SurfaceParticles, WallBoundary>(resolution_ref, BW);

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ObserverParticles>(observer_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_inner(shell_body);
    ContactRelationFromShellToFluid water_shell_contact(water_block, {&shell_body}, {false});
    //ContactRelationFromFluidToShell shell_water_contact(shell_body, {&water_block}, {false});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell_body, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_shell_contact});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    // shell dynamics
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);

    // fluid dynamics
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_shell_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_shell_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_shell_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_shell_contact);

    //----------------------------------------------------------------------
    // Left buffer
    //----------------------------------------------------------------------
    Vec2d left_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.5 * DH);
    Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize;
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer);

    BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(Pi), Vec2d(left_buffer_translation)), left_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer);
    //----------------------------------------------------------------------
    // Up buffer
    //----------------------------------------------------------------------
    Vec2d up_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.75);
    Vec2d up_buffer_translation = Vec2d(0.5 * (DL + DL1), 2.0 * DH - 0.5 * buffer_width);
    BodyAlignedBoxByCell up_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(-0.5 * Pi), Vec2d(up_buffer_translation)), up_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<UpOutflowPressure> right_up_emitter_inflow_injection(up_emitter, in_outlet_particle_buffer);

    BodyAlignedBoxByCell up_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(0.5 * Pi), Vec2d(up_buffer_translation)), up_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_up_disposer_outflow_deletion(up_disposer);
    //----------------------------------------------------------------------
    // Down buffer
    //----------------------------------------------------------------------
    Vec2d down_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.75);
    Vec2d down_buffer_translation = Vec2d(0.5 * (DL + DL1), -DH + 0.5 * buffer_width);
    BodyAlignedBoxByCell down_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(0.5 * Pi), Vec2d(down_buffer_translation)), down_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<DownOutflowPressure> right_down_emitter_inflow_injection(down_emitter, in_outlet_particle_buffer);

    BodyAlignedBoxByCell down_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(-0.5 * Pi), Vec2d(down_buffer_translation)), down_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_down_disposer_outflow_deletion(down_disposer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_shell_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<UpOutflowPressure>> up_inflow_pressure_condition(up_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<DownOutflowPressure>> down_inflow_pressure_condition(down_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
    body_states_recording.addToWrite<Vecd>(shell_body, "NormalDirection");
    body_states_recording.addToWrite<Real>(shell_body, "Average1stPrincipleCurvature");
    body_states_recording.addToWrite<Real>(shell_body, "Average2ndPrincipleCurvature");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    shell_average_curvature.exec();
    water_block_complex.updateConfiguration();
    boundary_indicator.exec();
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_up_emitter_inflow_injection.tag_buffer_particles.exec();
    right_down_emitter_inflow_injection.tag_buffer_particles.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
     size_t number_of_iterations = sph_system.RestartStep();
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
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
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
                up_inflow_pressure_condition.exec(dt);
                down_inflow_pressure_condition.exec(dt);
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
                body_states_recording.writeToFile();
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_up_emitter_inflow_injection.injection.exec();
            right_down_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_up_disposer_outflow_deletion.exec();
            right_down_disposer_outflow_deletion.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();

            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_up_emitter_inflow_injection.tag_buffer_particles.exec();
            right_down_emitter_inflow_injection.tag_buffer_particles.exec();
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
