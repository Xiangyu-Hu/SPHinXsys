/**
 * @file 	T_pipe_VIPO_shell.cpp
 * @brief 	This is a test of fluid shell interaction in a T-shaped pipe.
 * @details The flow with velocity inlet and pressure outlet boundary condition in a T-shaped pipe in 2D is considered.
 * @author 	Chenxi Zhao, Weiyi Kong, Xiangyu Hu
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
Real DL = 0.2;               /**< Reference length. */
Real DH = 0.1;               /**< Reference and the height of main channel. */
Real DL1 = 0.75 * DL;        /**< The length of the main channel. */
Real resolution_ref = 0.005; /**< Initial reference particle spacing. */
Real resolution_shell = resolution_ref;
Real BW = resolution_shell * 1.0;
Real buffer_width = resolution_ref * 4.0;                    /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;                        /**< Reference size of the emitter buffer to impose inflow condition. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * BW);
//----------------------------------------------------------------------
//	Global parameters on the fluid properties.
//----------------------------------------------------------------------
Real Outlet_pressure = 0;
Real rho0_f = 1.0;                                                    /**< Reference density of fluid. Using air density here.*/
Real Re = 100.0;                                                      /**< Reynolds number. */
Real U_f = 1.0;                                                       /**< Characteristic velocity. */
Real mu_f = rho0_f * U_f * DH / Re;                                   /**< Dynamics viscosity. */
Real c_f = 10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * (DL - DL1))); /** Reference sound speed needs to consider the flow speed in the narrow channels. */
//----------------------------------------------------------------------
//	Material parameters of the shell.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;         /** Normalized density. */
Real Youngs_modulus = 1.0e5; /** Normalized Youngs Modulus. */
Real poisson = 0.3;          /** Poisson ratio. */
//----------------------------------------------------------------------
//	Define geometry of SPH bodies.
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
namespace SPH
{
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
//----------------------------------------------------------------------
//	Particle generator of shell boundary.
//----------------------------------------------------------------------
class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real resolution_shell_;
    Real shell_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real resolution_shell, Real shell_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          resolution_shell_(resolution_shell), shell_thickness_(shell_thickness) {};
    void prepareGeometricData() override
    {
        auto particle_number_mid_surface_01 = int((DL1 + DL_sponge) / resolution_shell_);
        // std::cout << " particle_number_mid_surface_01 = " << particle_number_mid_surface_01 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_01 - 1; i++)
        {
            Real x = -DL_sponge + (Real(i) + 0.5) * resolution_shell_;
            // upper wall
            Real y1 = DH + 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction_1 = Vec2d(0, 1.0);
            addSurfaceProperties(normal_direction_1, shell_thickness_);
            // lower wall
            Real y2 = -0.5 * resolution_shell_; // lower wall
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            Vec2d normal_direction_2 = Vec2d(0, -1.0);
            addSurfaceProperties(normal_direction_2, shell_thickness_);
        }

        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_shell_, DH + 0.5 * resolution_shell_), resolution_shell_);
        addSurfaceProperties(Vec2d(-1.0, 1.0).normalized(), shell_thickness_);
        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_shell_, -0.5 * resolution_shell_), resolution_shell_);
        addSurfaceProperties(Vec2d(-1.0, -1.0).normalized(), shell_thickness_);

        auto particle_number_mid_surface_02 = int(DH / resolution_shell_);
        // std::cout << " particle_number_mid_surface_02 = " << particle_number_mid_surface_02 << std::endl;
        for (int i = 1; i < particle_number_mid_surface_02; i++)
        {
            // upper wall
            Real y1 = DH + (Real(i) + 0.5) * resolution_shell_;
            Real x = DL1 - 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction = Vec2d(-1.0, 0);
            addSurfaceProperties(normal_direction, shell_thickness_);
            // lower wall
            Real y2 = -(Real(i) + 0.5) * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            addSurfaceProperties(normal_direction, shell_thickness_);
        }

        auto particle_number_mid_surface_03 = int(1.5 * DH / resolution_shell_);
        // std::cout << " particle_number_mid_surface_03 = " << particle_number_mid_surface_03 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_03; i++)
        {
            Real y1 = 0.5 * DH + (Real(i) + 0.5) * resolution_shell_;
            // upper wall
            Real x = DL + 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction = Vec2d(1.0, 0);
            addSurfaceProperties(normal_direction, shell_thickness_);
            // lower wall
            Real y2 = 0.5 * DH - (Real(i) + 0.5) * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            addSurfaceProperties(normal_direction, shell_thickness_);
        }
    }
};
//----------------------------------------------------------------------
//	Inflow velocity.
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

    Vecd operator()(Vecd &position, Vecd &velocity, Real physical_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = physical_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * physical_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Pressure boundary condition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real physical_time)
    {
        return p;
    }
};

struct UpOutflowPressure
{
    template <class BoundaryConditionType>
    UpOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real physical_time)
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

    Real operator()(Real p, Real physical_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
//----------------------------------------------------------------------
//	Constrain region of shell body.
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        return base_particles_.ParticlePositions()[index_i][0] < -DL_sponge + buffer_width ||
               base_particles_.ParticlePositions()[index_i][1] > 2.0 * DH - buffer_width ||
               base_particles_.ParticlePositions()[index_i][1] < -DH + buffer_width;
    };
};
} // namespace SPH
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
    sph_system.setGenerateRegressionData(false);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody shell_body(sph_system, makeShared<ShellShape>("ShellBody"));
    shell_body.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell_body.defineBodyLevelSetShape(level_set_refinement_ratio, UsageType::Surface)
        ->writeLevelSet(sph_system);
    shell_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    shell_body.generateParticles<SurfaceParticles, WallBoundary>(resolution_shell, BW);

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
    ContactRelationFromFluidToShell shell_water_contact(shell_body, {&water_block}, {false});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell_body, water_block);
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_shell_contact});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    // Shell dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell_body);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell_body);
    /** Exert constrain on shell. */
    BoundaryGeometry boundary_geometry(shell_body);
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(boundary_geometry);
    //----------------------------------------------------------------------
    // Fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_shell_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_shell_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_shell_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_shell_contact);

    // left buffer
    Vec2d left_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.5 * DH);
    Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize;
    AlignedBoxByCell left_emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    // up buffer
    Vec2d up_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.75);
    Vec2d up_buffer_translation = Vec2d(0.5 * (DL + DL1), 2.0 * DH - 0.5 * buffer_width);
    AlignedBoxByCell up_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(-0.5 * Pi), Vec2d(up_buffer_translation)), up_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<UpOutflowPressure> right_up_bidirection_buffer(up_emitter, in_outlet_particle_buffer);
    // down buffer
    Vec2d down_buffer_halfsize = Vec2d(0.5 * buffer_width, 0.75);
    Vec2d down_buffer_translation = Vec2d(0.5 * (DL + DL1), -DH + 0.5 * buffer_width);
    AlignedBoxByCell down_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(0.5 * Pi), Vec2d(down_buffer_translation)), down_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<DownOutflowPressure> right_down_bidirection_buffer(down_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_shell_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<UpOutflowPressure>> up_inflow_pressure_condition(up_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<DownOutflowPressure>> down_inflow_pressure_condition(down_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    // FSI.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_on_shell(shell_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_body);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");
    body_states_recording.addToWrite<Vecd>(shell_body, "NormalDirection");
    body_states_recording.addToWrite<Vecd>(shell_body, "PressureForceFromFluid");
    body_states_recording.addToWrite<Real>(shell_body, "Average1stPrincipleCurvature");
    body_states_recording.addToWrite<Real>(shell_body, "Average2ndPrincipleCurvature");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_centerline_velocity("Velocity", velocity_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    shell_corrected_configuration.exec();
    shell_average_curvature.exec();
    constrain_holder.exec();
    water_block_complex.updateConfiguration();
    shell_water_contact.updateConfiguration();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_up_bidirection_buffer.tag_buffer_particles.exec();
    right_down_bidirection_buffer.tag_buffer_particles.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 15.0;              /**< End time. */
    Real Output_Time = end_time / 150; /**< Time stamps for output of body states. */
    Real dt = 0.0;                     /**< Default acoustic time step sizes. */
    Real dt_s = 0.0;                   /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
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
    while (sv_physical_time->getValue() < end_time)
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
            /** FSI for viscous force. */
            viscous_force_on_shell.exec();

            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_on_shell.exec();

                // boundary condition implementation
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                up_inflow_pressure_condition.exec(dt);
                down_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();

                density_relaxation.exec(dt);

                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    dt_s = shell_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    shell_stress_relaxation_first.exec(dt_s);

                    constrain_holder.exec(dt_s);

                    shell_stress_relaxation_second.exec(dt_s);
                    dt_s_sum += dt_s;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
                write_centerline_velocity.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_bidirection_buffer.injection.exec();
            right_up_bidirection_buffer.injection.exec();
            right_down_bidirection_buffer.injection.exec();

            left_bidirection_buffer.deletion.exec();
            right_up_bidirection_buffer.deletion.exec();
            right_down_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            shell_update_normal.exec();
            shell_body.updateCellLinkedList();
            shell_curvature_inner.updateConfiguration();
            shell_average_curvature.exec();
            shell_water_contact.updateConfiguration();
            water_block_complex.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();

            left_bidirection_buffer.tag_buffer_particles.exec();
            right_up_bidirection_buffer.tag_buffer_particles.exec();
            right_down_bidirection_buffer.tag_buffer_particles.exec();
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