/**
 * @file 	channel_flow_shell.cpp
 * @brief 	This is a basic test of fluid-rigid shell interaction.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real DL = 10.0; /**< Channel length. */
const Real DH = 2.0;  /**< Channel height. */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1.0;                  /**< Density. */
const Real U_f = 1.0;                     /**< Characteristic velocity. */
const Real c_f = 10.0 * U_f;              /**< Speed of sound. */
const Real Re = 100.0;                    /**< Reynolds number. */
const Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
namespace SPH
{
class WaterBlock : public MultiPolygonShape
{
  public:
    class FluidAxialObserver;
    explicit WaterBlock(const std::vector<Vecd> &shape, const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(shape, ShapeBooleanOps::add);
    }
};

/** Particle generator and constraint boundary for shell baffle. */
class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real DL_sponge_;
    Real BW_;
    Real resolution_ref_;
    Real wall_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real resolution_ref, Real wall_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          DL_sponge_(20 * resolution_ref), BW_(4 * resolution_ref),
          resolution_ref_(resolution_ref), wall_thickness_(wall_thickness) {};
    void prepareGeometricData() override
    {
        auto particle_number_mid_surface = int((DL + DL_sponge_ + 2 * BW_) / resolution_ref_);
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real x = -DL_sponge_ - BW_ + (Real(i) + 0.5) * resolution_ref_;
            // upper wall
            Real y1 = DH + 0.5 * resolution_ref_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_ref_);
            Vec2d normal_direction_1 = Vec2d(0, 1.0);
            addSurfaceProperties(normal_direction_1, wall_thickness_);
            // lower wall
            Real y2 = -0.5 * resolution_ref_; // lower wall
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_ref_);
            Vec2d normal_direction_2 = Vec2d(0, -1.0);
            addSurfaceProperties(normal_direction_2, wall_thickness_);
        }
    }
};

//----------------------------------------------------------------------
//	Inflow velocity
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
        if (aligned_box_.checkInBounds(position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};

StdVec<Vecd> createFluidAxialObservationPoints(Real resolution_ref)
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the entrance of the channel. */
    size_t number_observation_points = 51;
    Real range_of_measure = DL - resolution_ref * 4.0;
    Real start_of_measure = resolution_ref * 2.0;
    Real y_position = 0.5 * DH;
    /** the measuring locations */
    for (size_t i = 0; i < number_observation_points; ++i)
    {
        Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure, y_position);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};

StdVec<Vecd> createFluidRadialObservationPoints(Real resolution_ref)
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the entrance of the channel. */
    size_t number_observation_points = 21;
    Real range_of_measure = DH - resolution_ref * 4.0;
    Real start_of_measure = resolution_ref * 2.0;
    Real x_position = 0.5 * DL;
    /** the measuring locations */
    for (size_t i = 0; i < number_observation_points; ++i)
    {
        Vec2d point_coordinate(x_position, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
} // namespace SPH

void channel_flow_shell(const Real resolution_ref, const Real wall_thickness)
{
    Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
    Real BW = resolution_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */

    /** Domain bounds of the system. */
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -wall_thickness), Vec2d(DL + BW, DH + wall_thickness));
    //----------------------------------------------------------------------
    //	define geometry of SPH bodies
    //----------------------------------------------------------------------
    /** create a water block shape */
    auto createWaterBlockShape = [&]()
    {
        // geometry
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
        water_block_shape.push_back(Vecd(-DL_sponge, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

        return water_block_shape;
    };
    Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>(createWaterBlockShape(), "WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<SurfaceParticles, WallBoundary>(resolution_ref, wall_thickness);

    ObserverBody fluid_axial_observer(sph_system, "FluidAxialObserver");
    fluid_axial_observer.generateParticles<ObserverParticles>(createFluidAxialObservationPoints(resolution_ref));

    ObserverBody fluid_radial_observer(sph_system, "FluidRadialObserver");
    fluid_radial_observer.generateParticles<ObserverParticles>(createFluidRadialObservationPoints(resolution_ref));
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_boundary_inner(wall_boundary);
    ShellInnerRelationWithContactKernel shell_curvature_inner(wall_boundary, water_block);
    // shell normal should point from fluid to shell
    // normal corrector set to false if shell normal is already pointing from fluid to shell
    ContactRelationFromShellToFluid water_block_contact(water_block, {&wall_boundary}, {false});
    ContactRelation fluid_axial_observer_contact(fluid_axial_observer, {&water_block});
    ContactRelation fluid_radial_observer_contact(fluid_radial_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Pressure relaxation using Verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, 1.5 * U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Inflow boundary condition. */
    AlignedBoxByCell inflow_buffer(
        water_block, AlignedBox(xAxis, Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    /** Periodic BCs in x direction. */
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    // Curvature calculation
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(wall_boundary, "Average1stPrincipleCurvature");
    ObservedQuantityRecording<Vecd> write_fluid_axial_velocity("Velocity", fluid_axial_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_radial_velocity("Velocity", fluid_radial_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();

    /** initial curvature*/
    shell_curvature.exec();
    water_block_complex.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 10.0;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_correction.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** velocity */
                parabolic_inflow.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        write_real_body_states.writeToFile();
        fluid_axial_observer_contact.updateConfiguration();
        fluid_radial_observer_contact.updateConfiguration();
        write_fluid_axial_velocity.writeToFile(number_of_iterations);
        write_fluid_radial_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    /**
     * @brief 	Gtest start from here.
     */
    /* Define analytical solution of the inflow velocity.*/
    std::function<Vec2d(Vec2d)> inflow_velocity = [&](Vec2d pos)
    {
        Real y = 2 * pos[1] / DH - 1;
        return Vec2d(1.5 * U_f * (1 - y * y), 0);
    };
    /* Compare all simulation to the analytical solution. */
    // Axial direction.
    BaseParticles &fluid_axial_particles = fluid_axial_observer.getBaseParticles();
    Vecd *pos_axial = fluid_axial_particles.ParticlePositions();
    Vecd *vel_axial = fluid_axial_particles.getVariableDataByName<Vecd>("Velocity");
    for (size_t i = 0; i < fluid_axial_particles.TotalRealParticles(); i++)
    {
        EXPECT_NEAR(inflow_velocity(pos_axial[i])[1], vel_axial[i][1], U_f * 5e-2);
    }
    // Radial direction
    BaseParticles &fluid_radial_particles = fluid_radial_observer.getBaseParticles();
    Vecd *pos_radial = fluid_radial_particles.ParticlePositions();
    Vecd *vel_radial = fluid_radial_particles.getVariableDataByName<Vecd>("Velocity");
    for (size_t i = 0; i < fluid_radial_particles.TotalRealParticles(); i++)
    {
        EXPECT_NEAR(inflow_velocity(pos_radial[i])[1], vel_radial[i][1], U_f * 5e-2);
    }
}

TEST(channel_flow_shell, thickness_10x)
{ // for CI
    const Real resolution_ref = 0.05;
    const Real wall_thickness = 10 * resolution_ref;
    channel_flow_shell(resolution_ref, wall_thickness);
}

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
