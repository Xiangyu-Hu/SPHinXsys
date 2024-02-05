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
Real DL = 10.0;                            /**< Channel length. */
Real DH = 2.0;                             /**< Channel height. */
Real resolution_ref = 0.05;                /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0;    /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;            /**< Boundary width, determined by specific layer of boundary particles. */
Real wall_thickness = 10 * resolution_ref; /*<Thickness of wall boundary, same as global resolution>*/

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -wall_thickness), Vec2d(DL + BW, DH + wall_thickness));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                  /**< Density. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Speed of sound. */
Real Re = 100.0;                    /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int((DL + DL_sponge + 2 * BW) / resolution_ref);
class WallBoundaryParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit WallBoundaryParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real x = -DL_sponge - BW + (Real(i) + 0.5) * resolution_ref;
            // upper wall
            Real y1 = DH + 0.5 * resolution_ref;
            initializePositionAndVolumetricMeasure(Vecd(x, y1), resolution_ref);
            Vec2d normal_direction_1 = Vec2d(0, 1.0);
            initializeSurfaceProperties(normal_direction_1, wall_thickness);
            // lower wall
            Real y2 = -0.5 * resolution_ref; // lower wall
            initializePositionAndVolumetricMeasure(Vecd(x, y2), resolution_ref);
            Vec2d normal_direction_2 = Vec2d(0, -1.0);
            initializeSurfaceProperties(normal_direction_2, wall_thickness);
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
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        if (aligned_box_.checkInBounds(0, position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};
/** fluid observer particle generator */
class FluidAxialObserverParticleGenerator : public ParticleGeneratorObserver
{
  public:
    explicit FluidAxialObserverParticleGenerator(SPHBody &sph_body) : ParticleGeneratorObserver(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 51;
        Real range_of_measure = DL - resolution_ref * 4.0;
        Real start_of_measure = resolution_ref * 2.0;
        Real y_position = 0.5 * DH;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure, y_position);
            positions_.push_back(point_coordinate);
        }
    }
};
class FluidRadialObserverParticleGenerator : public ParticleGeneratorObserver
{
  public:
    explicit FluidRadialObserverParticleGenerator(SPHBody &sph_body) : ParticleGeneratorObserver(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 21;
        Real range_of_measure = DH - resolution_ref * 4.0;
        Real start_of_measure = resolution_ref * 2.0;
        Real x_position = 0.5 * DL;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vec2d point_coordinate(x_position, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0); // dummy material parameters
    wall_boundary.generateParticles<WallBoundaryParticleGenerator>();

    ObserverBody fluid_axial_observer(sph_system, "FluidAxialObserver");
    fluid_axial_observer.generateParticles<FluidAxialObserverParticleGenerator>();

    ObserverBody fluid_radial_observer(sph_system, "FluidRadialObserver");
    fluid_radial_observer.generateParticles<FluidRadialObserverParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_boundary_inner(wall_boundary);
    ShellInnerRelationWithContactKernel shell_curvature_inner(wall_boundary, water_block);
    ContactRelationToShell water_block_contact(water_block, {&wall_boundary});
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
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, 1.5 * U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Inflow boundary condition. */
    BodyAlignedBoxByCell inflow_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, water_block.getBodyShapeBounds(), xAxis);
    // Curvature calculation
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    wall_boundary.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
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
    while (GlobalStaticVariables::physical_time_ < end_time)
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
                GlobalStaticVariables::physical_time_ += dt;

                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
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
     * @brief 	Gtest strat from here.
     */
    /* Define analytical solution of the inflow velocity.*/
    std::function<Vec2d(Vec2d)> inflow_velocity = [&](Vec2d pos)
    {
        Real y = 2 * pos[1] / DH - 1;
        return Vec2d(1.5 * U_f * (1 - y * y), 0);
    };
    /* Compare all simulation to the anayltical solution. */
    // Axial direction.
    for (size_t i = 0; i < fluid_axial_observer.getBaseParticles().pos_.size(); i++)
    {
        EXPECT_NEAR(inflow_velocity(fluid_axial_observer.getBaseParticles().pos_[i])[1],
                    fluid_axial_observer.getBaseParticles().vel_[i][1],
                    U_f * 5e-2);
    }
    // Radial direction
    for (size_t i = 0; i < fluid_radial_observer.getBaseParticles().pos_.size(); i++)
    {
        EXPECT_NEAR(inflow_velocity(fluid_radial_observer.getBaseParticles().pos_[i])[1],
                    fluid_radial_observer.getBaseParticles().vel_[i][1],
                    U_f * 5e-2);
    }

    return 0;
}
