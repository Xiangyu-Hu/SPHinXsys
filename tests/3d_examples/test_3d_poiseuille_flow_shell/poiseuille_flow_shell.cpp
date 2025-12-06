/**
 * @file 	poiseuille_flow_shell.cpp
 * @brief 	3D poiseuille flow interaction with shell example
 * @details This is the one of the basic test cases for validating fluid-rigid shell interaction
 * @author  Weiyi Kong
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real diameter = 6.35 * scale;
const Real fluid_radius = 0.5 * diameter;
const Real full_length = 10 * fluid_radius;
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
const Real rho0_f = 1050.0; /**< Reference density of fluid. */
const Real mu_f = 3.6e-3;   /**< Viscosity. */
const Real Re = 100;
/**< Characteristic velocity. Average velocity */
const Real U_f = Re * mu_f / rho0_f / diameter;
const Real U_max = 2.0 * U_f;  // parabolic inflow, Thus U_max = 2*U_f
const Real c_f = 10.0 * U_max; /**< Reference sound speed. */

namespace SPH
{
StdVec<Vecd> createAxialObservationPoints(
    Real full_length, Vec3d translation = Vec3d(0.0, 0.0, 0.0))
{
    StdVec<Vecd> observation_points;
    int ny = 51;
    for (int i = 0; i < ny; i++)
    {
        Real y = full_length / (ny - 1) * i;
        Vec3d point_coordinate(0.0, y, 0.0);
        observation_points.emplace_back(point_coordinate + translation);
    }
    return observation_points;
};

StdVec<Vecd> createRadialObservationPoints(
    Real full_length, Real diameter, int number_of_particles,
    Vec3d translation = Vec3d(0.0, 0.0, 0.0))
{
    StdVec<Vecd> observation_points;
    int n = number_of_particles + 1;
    Real y = full_length / 2.0;
    for (int i = 0; i < n - 1; i++) // we leave out the point close to the boundary as the
                                    // interpolation there is incorrect
                                    // TODO: fix the interpolation
    {
        Real z = diameter / 2.0 * i / Real(n);
        observation_points.emplace_back(Vec3d(0.0, y, z) + translation);
        observation_points.emplace_back(Vec3d(0.0, y, -z) + translation);
    }
    return observation_points;
};

class ShellBoundary;
template <>
class ParticleGenerator<SurfaceParticles, ShellBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real resolution_shell_;
    Real wall_thickness_;
    Real shell_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real resolution_shell, Real wall_thickness, Real shell_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          resolution_shell_(resolution_shell),
          wall_thickness_(wall_thickness), shell_thickness_(shell_thickness) {};
    void prepareGeometricData() override
    {
        const Real radius_mid_surface = fluid_radius + resolution_shell_ * 0.5;
        const auto particle_number_mid_surface =
            int(2.0 * radius_mid_surface * Pi / resolution_shell_);
        const auto particle_number_height =
            int((full_length + 2.0 * wall_thickness_) / resolution_shell_);
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            for (int j = 0; j < particle_number_height; j++)
            {
                Real theta = (i + 0.5) * 2 * Pi / (Real)particle_number_mid_surface;
                Real x = radius_mid_surface * cos(theta);
                Real y = -wall_thickness_ + (full_length + 2 * wall_thickness_) * j / (Real)particle_number_height + 0.5 * resolution_shell_;
                Real z = radius_mid_surface * sin(theta);
                addPositionAndVolumetricMeasure(Vec3d(x, y, z),
                                                resolution_shell_ * resolution_shell_);
                Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
                addSurfaceProperties(n_0, shell_thickness_);
            }
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
    Vec3d halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(1.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vec3d operator()(Vec3d &position, Vec3d &velocity, Real current_time)
    {
        Vec3d target_velocity = Vec3d(0, 0, 0);
        target_velocity[1] = SMAX(2.0 * U_f *
                                      (1.0 - (position[0] * position[0] + position[2] * position[2]) /
                                                 fluid_radius / fluid_radius),
                                  0.0);
        return target_velocity;
    }
};
} // namespace SPH

//----------------------------------------------------------------------
//	Test case function
//----------------------------------------------------------------------
void poiseuille_flow(const Real resolution_ref, const Real resolution_shell, const Real shell_thickness)
{
    //----------------------------------------------------------------------
    //	Geometry parameters for shell.
    //----------------------------------------------------------------------
    const int number_of_particles = 10;
    const Real inflow_length = resolution_ref * 10.0; // Inflow region
    const Real wall_thickness = resolution_ref * 4.0;
    const int SimTK_resolution = 20;
    const Vec3d translation_fluid(0., full_length * 0.5, 0.);

    //----------------------------------------------------------------------
    //	Geometry parameters for boundary condition.
    //----------------------------------------------------------------------
    const Vec3d emitter_halfsize(fluid_radius, resolution_ref * 2, fluid_radius);
    const Vec3d emitter_translation(0., resolution_ref * 2, 0.);
    const Vec3d emitter_buffer_halfsize(fluid_radius, inflow_length * 0.5, fluid_radius);
    const Vec3d emitter_buffer_translation(0., inflow_length * 0.5 - 2 * resolution_ref, 0.);
    const Vec3d disposer_halfsize(fluid_radius * 1.1, resolution_ref * 2, fluid_radius * 1.1);
    const Vec3d disposer_translation(0., full_length - disposer_halfsize[1], 0.);

    //----------------------------------------------------------------------
    //	Domain bounds of the system.
    //----------------------------------------------------------------------
    const BoundingBoxd system_domain_bounds(Vec3d(-0.5 * diameter, 0, -0.5 * diameter) -
                                               Vec3d(shell_thickness, wall_thickness,
                                                     shell_thickness),
                                           Vec3d(0.5 * diameter, full_length, 0.5 * diameter) +
                                               Vec3d(shell_thickness, wall_thickness,
                                                     shell_thickness));
    //----------------------------------------------------------------------
    //  Define water shape
    //----------------------------------------------------------------------
    auto water_block_shape = makeShared<ComplexShape>("WaterBody");
    water_block_shape->add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius,
                                                      full_length * 0.5, SimTK_resolution,
                                                      translation_fluid);

    //----------------------------------------------------------------------
    //  Build up -- a SPHSystem --
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, water_block_shape);
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody shell_boundary(system, makeShared<DefaultShape>("Shell"));
    shell_boundary.defineAdaptation<SPH::SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell_boundary.defineMaterial<LinearElasticSolid>(1, 1e3, 0.45);
    shell_boundary.generateParticles<SurfaceParticles, ShellBoundary>(resolution_shell, wall_thickness, shell_thickness);

    ObserverBody observer_axial(system, "fluid_observer_axial");
    observer_axial.generateParticles<ObserverParticles>(createAxialObservationPoints(full_length));
    ObserverBody observer_radial(system, "fluid_observer_radial");
    observer_radial.generateParticles<ObserverParticles>(createRadialObservationPoints(full_length, diameter, number_of_particles));
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_boundary_inner(shell_boundary);
    ShellInnerRelationWithContactKernel wall_curvature_inner(shell_boundary, water_block);
    // shell normal should point from fluid to shell
    // normal corrector set to false if shell normal is already pointing from fluid to shell
    ContactRelationFromShellToFluid water_block_contact(water_block, {&shell_boundary}, {false});
    ContactRelation observer_contact_axial(observer_axial, {&water_block});
    ContactRelation observer_contact_radial(observer_radial, {&water_block});
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
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> wall_corrected_configuration(shell_boundary_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(wall_curvature_inner);

    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_block_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    //	Boundary conditions. Inflow & Outflow in Y-direction
    //----------------------------------------------------------------------
    AlignedBoxByParticle emitter(water_block, AlignedBox(yAxis, Transform(Vec3d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);
    AlignedBoxByCell emitter_buffer(water_block, AlignedBox(yAxis, Transform(Vec3d(emitter_buffer_translation)), emitter_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);
    AlignedBoxByCell disposer(water_block, AlignedBox(yAxis, Transform(Vec3d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(system);
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(shell_boundary, "Average1stPrincipleCurvature");
    body_states_recording.addToWrite<Real>(shell_boundary, "Average2ndPrincipleCurvature");
    ObservedQuantityRecording<Vec3d> write_fluid_velocity_axial("Velocity", observer_contact_axial);
    ObservedQuantityRecording<Vec3d> write_fluid_velocity_radial("Velocity", observer_contact_radial);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_corrected_configuration.exec();
    shell_curvature.exec();
    water_block_complex.updateConfiguration();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 2.0;               /**< End time. */
    Real Output_Time = end_time / 100; /**< Time stamps for output of body states. */
    Real dt = 0.0;                     /**< Default acoustic time step sizes. */

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
    body_states_recording.writeToFile(0);

    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            /** Acceleration due to viscous force and gravity. */
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(),
                          Dt - relaxation_time);
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                emitter_buffer_inflow_condition.exec();
            }
            interval_computing_pressure_relaxation +=
                TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << "	Time = " << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;

        /** Update observer and write output of observer. */
        observer_contact_axial.updateConfiguration();
        observer_contact_radial.updateConfiguration();
        write_fluid_velocity_axial.writeToFile(number_of_iterations);
        write_fluid_velocity_radial.writeToFile(number_of_iterations);
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9)
              << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    //----------------------------------------------------------------------
    //	Gtest starts from here.
    //----------------------------------------------------------------------
    /* Define analytical solution of the inflow velocity.*/
    std::function<Vec3d(Vec3d)> inflow_velocity = [&](Vec3d pos)
    {
        return Vec3d(0.0,
                     2.0 * U_f * (1.0 - (pos[0] * pos[0] + pos[2] * pos[2]) / (diameter * 0.5) / (diameter * 0.5)),
                     0.0);
    };
    /* Compare all simulation to the analytical solution. */
    // Axial direction.
    BaseParticles &axial_particles = observer_axial.getBaseParticles();
    Vecd *pos_axial = axial_particles.ParticlePositions();
    Vecd *vel_axial = axial_particles.getVariableDataByName<Vecd>("Velocity");
    for (size_t i = 0; i < axial_particles.TotalRealParticles(); i++)
    {
        EXPECT_NEAR(inflow_velocity(pos_axial[i])[1], vel_axial[i][1], U_max * 10e-2); // it's below 5% but 10% for CI
    }
    // Radial direction
    BaseParticles &radial_particles = observer_axial.getBaseParticles();
    Vecd *pos_radial = radial_particles.ParticlePositions();
    Vecd *vel_radial = radial_particles.getVariableDataByName<Vecd>("Velocity");
    for (size_t i = 0; i < radial_particles.TotalRealParticles(); i++)
    {
        EXPECT_NEAR(inflow_velocity(pos_radial[i])[1], vel_radial[i][1], U_max * 10e-2); // it's below 5% but 10% for CI
    }
}

TEST(poiseuille_flow, 10_particles)
{ // for CI
    const int number_of_particles = 10;
    const Real resolution_ref = diameter / number_of_particles;
    const Real resolution_shell = 0.5 * resolution_ref;
    const Real shell_thickness = 0.5 * resolution_shell;
    poiseuille_flow(resolution_ref, resolution_shell, shell_thickness);
}

TEST(DISABLED_poiseuille_flow, 20_particles)
{ // for CI
    const int number_of_particles = 20;
    const Real resolution_ref = diameter / number_of_particles;
    const Real resolution_shell = 0.5 * resolution_ref;
    const Real shell_thickness = 0.5 * resolution_shell;
    poiseuille_flow(resolution_ref, resolution_shell, shell_thickness);
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
