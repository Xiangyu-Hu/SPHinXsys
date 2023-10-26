/**
 * @file 	elastic_tube.cpp
 * @brief 	3D pusatile flow in an elastic tube
 * @details This is the one of the basic test cases for validating fluid-solid interaction
 * @author
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include <gtest/gtest.h>

#include "base_data_type.h"
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Namespace cite here.
 */
using namespace SPH;
// unit change
const Real unit_length = 0.01;   // cm to m
const Real unit_force = 1e-5;    // dynes to N
const Real unit_viscosity = 0.1; // poise to 0.1 Pa.s
const Real unit_density = 1e3;   // g.cm-3 to kg.m-3
// time cycle
const int number_of_cycle = 2;
const Real time_cycle = 3e-3;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real fluid_radius = 0.5 * unit_length;
const Real diameter = 2.0 * fluid_radius;
const Real wall_thickness = 0.1 * unit_length;
const Real wall_radius = fluid_radius + wall_thickness;
const int number_of_particles = 40;
const Real resolution_fluid = diameter / Real(number_of_particles);
const Real emitter_length = resolution_fluid * 4.0;
const Real inflow_length = resolution_fluid * 10.0; // Region to impose velocity
const Real resolution_solid = wall_thickness / 4.0; // 4 layers of wall particles
const Real full_length = 5 * unit_length;
const int simtk_resolution = 20;
const Vec3d translation_fluid(0., full_length * 0.5, 0.);
/**
 * @brief Geometry parameters for boundary condition.
 */
const Vec3d emitter_halfsize(fluid_radius, emitter_length * 0.5, fluid_radius);
const Vec3d emitter_translation(0., emitter_length * 0.5, 0.);
const Vec3d buffer_halfsize(fluid_radius, inflow_length * 0.5, fluid_radius);
const Vec3d buffer_translation(0., inflow_length * 0.5, 0.);
const Vec3d disposer_halfsize(fluid_radius * 1.1, resolution_fluid * 2, fluid_radius * 1.1);
const Vec3d disposer_translation(0., full_length - disposer_halfsize[1], 0.);

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-diameter, 0, -diameter) - Vec3d(wall_thickness, wall_thickness, wall_thickness),
                                 Vec3d(diameter, full_length, diameter) + Vec3d(wall_thickness, wall_thickness, wall_thickness));
/**
 * @brief Material properties of the fluid.
 */
const Real rho0_f = 1.0 * unit_density;  /**< Reference density of fluid. */
const Real mu_f = 0.03 * unit_viscosity; /**< Viscosity. */
const Real U_f = 1.0;
const Real U_max = 5.0 * U_f;  // parabolic inflow, Thus U_max = 2*U_f
const Real c_f = 10.0 * U_max; /**< Reference sound speed. */
/**
 * @brief Material properties of the solid.
 */
const Real rho0_s = 1.0 * unit_density; /**< Reference density.*/
const Real poisson = 0.3;               /**< Poisson ratio.*/
const Real Youngs_modulus = 3e6 * unit_force / unit_length / unit_length;
/**
 * @brief Define water shape
 */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius, full_length * 0.5, simtk_resolution, translation_fluid);
    }
};
/**
 * @brief Define wall shape
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), wall_radius, full_length * 0.5 + wall_thickness, simtk_resolution, translation_fluid);
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(0., 1., 0.), fluid_radius, full_length * 0.5 + wall_thickness, simtk_resolution, translation_fluid);
    }
};
/** create the wall constrain shape. */
class FixedShape : public ComplexShape
{
  public:
    explicit FixedShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d inlet_fixation_halfsize(wall_radius, 0.5 * (wall_thickness + inflow_length), wall_radius);
        Vec3d inlet_fixation_translation(0, inflow_length - inlet_fixation_halfsize[1], 0);

        Vec3d outlet_fixation_halfsize(wall_radius, wall_thickness, wall_radius);
        Vec3d outlet_fixation_translation(0, full_length, 0);

        add<TriangleMeshShapeBrick>(inlet_fixation_halfsize, simtk_resolution, inlet_fixation_translation);
        add<TriangleMeshShapeBrick>(outlet_fixation_halfsize, simtk_resolution, outlet_fixation_translation);
    }
};

/**
 * @brief Inflow velocity
 */
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    int number_of_cycle_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(time_cycle), number_of_cycle_(number_of_cycle),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vec3d(0, 0, 0);
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < number_of_cycle_ * t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / (0.5 * t_ref_))) : 0;
        target_velocity[1] = 2.0 * u_ave * (1.0 - (position[0] * position[0] + position[2] * position[2]) / fluid_radius / fluid_radius);
        return target_velocity;
    }
};

/**
 * @brief 	Main program starts here.
 */
int main()
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem system(system_domain_bounds, resolution_fluid);
    system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    system.setReloadParticles(true);        // Tag for computation with save particles distribution
    IOEnvironment io_environment(system);
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineAdaptation<SPH::SPHAdaptation>(1.15, resolution_fluid / resolution_solid);
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    wall_boundary.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&wall_boundary});
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_wall_boundary_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
    }
    /** topology */
    InnerRelation wall_boundary_inner(wall_boundary);
    ComplexRelation water_block_complex(water_block, {&wall_boundary});
    ContactRelation wall_boundary_contact(wall_boundary, {&water_block});
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
    /** time-space method to detect surface particles. (Important for DensitySummationFreeSurfaceComplex work correctly.)*/
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex> inlet_outlet_surface_particle_indicator(water_block_complex);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    /** Impose transport velocity. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_complex);
    /**
     * @brief 	Boundary conditions. Inflow & Outflow in Y-direction
     */
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec3d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, yAxis);
    BodyAlignedBoxByCell inflow_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec3d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inflow_buffer);
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec3d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, yAxis);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Corrected configuration for the elastic insert body. */
    InteractionWithUpdate<CorrectedConfigurationInner> wall_boundary_corrected_configuration(wall_boundary_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(wall_boundary_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid>
        fluid_force_on_solid_update(wall_boundary_contact, viscous_force_on_solid);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(wall_boundary);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    /** Compute time step size of elastic solid. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> wall_boundary_computing_time_step_size(wall_boundary);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> wall_boundary_stress_relaxation_first_half(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_boundary_stress_relaxation_second_half(wall_boundary_inner);
    /** Constrain region of the inserted body. */
    auto wall_fixed_shape = makeShared<FixedShape>("wall_fixation_shape");
    BodyRegionByParticle wall_fixed(wall_boundary, wall_fixed_shape);
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constrain_wall(wall_fixed);
    /** Update norm .*/
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_boundary_update_normal(wall_boundary);
    // Damping the velocity field for quasi-static solution
    const Real shape_constant = 0.4;
    const Real physical_viscosity = shape_constant / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * wall_thickness;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        wall_damping(0.2, wall_boundary_inner, "Velocity", physical_viscosity);
    /**
     * @brief Output.
     */
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<Real>("MassiveMeasure");
    water_block.addBodyStateForRecording<Real>("VolumetricMeasure");
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(io_environment, system.real_bodies_);
    /**
     * @brief Setup geometry and initial conditions.
     */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    wall_boundary_corrected_configuration.exec();
    /** Output the start states of bodies. */
    body_states_recording.writeToFile(0);
    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 5;
    Real end_time = time_cycle * (number_of_cycle + 1); /**< End time. */
    Real Output_Time = end_time / 100;                  /**< Time stamps for output of body states. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * @brief 	Main loop starts here.
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_solid.exec();
            /** Update normal direction on elastic body.*/
            wall_boundary_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                fluid_force_on_solid_update.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(wall_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    wall_boundary_stress_relaxation_first_half.exec(dt_s);
                    constrain_wall.exec();
                    wall_damping.exec(dt_s);
                    constrain_wall.exec();
                    wall_boundary_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inflow_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;

            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            wall_boundary.updateCellLinkedList();

            water_block_complex.updateConfiguration();
            wall_boundary_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
