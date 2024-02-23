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
/** Geometry parameters. */
Real scale = 0.001;
Real diameter = 6.35 * scale;
Real radius = diameter * 0.5;
Real full_length = 7.5 * diameter;
size_t number_of_particles = 10;
Real resolution_ref = diameter / number_of_particles;
Real wall_thickness = resolution_ref * 4;
auto inlet_normal = Vec3d::UnitX();
auto SimTK_inlet_normal = SimTK::UnitVec3(inlet_normal[0], inlet_normal[1], inlet_normal[2]);

//----------------------------------------------------------------------
//	Boundary related parameters
//----------------------------------------------------------------------
auto inlet_center = Vec3d::Zero();
auto outlet_center = inlet_center + full_length * inlet_normal;
size_t axial_direction = 0;
Real buffer_half_width = 3 * resolution_ref;
auto buffer_half_size = Vec3d(buffer_half_width, radius, radius);
auto left_bc_translation = inlet_center + inlet_normal * buffer_half_width;
auto left_bc_transform = Transform(left_bc_translation);

auto left_disposer_transform = Transform(Rotation3d(M_PI, Vec3d::UnitY()), left_bc_translation);

auto right_bc_translation = outlet_center - inlet_normal * buffer_half_width;
Rotation3d right_bc_rotation = Rotation3d(M_PI, Vec3d::UnitY());
auto right_bc_transform = Transform(right_bc_rotation, right_bc_translation);

auto right_disposer_translation = outlet_center - inlet_normal * buffer_half_width;
auto right_disposer_transform = Transform(right_disposer_translation);
Real compute_pressure(double p)
{
    // Real run_time = GlobalStaticVariables::physical_time_;
    // Real pressure = run_time < 0.25 ? 0.5 * p * (1.0 - cos(Pi * run_time / 0.25)) : p;
    Real pressure = p;
    return pressure;
};

//----------------------------------------------------------------------
//	Materials\ Physics properties
//----------------------------------------------------------------------
Gravity gravity(Vec3d(0.0, 0.0, 0.0));
Real left_pressure = 5.0;
Real right_pressure = 0.0;
Real delta_p = left_pressure - right_pressure;
Real mu_f = 3.5e-3;
Real U_f = delta_p * radius * radius / (8 * full_length * mu_f);
Real U_max = 2.5 * U_f; // analytical maximum velocity will be 2*U_f, using 2.5 for safty here
Real c_f = 10 * U_max;
Real rho0_f = 1050;
int simtk_resolution = 20;

class BidirectionalBufferCondition : public fluid_dynamics::BidirectionalBuffer
{
  private:
    double inlet_pressure_ = 0;

  public:
    BidirectionalBufferCondition(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction, double inlet_pressure)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction),
          inlet_pressure_(inlet_pressure) {}
    Real getTargetPressure(Real dt) override
    {
        return compute_pressure(inlet_pressure_);
    }
};
class InflowPressure : public fluid_dynamics::FlowPressureBuffer
{
  private:
    double inlet_pressure_ = 0;
    string name_;

  public:
    InflowPressure(BodyPartByCell &constrained_region, Vecd normal_vector, double inlet_pressure, string name)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector), inlet_pressure_(inlet_pressure), name_(name)
    {
    }

    Real getTargetPressure(Real dt) override
    {
        return compute_pressure(inlet_pressure_);
    }

    void setupDynamics(Real dt = 0.0) override {}
};
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-wall_thickness, -radius - wall_thickness, -radius - wall_thickness);
Vec3d domain_upper_bound(full_length + wall_thickness, radius + wall_thickness, radius + wall_thickness);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Geomtrey
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        std::cout << "WaterBlock" << std::endl;
        add<TriangleMeshShapeCylinder>(SimTK_inlet_normal, radius,
                                       full_length * 0.5, simtk_resolution, Vec3d(full_length * 0.5, 0, 0));
    }
};

class WallBlock : public ComplexShape
{
  public:
    explicit WallBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(SimTK_inlet_normal, radius + wall_thickness,
                                       full_length * 0.5 + wall_thickness, simtk_resolution, Vec3d(full_length * 0.5, 0, 0));
        subtract<TriangleMeshShapeCylinder>(SimTK_inlet_normal, radius,
                                            full_length * 0.5 + wall_thickness, simtk_resolution, Vec3d(full_length * 0.5, 0, 0));
    }
};

StdVec<Vec3d> generate_observation_location_axial()
{
    StdVec<Vec3d> pos_;
    size_t number_of_oberver = 50;
    double dx = full_length / number_of_oberver;
    Vec3d pos = inlet_center;
    while (pos[axial_direction] < full_length)
    {

        pos_.push_back(pos);
        pos[axial_direction] += dx;
    }
    return pos_;
};

StdVec<Vec3d> generate_observation_location_radial()
{
    StdVec<Vec3d> pos_;
    size_t number_of_oberver = 25;
    double x = full_length * 0.5;
    double dz = radius * 2 / number_of_oberver;
    Vec3d pos = Vec3d(x, 0, -radius);
    while (pos[2] <= radius)
    {
        pos_.push_back(pos);
        pos[2] += dz;
    }
    return pos_;
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
    std::cout << "U_f = " << U_f << std::endl;
    std::cout << "Reynolds number = " << diameter * 2 * U_f * rho0_f / mu_f << std::endl;
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setGenerateRegressionData(false);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);

    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBlock>("Wall"));
    wall_boundary.defineAdaptationRatios(1.3, 1.0);
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    (sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody velocity_observer_axial(sph_system, "VelocityObserverAxial");
    StdVec<Vec3d> observer_location_axial = generate_observation_location_axial();
    velocity_observer_axial.generateParticles<ParticleGeneratorObserver>(observer_location_axial);

    ObserverBody velocity_observer_radial(sph_system, "VelocityObserverRadial");
    StdVec<Vec3d> observer_location_radial = generate_observation_location_radial();
    velocity_observer_radial.generateParticles<ParticleGeneratorObserver>(observer_location_radial);
    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_axial_contact(velocity_observer_axial, {&water_block});
    ContactRelation velocity_observer_radial_contact(velocity_observer_radial, {&water_block});

    ObservedQuantityRecording<Vecd> write_fluid_velocity_axial("Velocity", velocity_observer_axial_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity_radial("Velocity", velocity_observer_radial_contact);
    ObservedQuantityRecording<Real> write_fluid_pressure_axial("Pressure", velocity_observer_axial_contact);
    ObservedQuantityRecording<Real> write_fluid_pressure_radial("Pressure", velocity_observer_radial_contact);

    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	check whether run particle relaxation for body fitted particle distribution.
    //----------------------------------------------------------------------

    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;

        /** Define inner relationship. */
        InnerRelation wall_inner(wall_boundary);
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_wall_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_wall_to_vtp(wall_boundary);
        /** Write the particle reload files. */

        ReloadParticleIO write_particle_reload_files(wall_boundary);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(wall_inner);
        /**
         * @brief 	Particle relaxation starts here.
         */
        random_wall_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();

        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the wall body N = " << ite_p << "\n";
                write_wall_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0.0);
        return 0;
    }
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** delete outflow particles */
    BodyAlignedBoxByCell left_disposer(
        water_block, makeShared<AlignedBoxShape>(left_disposer_transform, buffer_half_size));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer, xAxis);
    BodyAlignedBoxByCell right_disposer(
        water_block, makeShared<AlignedBoxShape>(right_disposer_transform, buffer_half_size));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        boundary_indicator(water_block_inner, water_block_contact);
    /** bidrectional buffer */
    BidirectionalBufferCondition left_emitter_inflow_injection(
        water_block, makeShared<AlignedBoxShape>(left_bc_transform, buffer_half_size), 10, xAxis, left_pressure);
    BidirectionalBufferCondition right_emitter_inflow_injection(
        water_block, makeShared<AlignedBoxShape>(right_bc_transform, buffer_half_size), 10, xAxis, right_pressure);
    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /** pressure boundary condition. */
    BodyRegionByCell left_pressure_region(water_block, makeShared<AlignedBoxShape>(left_bc_transform, buffer_half_size));
    SimpleDynamics<InflowPressure> left_pressure_condition(left_pressure_region, inlet_normal, left_pressure, "left pressure region");
    BodyRegionByCell right_pressure_region(water_block, makeShared<AlignedBoxShape>(right_bc_transform, buffer_half_size));
    SimpleDynamics<InflowPressure> right_pressure_condition(right_pressure_region, inlet_normal, right_pressure, "right pressure region");

    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** momentum equation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** mass equation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>>
        transport_velocity_correction(water_block_inner, water_block_contact);

    /** output parameters */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");
    water_block.addBodyStateForRecording<Vec3d>("KernelSummation");
    velocity_observer_axial.addBodyStateForRecording<double>("Pressure");
    velocity_observer_radial.addBodyStateForRecording<double>("Pressure");
    BodyStatesRecordingToVtp body_states_recording(sph_system.sph_bodies_);
    /**
     * @brief Output.
     */
    /** Output the body states. */

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
    Real end_time = 5.0;     /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;           /**< Default acoustic time step sizes. */
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

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                dt = SMIN(dt, Dt - relaxation_time);
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

                velocity_observer_axial_contact.updateConfiguration();
                velocity_observer_radial_contact.updateConfiguration();
                write_fluid_velocity_axial.writeToFile();
                write_fluid_velocity_radial.writeToFile();
                write_fluid_pressure_axial.writeToFile();
                write_fluid_pressure_radial.writeToFile();

                body_states_recording.writeToFile();
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();
            // water_block.updateCellLinkedListWithParticleSort(100);
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();
        }
        TickCount t2 = TickCount::now();

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
