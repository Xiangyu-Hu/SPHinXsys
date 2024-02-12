/**
 * @file 	Carotid_flow.cpp
 * @brief 	3D pulsatile Carotid flow example
 * @details This is the one of the basic test cases for pressure boundary condition and bidirectional buffer.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h" 
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
#include "bidirectional_buffer.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;

std::string full_path_to_stl_file_wall = "./input/wall.stl";
std::string full_path_to_stl_file_fluid = "./input/fluid.stl";
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-7.0E-3, -5.0E-3, -33.0E-3), Vecd(13.0E-3, 11.0E-3, 25.0E-3));
/**
 * @brief Material properties of the fluid.
 */
Real Inlet_pressure = 0.0;
Real Outlet_pressure = 2666.44;
Real rho0_f = 1060.0;                   
Real mu_f = 0.00355;
Real U_f = 1.0;
Real c_f = 10.0*U_f;

Real resolution_ref = 0.15E-3;
Real scaling = 0.001;
Vecd translation(0.0, 0.0, 0.0);
StdVec<Vecd> observer_location = {Vecd(3.5E-3, 4.1E-3,13.2E-3)};
// inlet pipe
Vecd buffer_halfsize_inlet = Vecd(3.0 * resolution_ref, 4.0E-3, 4.0E-3);
Vecd buffer_translation_inlet = Vecd(1.7284E-3, 6.1295E-3, -30.4405E-3);
Vecd normal_vector_inlet = Vecd(0.1039, -0.0458, 0.9935);
Vecd vector_inlet = Vecd(0, -0.9989, -0.0460);
// thick outlet pipe
Vecd buffer_halfsize_1 = Vecd(3.0 * resolution_ref, 6.0E-3, 6.0E-3);
Vecd buffer_translation_1 = Vecd(-3.7724E-3, 1.0088E-3, 20.9530E-3);
Vecd normal_vector_1 = Vecd(-0.3162, 0.0, 0.9487);
Vecd vector_1 = Vecd(0.0, -1.0, 0.0);
//thin outlet pipe
Vecd buffer_halfsize_2 = Vecd(3.0 * resolution_ref, 5.0E-3, 5.0E-3);
Vecd buffer_translation_2 = Vecd(8.1739E-3, -0.0297E-3, 18.17418E-3);
Vecd normal_vector_2 = Vecd(0.0, 0.0, 1.0);
Vecd vector_2 = Vecd(0, -1.0, 0.0);

/**
 * @brief 	water body definition.
 */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_fluid, translation, scaling);
        Vecd halfsize_water(1.0E-3, 1.0E-3, 1.0E-3);
        Transform translation_water(Vecd(5.6E-3, 4.7E-3, 6.5E-3));
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
        subtract<TriangleMeshShapeSTL>(full_path_to_stl_file_wall, translation, scaling);
    }
};
/**
 * @brief 	wall body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_stl_file_wall, translation, scaling);
    }
};
/**
 * @brief 	Bidirectional buffer definition.
 */
class BidirectionalBufferConditionInlet : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferConditionInlet(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                      size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
       Real pressure = Inlet_pressure;
       return pressure;
    }
};
class BidirectionalBufferConditionOutlet1 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferConditionOutlet1(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real pressure = run_time < 0.1 ? 0.5 * Outlet_pressure * (1.0 - cos(Pi * run_time / 0.1)) : Outlet_pressure;
        return pressure;
    }
};
class BidirectionalBufferConditionOutlet2 : public fluid_dynamics::BidirectionalBuffer
{
  public:
    BidirectionalBufferConditionOutlet2(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                                 size_t body_buffer_width, int axis_direction)
        : fluid_dynamics::BidirectionalBuffer(real_body, shape_ptr,
                                              body_buffer_width, axis_direction) {}
    Real getTargetPressure(Real dt) override
    {       
        Real run_time = GlobalStaticVariables::physical_time_;
        Real pressure = run_time < 0.1 ? 0.5 * Outlet_pressure * (1.0 - cos(Pi * run_time / 0.1)) : Outlet_pressure;
        return pressure;
    }
};
/**
 * @brief 	inflow velocity definition.
 */
struct InflowVelocity
{
    Real u_ref_, t_ref_, u_ave;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(0.1), t_ref_(0.1), u_ave(0.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;
        if (run_time < t_ref_)
        {
            u_ave = u_ref_;  
        }
        else if (((run_time - t_ref_) - int((run_time - t_ref_) / 0.5) * 0.5) >= 0.0 
            && ((run_time - t_ref_) - int((run_time - t_ref_) / 0.5) * 0.5) <= 0.218)
        {
            u_ave = 0.5 * sin(4.0 * Pi * ((run_time - t_ref_) + 0.0160236));
        }
        else
        {
            u_ave = u_ref_;        
        }

        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;
       
        return target_velocity;            
    }
};
/**
 * @brief 	outlet pressure boundary definition.
 */
class OutflowPressure_1 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_1(BodyPartByCell &constrained_region, Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector) {}
    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real pressure = run_time < 0.1 ? 0.5 * Outlet_pressure * (1.0 - cos(Pi * run_time / 0.1)) : Outlet_pressure;
        return pressure;
    }
    void setupDynamics(Real dt = 0.0) override {}
};
class OutflowPressure_2 : public fluid_dynamics::FlowPressureBuffer
{
  public:
    OutflowPressure_2(BodyPartByCell &constrained_region,Vecd normal_vector)
        : fluid_dynamics::FlowPressureBuffer(constrained_region, normal_vector) {}
    Real getTargetPressure(Real dt) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real pressure = run_time < 0.1 ? 0.5 * Outlet_pressure * (1.0 - cos(Pi * run_time / 0.1)) : Outlet_pressure;
        return pressure;
    }
    void setupDynamics(Real dt = 0.0) override {}
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
    sph_system.setGenerateRegressionData(true);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape()->correctLevelSetSign();
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<ParticleGeneratorReload>(water_block.getName())
        : water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody velocity_observer(sph_system, "VelocityObserver");
    velocity_observer.generateParticles<ParticleGeneratorObserver>(observer_location);
    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation velocity_observer_contact(velocity_observer, {&water_block});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation water_block_inner(water_block);
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> water_block_random_inserted_body_particles(water_block);
        SimpleDynamics<RandomizeParticlePosition> wall_boundary_random_inserted_body_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp water_block_write_inserted_body_to_vtp({&water_block});
        BodyStatesRecordingToVtp wall_boundary_write_inserted_body_to_vtp({&wall_boundary});
        /** Write the particle reload files. */
        ReloadParticleIO water_block_write_particle_reload_files({&water_block});
        ReloadParticleIO wall_boundary_write_particle_reload_files({&wall_boundary});
        /** A  Physics relaxation step. */
        RelaxationStepInner water_block_relaxation_step_inner(water_block_inner);
        RelaxationStepInner wall_boundary_relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        water_block_random_inserted_body_particles.exec(0.25);
        wall_boundary_random_inserted_body_particles.exec(0.25);
        water_block_relaxation_step_inner.SurfaceBounding().exec();
        wall_boundary_relaxation_step_inner.SurfaceBounding().exec();
        water_block_write_inserted_body_to_vtp.writeToFile(0);
        wall_boundary_write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 3000)
        {
            water_block_relaxation_step_inner.exec();
            wall_boundary_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                water_block_write_inserted_body_to_vtp.writeToFile(ite_p);
                wall_boundary_write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        water_block_write_particle_reload_files.writeToFile(0);
        wall_boundary_write_particle_reload_files.writeToFile(0);
        return 0;
    }
    /**
     * @brief 	Define all numerical methods which are used in this case.
     */
    /**
     * @brief 	Methods used for time stepping.
     */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** delete outflow particles */
    BodyAlignedBoxByCell disposer_inlet(water_block, makeShared<AlignedBoxShape>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_inlet)),
        -vector_inlet),Vecd(buffer_translation_inlet)),buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_inlet(disposer_inlet, xAxis);
    BodyAlignedBoxByCell disposer_1(water_block, makeShared<AlignedBoxShape>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_1)), 
        vector_1),Vecd(buffer_translation_1)),buffer_halfsize_1));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_1(disposer_1, xAxis);
    BodyAlignedBoxByCell disposer_2(water_block, makeShared<AlignedBoxShape>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_2)),
        vector_2),Vecd(buffer_translation_2)),buffer_halfsize_2));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion_2(disposer_2, xAxis);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        free_stream_surface_indicator(water_block_inner, water_block_contact);
    /** bidrectional buffer */
    BidirectionalBufferConditionInlet emitter_inflow_injection_inlet(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)), 
        vector_inlet),(buffer_translation_inlet)), buffer_halfsize_inlet), 100, xAxis);
    BidirectionalBufferConditionOutlet1 emitter_inflow_injection1(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_1)), 
        -vector_1),(buffer_translation_1)), buffer_halfsize_1), 100, xAxis);
    BidirectionalBufferConditionOutlet2 emitter_inflow_injection2(
        water_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(-normal_vector_2)), 
        -vector_2),(buffer_translation_2)), buffer_halfsize_2), 100, xAxis);
    /** output parameters */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("BufferParticleIndicator");  
    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /* Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** pressure boundary condition. */ 
    BodyAlignedBoxByCell inflow_pressure_region(water_block, makeShared<AlignedBoxShape>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_inlet)),vector_inlet),
        (buffer_translation_inlet)),buffer_halfsize_inlet));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(inflow_pressure_region);
    BodyRegionByCell outflow_pressure_region1(water_block, makeShared<TransformShape<GeometricShapeBox>>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_1)), vector_1),
        (buffer_translation_1)),buffer_halfsize_1));
    SimpleDynamics<OutflowPressure_1> outflow_pressure_condition1(outflow_pressure_region1, normal_vector_1);
    BodyRegionByCell outflow_pressure_region2(water_block, makeShared<TransformShape<GeometricShapeBox>>
        (Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(normal_vector_2)), vector_2),
        (buffer_translation_2)),buffer_halfsize_2));
    SimpleDynamics<OutflowPressure_2> outflow_pressure_condition2(outflow_pressure_region2, normal_vector_2);
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
    free_stream_surface_indicator.exec();
    emitter_inflow_injection_inlet.tag_buffer_particles.exec();
    emitter_inflow_injection1.tag_buffer_particles.exec();
    emitter_inflow_injection2.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();

    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = 0.0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 0.5;   /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;

    /** Output the start states of bodies. */
    body_states_recording.writeToFile(0);
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
                emitter_buffer_inflow_condition.exec();
                outflow_pressure_condition1.exec(dt); 
                outflow_pressure_condition2.exec(dt);           
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
            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            /** Water block configuration and periodic condition. */
            emitter_inflow_injection_inlet.injection.exec();
            emitter_inflow_injection1.injection.exec();
            emitter_inflow_injection2.injection.exec();
            disposer_outflow_deletion_inlet.exec();
            disposer_outflow_deletion_1.exec();
            disposer_outflow_deletion_2.exec();
            water_block.updateCellLinkedList();
            water_block_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            free_stream_surface_indicator.exec();
            emitter_inflow_injection_inlet.tag_buffer_particles.exec();
            emitter_inflow_injection1.tag_buffer_particles.exec();
            emitter_inflow_injection2.tag_buffer_particles.exec();
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
