/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
// #include "fsi2.h" // case file to setup the test case
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
// Real rho0_s = 10.0; /**< Reference density.*/
// Real poisson = 0.4; /**< Poisson ratio.*/
// Real Ae = 1.4e5;    /**< Normalized Youngs Modulus. */
// Real Youngs_modulus = Ae;
Real DL = 8.0;                                             /**< Channel length. */
Real DH = 2.0;                                             /**< Channel height. */
Real resolution_ref = DH / 50.0;  //更改分辨率也会影响计算出的速度最大值，/20的时候最最大值0.68左右，/100时最大值0.73左右   /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 3.0;
Real Outlet_pressure = 0.0;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
//Real mu_f = 1.0e-2;
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
//Real U_f = 1.0;
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL - 2.5 * resolution_ref, 0.5 * DH);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        Real pressure = Inlet_pressure;
        return pressure;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

class ConstantVelocityConstraint : public MotionConstraint<SPHBody>
{
public:
    /** 构造函数：传入一个粒子集合（BodyPartByParticle），例如整个 water_block */
    ConstantVelocityConstraint(SPHBody &body_part)
        : MotionConstraint<SPHBody>(body_part),
        physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")){};
    virtual ~ConstantVelocityConstraint() {};

    /** 每个 time step 自动对 body_part 中的粒子调用一次这里 */
    void update(size_t index_i, Real dt)
    {
        Real t = *physical_time_;
        vel_[index_i][0] = 0.04564;  // x 方向速度
        vel_[index_i][1] = 0.0;  // y 方向速度
    }
    protected:
    Real *physical_time_;
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
// struct InflowVelocity
// {
//     Real u_ave;

//     template <class BoundaryConditionType>
//     InflowVelocity(BoundaryConditionType &boundary_condition)
//         : u_ave(0.0) {}

//     Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
//     {
//         Vecd target_velocity = Vecd::Zero();

        
//         target_velocity[0] = 1.0;
//         target_velocity[1] = 0.0;

//         return target_velocity;
//     }
// };

//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
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

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------

class Vessel : public MultiPolygonShape
{
  public:
    explicit Vessel(const std::string &shape_name) : MultiPolygonShape(shape_name)
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

// //----------------------------------------------------------------------
// //	Constrain Shape definition.
// //----------------------------------------------------------------------

// std::vector<Vecd> createLeftwallshape()
// {
//     std::vector<Vecd> left_wall_shape;
//     left_wall_shape.push_back(Vecd(-0.0001, -0.0003));
//     left_wall_shape.push_back(Vecd(-0.0001, 0.002));
//     left_wall_shape.push_back(Vecd(0.0001, 0.002));
//     left_wall_shape.push_back(Vecd(0.0001, -0.0003));
//     left_wall_shape.push_back(Vecd(-0.0001, -0.0003));
//     return left_wall_shape;
// }
// std::vector<Vecd> createRightwallshape()
// {
//     std::vector<Vecd> right_wall_shape;
//     right_wall_shape.push_back(Vecd(0.00375, -0.003));
//     right_wall_shape.push_back(Vecd(0.00375, 0.002));
//     right_wall_shape.push_back(Vecd(0.005, 0.002));
//     right_wall_shape.push_back(Vecd(0.005, -0.003));
//     right_wall_shape.push_back(Vecd(0.00375, -0.003));
//     return right_wall_shape;
// }   

// /** create constrain shape. */
// MultiPolygon createBeamBaseShape()
// {
//     MultiPolygon multi_polygon;
//     multi_polygon.addAPolygon(createLeftwallshape(), ShapeBooleanOps::add);
//     multi_polygon.addAPolygon(createRightwallshape(), ShapeBooleanOps::add);
//     return multi_polygon;
// }

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the entrance of the channel. */
    size_t number_observation_points = 80;
    Real range_of_measure = DH;         // 全高度
    Real start_of_measure = 0.0;        // 从下壁面开始
    /** the measuring locations */
    for (size_t i = 0; i < number_observation_points; ++i)
    {
        Vec2d point_coordinate(
            resolution_ref * 75.0, 
            range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure
        );
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    // sph_system.setReloadParticles(true);         // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();; // handle command line arguments
    //IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(1.0);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);
    
    // SolidBody vessel_boundary(sph_system, makeShared<Vessel>("Vessel"));
    // // vessel_boundary.defineAdaptationRatios(1.15, 2.0);
    // vessel_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    // vessel_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    // (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    //     ? vessel_boundary.generateParticles<BaseParticles, Reload>(vessel_boundary.getName())
    //     : vessel_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody vessel_boundary(sph_system, makeShared<Vessel>("Vessel"));
    vessel_boundary.defineMaterial<Solid>();
    vessel_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.defineAdaptationRatios(0.25, 1.0);
    fluid_observer.generateParticles<ObserverParticles>(createObservationPoints());

    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    // if (sph_system.RunParticleRelaxation())
    // {
    //     //----------------------------------------------------------------------
    //     //	Define body relation map used for particle relaxation.
    //     //----------------------------------------------------------------------
    //     InnerRelation insert_body_inner(vessel_boundary);
    //     //----------------------------------------------------------------------
    //     //	Methods used for particle relaxation.
    //     //----------------------------------------------------------------------
    //     using namespace relax_dynamics;
    //     SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(vessel_boundary);
    //     //RelaxationStepInner relaxation_step_inner(insert_body_inner);
    //     BodyStatesRecordingToVtp write_insert_body_to_vtp(vessel_boundary);    
    //     BodyStatesRecordingToVtp write_real_body_states(sph_system);
    //     ReloadParticleIO write_particle_reload_files(vessel_boundary);
    //     //----------------------------------------------------------------------
    //     //	Particle relaxation starts here.
    //     //----------------------------------------------------------------------
    //     random_insert_body_particles.exec(0.25);
    //     //relaxation_step_inner.SurfaceBounding().exec();
    //     write_insert_body_to_vtp.writeToFile(0);
    //     //----------------------------------------------------------------------
    //     //	Relax particles of the insert body.
    //     //----------------------------------------------------------------------
    //     int ite_p = 0;
    //     while (ite_p < 1000)
    //     {
    //         //relaxation_step_inner.exec();
    //         ite_p += 1;
    //         if (ite_p % 200 == 0)
    //         {
    //             std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    //             write_insert_body_to_vtp.writeToFile(ite_p);
    //         }
    //     }
    //     std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    //     /** Output results. */
    //     write_particle_reload_files.writeToFile(0);
    //     write_real_body_states.writeToFile(0);
    //     return 0;
    // }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    //InnerRelation vessel_boundary_inner(vessel_boundary);
    ContactRelation water_block_contact(water_block, {&vessel_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //ContactRelation insert_body_contact(vessel_boundary, {&water_block});
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
    SimpleDynamics<NormalDirectionFromBodyShape> vessel_boundary_normal_direction(vessel_boundary);
    SimpleDynamics<ConstantVelocityConstraint> update_water_velocity(water_block);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    // InteractionWithUpdate<LinearGradientCorrectionMatrixInner> vessel_boundary_corrected_configuration(vessel_boundary_inner);

    // Dynamics1Level<solid_dynamics::Integration1stHalfPK2> vessel_boundary_stress_relaxation_first_half(vessel_boundary_inner);
    // Dynamics1Level<solid_dynamics::Integration2ndHalf> vessel_boundary_stress_relaxation_second_half(vessel_boundary_inner);

    //ReduceDynamics<solid_dynamics::AcousticTimeStep> vessel_boundary_computing_time_step_size(vessel_boundary);
    // BodyRegionByParticle beam_base(vessel_boundary, makeShared<MultiPolygonShape>(createBeamBaseShape()));
    // SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);

    //SimpleDynamics<VonMisesStress> calculate_stress(vessel_boundary);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

    AlignedBoxPartByCell left_emitter(water_block, AlignedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    AlignedBoxPartByCell right_emitter(water_block, AlignedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    // InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_density_by_summation(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    //SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    // PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    // PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    // //----------------------------------------------------------------------
    // //	Algorithms of FSI.
    // //----------------------------------------------------------------------
    // solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(vessel_boundary);
    // SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(vessel_boundary);
    // InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    // InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations and regression tests of the simulation.
    //----------------------------------------------------------------------
     BodyStatesRecordingToVtp body_states_recording(sph_system);
     body_states_recording.addToWrite<Real>(water_block, "Pressure");
     body_states_recording.addToWrite<int>(water_block, "Indicator");
     body_states_recording.addToWrite<Real>(water_block, "Density");
     body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
     ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    //body_states_recording.addToWrite<Vecd>(vessel_boundary, "Velocity");
    //body_states_recording.addToWrite<Matd>(vessel_boundary, "DeformationGradient");
    //body_states_recording.addToWrite<Real>(vessel_boundary, "VonMisesStress");
    // RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_centerline_velocity("Velocity", velocity_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    vessel_boundary_normal_direction.exec();
    update_water_velocity.exec();
    //vessel_boundary_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 1000.0;
    Real output_interval =0.1;
    Real dt = 0.0;
    int observation_sample_interval = screen_output_interval * 2;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    // damping ratio for vessel
    // Real physical_viscosity_wall = 100.0;
    // DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
    // wall_velocity_damping(0.6, vessel_boundary_inner, "Velocity", physical_viscosity_wall);
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            // Real Dt = get_fluid_advection_time_step_size.exec();
            // // update_density_by_summation.exec();
            // // viscous_force.exec();
            // // transport_correction.exec();
            // /** Update normal direction on elastic body.*/
            // //insert_body_update_normal.exec();
            // size_t inner_ite_dt = 0;
            // size_t inner_ite_dt_s = 0;
            // Real relaxation_time = 0.0;
            // while (relaxation_time < Dt)
            // {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            //Real Dt = 0.001;
            update_density_by_summation.exec();
            //kernel_correction_complex.exec();
            viscous_force.exec();
            transport_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                //dt = 0.001;
                pressure_relaxation.exec(dt);
               //update_water_velocity.exec();
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                //inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
            time_instance = TickCount::now();
            
            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();
            
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            //boundary_indicator.exec();
            water_block.updateCellLinkedList();
            //boundary_indicator.exec();
            water_block_complex.updateConfiguration();
            //boundary_indicator.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
            //vessel_boundary.updateCellLinkedList();
            // //insert_body_contact.updateConfiguration();
            // // interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            //boundary_indicator.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
            //boundary_indicator.exec();
            // //water_block_complex.updateConfiguration();
            // //interval_updating_configuration += TickCount::now() - time_instance;
        }

        TickCount t2 = TickCount::now();
        //calculate_stress.exec();
        body_states_recording.writeToFile();
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
      std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";


    return 0;
}