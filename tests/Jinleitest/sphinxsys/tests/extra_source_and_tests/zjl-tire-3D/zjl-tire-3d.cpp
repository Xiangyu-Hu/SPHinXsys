/**
 * @file 	poiseuille_flow_shell.cpp
 * @brief 	3D poiseuille flow interaction with shell example
 * @details This is the one of the basic test cases for validating fluid-rigid shell interaction
 * @author  Weiyi Kong
 */

#include "sphinxsys.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
#include "base_general_dynamics.h"
#include "base_data_package.h"

using namespace SPH;
/**
 * @brief 3D Tire-water interaction simulation
 * @author Jinlei Zhou
 * @date 2025-06
 */

//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string tire_path = "/Users/chloe-jinlei/Desktop/Jinlei_Tire/sphinxsys/tests/extra_source_and_tests/zjl-tire-3D/data/tirereshape.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Vec3d translation_tire(0.0, 0.0, 0.0);                                   /**< Initial translation of the tire, unit: meters (m) */
Real length_scale_tire = pow(10, -3);                                    /**< Length scale factor, dimensionless (unitless) */
Real length_scale = pow(10, -3);                                          /**< Length scale factor, dimensionless (unitless) */
Real resolution_ref = 0.2 * length_scale;                                 /**< Initial reference particle spacing, unit: meters (m) */
Vec3d domain_lower_bound = Vec3d(-50.0, -50.0, -50.0) * length_scale;      /**< Lower boundary of the system domain, unit: meters (m) */
Vec3d domain_upper_bound = Vec3d(50.0, 50.0, 50.0) * length_scale;        /**< Upper boundary of the system domain, unit: meters (m) */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound); /**< Defines the bounding box of the system domain */
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class Tire : public ComplexShape
{
  public:
    explicit Tire(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(tire_path, translation_tire, length_scale_tire);
    }
};
int main(int ac, char *av[]) {
    //--------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    //-------------------------------------------------------------------
    SolidBody road(sph_system, makeShared<Tire>("Tire"));
    road.defineMaterial<Solid>();
    road.generateParticles<BaseParticles, Lattice>();
    
    BodyStatesRecordingToVtp write_stent_to_vtp(road);
    write_stent_to_vtp.writeToFile(0);

  return 0;
}























// //----------------------------------------------------------------------
// // Geometry and resolution (2D)
// //----------------------------------------------------------------------
// Real resolution_ref = 0.0025;
// // Real R_outer = 0.300;
// // Real R_inner = 0.200;
// Real rim_thickness = 0.1;
// // Real road_length = 0.400;
// // Real road_height = 0.05;
// //Real gravity1_ = 30.81;
// Real gravity2_ = 30.81;
// Vec2d rim_upper_point(road_length*0.5,road_height+R_outer+0.5*rim_thickness);
// //Real rho0_s = 1200.0;//轮胎密度 整体给还是轮毂和橡胶分开给？
// Real Youngs_modulus = 1e4;//20-80MPa
// //Real poisson = 0.45;
// Real physical_viscosity = 1e6;
// Real U_f = 1.0;                                               /**< Characteristic velocity. */
// Real rho0_f = 1000.0;
// Real c_f = 10.0;
// Real mu_f = 0.001;
// Real water_length = 1.0;
// Real water_height = 0.01;

// Real gate_moving_time = 0.30;
// Real contact_time = 0.25;
// Real angular_0 = 10.0; // rad/s

// BoundingBox system_domain_bounds(Vec2d(-0.8, -0.5), Vec2d(1.3, 0.7));

// //----------------------------------------------------------------------
// // Tire geometry (2D circular ring)
// //----------------------------------------------------------------------
// class Tire : public MultiPolygonShape {
// public:
//     explicit Tire(const std::string &shape_name) : MultiPolygonShape(shape_name) {
//         multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_outer, 50, ShapeBooleanOps::add);
//         multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
//     }
// };

// //----------------------------------------------------------------------
// // road (solid domain)
// //----------------------------------------------------------------------
// class Road : public MultiPolygonShape {
// public:
//     explicit Road(const std::string &shape_name) : MultiPolygonShape(shape_name) {
//         std::vector<Vecd> road_shape = {
//             Vecd(0.0, 0.0), Vecd(0.0, road_height),
//             Vecd(road_length+1.0, road_height), Vecd(road_length+1.0, 0.0), Vecd(0.0, 0.0)
//         };
//         multi_polygon_.addAPolygon(road_shape, ShapeBooleanOps::add);
//     }
// };

// //----------------------------------------------------------------------
// // Water block (fluid domain)
// //----------------------------------------------------------------------
// class WaterBlock : public MultiPolygonShape {
// public:
//     explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name) {
//         std::vector<Vecd> water_shape = {
//             Vecd(0.5, road_height), Vecd(0.5, road_height+water_height),
//             Vecd(water_length+0.2, road_height+water_height), Vecd(water_length+0.2, road_height), Vecd(0.5,road_height)
//         };
//         multi_polygon_.addAPolygon(water_shape, ShapeBooleanOps::add);
//     }
// };

// //----------------------------------------------------------------------
// // Wheel rim geometry (2D circular ring)
// //----------------------------------------------------------------------
// /** create the Wheel rim as constrain shape. */
// MultiPolygon createRimBaseShape()
// {
//     MultiPolygon multi_polygon;
//     multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner, 50, ShapeBooleanOps::add);
//     multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
//     return multi_polygon;
// }
// //----------------------------------------------------------------------
// // define Tire Motion
// //----------------------------------------------------------------------
// /**
//  * @class	MotionConstraint
//  * @brief	Base class for constraining with prescribed motion.
//  * 			Exact motion function will be defined in derive class.
//  * 			Note that, we do not impose acceleration, so that this constraint
//  * 			must be imposed after updating particle velocity by forces
//  * 			and before updating particle position.
//  */
// //在tire-presim-2d里得到经过0.1s的车身自重施加在轮毂上后整个轮胎圆心下降了0.003-0.004，假设一开始圆心下降的时间是0.05s，那么圆心匀速下降的速度为0.0035/0.05=0.07m/s
//  class TireMotionConstraint : public MotionConstraint<BodyPartByParticle>
// {
//   public:
//     TireMotionConstraint(BodyPartByParticle &body)
//         : MotionConstraint<BodyPartByParticle>(body),
//           physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")){};
//     virtual ~TireMotionConstraint(){};
//     void update(size_t index_i, Real dt)
//     {   
//        // 当前时间
//     Real t = *physical_time_;

//     // 平动速度
//     Vecd v_trans;
//     if (t < 0.05)
//     {
//         v_trans = Vecd(0.0, -0.3);
//     }
//     else
//     {
//         v_trans = Vecd(1.0, 0.0);
//     }

//     // 圆心位置
//     Real cx = 0.5 * road_length;
//     Real cy = road_height + R_outer;
//     if (t < 0.05)
//     {
//         cy += -0.3 * t;
//     }
//     else
//     {
//         cy += -0.3 * 0.05;
//         cx += 1.0 * (t - 0.05);
//     }

//     // 旋转分量
//     Real x = pos_[index_i][0] - cx;
//     Real y = pos_[index_i][1] - cy;
//     Vecd v_rot(angular_0 * y,
//                -angular_0 * x);

//     // 合成
//     vel_[index_i] = v_rot + v_trans;
//     };

//   protected:
//     Real *physical_time_;
// };

// //----------------------------------------------------------------------
// // Main program
// //----------------------------------------------------------------------
// int main(int ac, char *av[]) {
//     //--------------------------------------------------------------------
//     SPHSystem sph_system(system_domain_bounds, resolution_ref);
//     sph_system.setRunParticleRelaxation(false);
//     sph_system.setReloadParticles(true);
//     sph_system.handleCommandlineOptions(ac, av);
//     IOEnvironment io_environment(sph_system);

//     //-------------------------------------------------------------------

//     SolidBody tire(sph_system, makeShared<Tire>("Tire"));
//     tire.defineBodyLevelSetShape()->writeLevelSet(sph_system);
//     //tire.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
//     tire.defineMaterial<TireBodyComposite>();
//         (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
//         ? tire.generateParticles<BaseParticles, Reload>(tire.getName())
//         : tire.generateParticles<BaseParticles, Lattice>();
    
//     FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
//     water_block.defineAdaptationRatios(1.15,1.0);
//     water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
//     water_block.generateParticles<BaseParticles, Lattice>();

//     SolidBody road(sph_system, makeShared<Road>("Road"));
//     road.defineMaterial<Solid>();
//     road.generateParticles<BaseParticles, Lattice>();

//     ObserverBody tire_observer(sph_system, "BallObserver");
//     tire_observer.generateParticles<ObserverParticles>(StdVec<Vec2d>{rim_upper_point});

//     //----------------------------------------------------------------------
//     //	Run particle relaxation for body-fitted distribution if chosen.
//     //----------------------------------------------------------------------
//     if (sph_system.RunParticleRelaxation())
//     {
//         //----------------------------------------------------------------------
//         //	Define body relation map used for particle relaxation.
//         //----------------------------------------------------------------------
//         InnerRelation insert_body_inner(tire);
//         //----------------------------------------------------------------------
//         //	Methods used for particle relaxation.
//         //----------------------------------------------------------------------
//         using namespace relax_dynamics;
//         SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(tire);
//         RelaxationStepInner relaxation_step_inner(insert_body_inner);
//         BodyStatesRecordingToVtp write_insert_body_to_vtp(tire);    
//         BodyStatesRecordingToVtp write_real_body_states(sph_system);
//         ReloadParticleIO write_particle_reload_files(tire);
//         //----------------------------------------------------------------------
//         //	Particle relaxation starts here.
//         //----------------------------------------------------------------------
//         random_insert_body_particles.exec(0.25);
//         relaxation_step_inner.SurfaceBounding().exec();
//         write_insert_body_to_vtp.writeToFile(0);
//         //----------------------------------------------------------------------
//         //	Relax particles of the insert body.
//         //----------------------------------------------------------------------
//         int ite_p = 0;
//         while (ite_p < 1000)
//         {
//             relaxation_step_inner.exec();
//             ite_p += 1;
//             if (ite_p % 200 == 0)
//             {
//                 std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
//                 write_insert_body_to_vtp.writeToFile(ite_p);
//             }
//         }
//         std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
//         /** Output results. */
//         write_particle_reload_files.writeToFile(0);
//         write_real_body_states.writeToFile(0);
//         return 0;
//     }
//     //----------------------------------------------------------------------
//     //	Define body relation map.
//     //	The contact map gives the topological connections between the bodies.
//     //	Basically the the range of bodies to build neighbor particle lists.
//     //  Generally, we first define all the inner relations, then the contact relations.
//     //----------------------------------------------------------------------
//     InnerRelation insert_body_inner(tire);
//     InnerRelation water_block_inner(water_block);
//     ContactRelation water_block_contact(water_block, {&road, &tire});
//     ContactRelation insert_body_contact(tire, {&water_block});
//     // ContactRelation tire_road_contact(tire, {&road});  
//     // ContactRelation road_tire_contact(road, {&tire});应该是只能用于observe body和SPH body之间
//     SurfaceContactRelation tire_road_contact(tire, {&road});
//     //SurfaceContactRelation road_tire_contact(road, {&tire});需要反过来在定义一次吗
//     //ContactRelation tire_observer_contact(tire_observer, {&tire});
//     ComplexRelation water_block_complex(water_block_inner, water_block_contact);
//     //----------------------------------------------------------------------
//     //	Define the main numerical methods used in the simulation.
//     //	Note that there may be data dependence on the constructors of these methods.
//     //----------------------------------------------------------------------
//     //SimpleDynamics<TireInitialCondition> tire_initial_rotation_velocity(tire);
//     //车身重力仅作用于轮毂算出来的加速度  a=F/m=5000N/(Area of wheel rim*1*density of wheel) = 5000/（0.377*1200）=11.052
//     BodyRegionByParticle wheel_rim(tire, makeShared<MultiPolygonShape>(createRimBaseShape()));
//     //Gravity gravity1(Vec2d(0.0,-gravity1_));
//     //SimpleDynamics<CarweightForce<Gravity>> constant_gravity(wheel_rim, gravity1);//改成部分，轮毂加速度
//     Gravity gravity2(Vec2d(0.0,-gravity2_));
//     SimpleDynamics<GravityForce<Gravity>> constant_gravity_water(water_block, gravity2);
//     SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(road);
//     SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(tire);
//     SimpleDynamics<TireMotionConstraint> update_tire_velocity(wheel_rim);
    
//     //SimpleDynamics<TireSpinningConstraint> update_tire_angular_velocity(tire);
//     InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
//     InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
//     InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);
//     SimpleDynamics<TireMaterialInitialization> composite_material_id(tire);
//     /** stress relaxation for the tire. */
//     Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
//     Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
//     ReduceDynamics<solid_dynamics::AcousticTimeStep> insert_body_computing_time_step_size(tire);
//     //BodyRegionByParticle beam_base(tire, makeShared<MultiPolygonShape>(createRimBaseShape()));
//     SimpleDynamics<VonMisesStress> calculate_stress(tire);
//     //----------------------------------------------------------------------
//     //	Algorithms of fluid dynamics.
//     //----------------------------------------------------------------------
//     Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
//     Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
//     InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
//     InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_block_contact);
//     InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

//     ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
//     ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
//     // InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
//     //----------------------------------------------------------------------
//     //	Algorithms of FSI.
//     //----------------------------------------------------------------------
//     solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(tire);
//     SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(tire);
//     InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
//     InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);
//     //----------------------------------------------------------------------
//     //	Define the configuration related particles dynamics.
//     //----------------------------------------------------------------------
//     ParticleSorting particle_sorting(water_block);
//     /** Algorithms for tire-road contact. */
//     InteractionDynamics<solid_dynamics::ContactFactorSummation> tire_update_contact_density(tire_road_contact);
//     InteractionWithUpdate<solid_dynamics::ContactForceFromWall> tire_compute_solid_contact_forces(tire_road_contact);
//     // damping ratio for tire
//     Real physical_viscosity_wall = 1000.0;
//     DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
//     wall_velocity_damping(0.4, insert_body_inner, "Velocity", physical_viscosity_wall);
//     //ReduceDynamics<solid_dynamics::AcousticTimeStep> tire_get_time_step_size(tire, 0.45);
    
//     //BodyRegionByParticle holder(wheel_rim, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_holder), halfsize_holder));
//     //SimpleDynamics<FixedInAxisDirection> constrain_holder(holder, Vecd(1.0, 0.0));
//     //----------------------------------------------------------------------
//     //	Define the methods for I/O operations and observations of the simulation
//     //----------------------------------------------------------------------
//     BodyStatesRecordingToVtp write_real_body_states(sph_system);
//     // BodyStatesRecordingToVtp write_ball_state(tire);
//     //ObservedQuantityRecording<Vecd>write_tire_displacement("Position", tire_observer_contact);
//     write_real_body_states.addToWrite<Vecd>(tire, "Position");
//     write_real_body_states.addToWrite<Real>(water_block, "Pressure");
//     write_real_body_states.addToWrite<Real>(tire, "VonMisesStress");
//     //----------------------------------------------------------------------
//     //	Prepare the simulation with cell linked list, configuration
//     //	and case specified initial condition if necessary.
//     //----------------------------------------------------------------------
//     /** initialize cell linked lists for all bodies. */
//     sph_system.initializeSystemCellLinkedLists();
//     sph_system.initializeSystemConfigurations();
//     //tire_initial_rotation_velocity.exec();
//     wall_boundary_normal_direction.exec();
//     insert_body_normal_direction.exec();
//     insert_body_corrected_configuration.exec();
//     constant_gravity_water.exec();
//     /** initialize material ids for the tire. */
//     composite_material_id.exec();
//     //----------------------------------------------------------------------
//     //	Setup for time-stepping control
//     //----------------------------------------------------------------------
//     Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
//     int ite = 0;
//     Real T0 = 1.0;
//     Real end_time = T0;
//     Real output_interval = T0/200;
//     Real Dt = 0.01 * output_interval;
//     Real dt = 0.0;
//     size_t number_of_iterations = 0;
//     int screen_output_interval = 1;
//     //----------------------------------------------------------------------
//     //	Statistics for CPU time
//     //----------------------------------------------------------------------
//     TickCount t1 = TickCount::now();
//     TimeInterval interval;
//     TimeInterval interval_computing_time_step;
//     TimeInterval interval_computing_pressure_relaxation;
//     TimeInterval interval_updating_configuration;
//     TickCount time_instance;
//     //----------------------------------------------------------------------
//     //	First output before the main loop.
//     //----------------------------------------------------------------------
//     write_real_body_states.writeToFile();
//     //----------------------------------------------------------------------
//     //	Main loop starts here.
//     //----------------------------------------------------------------------
//      while (physical_time < end_time)
//     {
//         Real integration_time = 0.0;
//         while (integration_time < output_interval)
//         {   
//             Real Dt = get_fluid_advection_time_step_size.exec();
//             update_density_by_summation.exec();
//             viscous_force.exec();
//             transport_correction.exec();

//             /** FSI for viscous force - only after 0.1s */
//             if (physical_time > 0.001)
//             {
//             viscous_force_from_fluid.exec();
//             }

//             /** Update normal direction on elastic body.*/
//             insert_body_update_normal.exec();
//             size_t inner_ite_dt = 0;
//             size_t inner_ite_dt_s = 0;
//             Real relaxation_time = 0.0;
//             while (relaxation_time < Dt)
//             {   
//                 Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
//                 /** Fluid pressure relaxation */
//                 pressure_relaxation.exec(dt);
                
//                 /** FSI for pressure force - only after 0.1s */
//                 if (physical_time > 0.001)
//                 {
//                 pressure_force_from_fluid.exec();
//                 }

//                 /** Fluid density relaxation */
//                 kernel_summation.exec();
//                 density_relaxation.exec(dt);

//                 //  /** Solid dynamics. */
//                 inner_ite_dt_s = 0;
//                 Real dt_s_sum = 0.0;
//                 average_velocity_and_acceleration.initialize_displacement_.exec();
//                 while (dt_s_sum < dt)
//                 {
//                     Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt - dt_s_sum);
//                     tire_update_contact_density.exec();
//                     tire_compute_solid_contact_forces.exec();
//                     insert_body_stress_relaxation_first_half.exec(dt_s);
//                     update_tire_velocity.exec();
//                     wall_velocity_damping.exec(dt_s);
//                     update_tire_velocity.exec();
//                     insert_body_stress_relaxation_second_half.exec(dt_s);
//                     tire.updateCellLinkedList();
//                     tire_road_contact.updateConfiguration();
//                     dt_s_sum += dt_s;
//                     inner_ite_dt_s++;
//                 }
                
//                 /** Update averages for FSI - only after 0.1s */
//                 if (physical_time > 0.001)
//                 {
//                 average_velocity_and_acceleration.update_averages_.exec(dt);
//                 }

//                 relaxation_time += dt;
//                 integration_time += dt;
//                 physical_time += dt;
//                 //boundary_indicator.exec();
//                 inner_ite_dt++;
//             }

//             if (number_of_iterations % screen_output_interval == 0)
//             {
//                 std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
//                           << physical_time
//                           << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
//                 //write_beam_tip_displacement.writeToFile(number_of_iterations);
//             }
//             number_of_iterations++;
//             time_instance = TickCount::now();

//             /** Water block configuration and periodic condition. */
//             //boundary_indicator.exec();
//             //periodic_condition.bounding_.exec();
//             if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
//             {
//                 particle_sorting.exec();
//             }
//             water_block.updateCellLinkedList();
//             //tire.updateCellLinkedList();
//             //tire_road_contact.updateConfiguration();
//             //periodic_condition.update_cell_linked_list_.exec();
//             //boundary_indicator.exec();
//             water_block_complex.updateConfiguration();
//             interval_updating_configuration += TickCount::now() - time_instance;
//             //tire.updateCellLinkedList();
//             insert_body_contact.updateConfiguration();
//             boundary_indicator.exec();
//         }
//         /** write run-time observation into file */
//         TickCount t2 = TickCount::now();
//         calculate_stress.exec();
//         //compute_vorticity.exec();
//         //write_ball_state.writeToFile();
//         write_real_body_states.writeToFile();
//         TickCount t3 = TickCount::now();
//         interval += t3 - t2;
//     }
//     TickCount t4 = TickCount::now();

//     TimeInterval tt;
//     tt = t4 - t1 - interval;
//     std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

//   return 0;
// }