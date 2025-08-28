/**
 * @file 	poiseuille_flow_shell.cpp
 * @brief 	3D poiseuille flow interaction with shell example
 * @details This is the one of the basic test cases for validating fluid-rigid shell interaction
 * @author  Weiyi Kong
 */

#include <gtest/gtest.h>
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
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e5;    /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae;
Real DH = 0.001;                                             /**< Channel height. */
Real resolution_ref = DH / 15.0;                             /**< Initial reference particle spacing. */
Real diameter = 0.001;
Real fluid_radius = 0.5 * diameter;//0.0005
Real outer_radius = fluid_radius+4*resolution_ref ;      // 外径（血管外壁半径）

Real full_length = 0.004;   // 原来的 DL = 0.004
const BoundingBox system_domain_bounds(
    Vec3d(-0.005,             
          -fluid_radius - 0.001,
          -fluid_radius - 0.001),
    Vec3d(full_length + 0.005, 
          fluid_radius + 0.001,
          fluid_radius + 0.001));

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.0;
Real Outlet_pressure = 0.0;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = 4.0e-3;
Real U_f = 0.01;
Real c_f = 10.0 * U_f;
Real rho = 1265; // kg/m^3



//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec3d bidirectional_buffer_halfsize(2 * resolution_ref,
                                    fluid_radius ,
                                    fluid_radius);

Vec3d left_bidirectional_translation(2 * resolution_ref, 0.0, 0.0); // 圆柱起点
Vec3d right_bidirectional_translation(full_length-2 * resolution_ref, 0.0, 0.0); // 圆柱终点


//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
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
        target_velocity[0] = 0.01;
        target_velocity[1] = 0.0;
        target_velocity[2] = 0.0;
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//  Define Water Shape
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
      public:
      explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
      {
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), fluid_radius, 0.5 * full_length, 50, Vec3d(0.5 * full_length,0.0, 0.0));
      }
};

//----------------------------------------------------------------------
//  Define Vessel Shape
//----------------------------------------------------------------------
class Vessel : public ComplexShape
{
      public:
      explicit Vessel(const std::string &shape_name) : ComplexShape(shape_name)
      {
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), outer_radius, 0.5 * full_length, 50, Vec3d(0.5 * full_length,0.0, 0.0));
        subtract<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1., 0., 0.), fluid_radius, 0.5 * full_length, 50, Vec3d(0.5 * full_length,0.0, 0.0));
      }
};



// template <>
// class ParticleGenerator<SurfaceParticles, Vessel> : public ParticleGenerator<SurfaceParticles>
// {
//      Real resolution_vessel_;
//      Real vessel_wall_half_thickness_;
//      Real vessel_surface_thickness_;


//   public:
//     explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
//                                Real resolution_vessel_,
//                                Real vessel_wall_half_thickness_,
//                                Real vessel_surface_thickness_ )
//         : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
//           resolution_vessel_(resolution_vessel),
//           vessel_wall_half_thickness_(vessel_wall_half_thickness),
//           vessel_surface_thickness_(vessel_surface_thickness) {}

//     void prepareGeometricData() override
//     {
//          const Real radius_mid_surface = (fluid_radius + 0.0007) * 0.5;
//          const int particle_number_circumference = int(2.0 * radius_mid_surface * Pi / resolution_vessel_);
//          const int particle_number_length = int((full_length + 2.0 * vessel_wall_half_thickness_) / resolution_vessel_);
//          const int particle_number_radial = 4;  // 想要的层数
//          const Real inner_radius = fluid_radius;
//          const Real outer_radius = 0.0007;  // 定义的外径
//          const Real dr = (outer_radius - inner_radius) / Real(particle_number_radial);

//         for (int k = 0; k < particle_number_radial; ++k)
//       {
//          Real radius = inner_radius + 0.5 * dr + Real(k) * dr;

//           for (int i = 0; i < particle_number_circumference; ++i)
//        {
//            for (int j = 0; j < particle_number_length; ++j)
//           {
//             Real theta = (i + 0.5) * 2.0 * Pi / Real(particle_number_circumference);
//             Real x = -vessel_wall_half_thickness_ + (full_length + 2.0 * vessel_wall_half_thickness_) * j / Real(particle_number_length) + 0.5 * resolution_vessel_;
//             Real y = radius * cos(theta);
//             Real z = radius * sin(theta);

//             Vec3d pos(x, y, z);
//             Vec3d normal(0.0, y / radius, z / radius);

//             this->addPositionAndVolumetricMeasure(pos, resolution_vessel_ * resolution_vessel_);
//             this->addSurfaceProperties(normal, vessel_surface_thickness_);
//            }
//         }
//       }
//     }
// };




//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);         // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(4.0);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);
    

    SolidBody vessel_boundary(sph_system, makeShared<Vessel>("Vessel"));
    vessel_boundary.defineAdaptationRatios(1.15,1.0);
    vessel_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    vessel_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
     (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? vessel_boundary.generateParticles<BaseParticles, Reload>(vessel_boundary.getName())
        : vessel_boundary.generateParticles<BaseParticles, Lattice>();
   
    // vessel_boundary.generateParticles<SurfaceParticles,Vessel>(
    // resolution_vessel, vessel_wall_half_thickness, vessel_surface_thickness); surfaceparticle 是用来生成shell的
    
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation insert_body_inner(vessel_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(vessel_boundary);
        RelaxationStepInner relaxation_step_inner(insert_body_inner);
        BodyStatesRecordingToVtp write_insert_body_to_vtp(vessel_boundary);    
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_particle_reload_files(vessel_boundary);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_insert_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_insert_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        write_real_body_states.writeToFile(0);
        return 0;
    }

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation vessel_boundary_inner(vessel_boundary);
    ContactRelation water_block_contact(water_block, {&vessel_boundary});
    ContactRelation insert_body_contact(vessel_boundary, {&water_block});
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
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> vessel_boundary_corrected_configuration(vessel_boundary_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> vessel_boundary_stress_relaxation_first_half(vessel_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> vessel_boundary_stress_relaxation_second_half(vessel_boundary_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> vessel_boundary_computing_time_step_size(vessel_boundary);
    // 左端固定区域
    Vecd half_size_left(4*resolution_ref, 0.001, 0.001);  // 0.001 应该大于 fluid_radius
    Vecd center_left(0.0, 0.0, 0.0);            // x=0 为圆柱左端

    BodyRegionByParticle holder_left(vessel_boundary, 
    makeShared<TriangleMeshShapeBrick>(half_size_left, 1, center_left));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder_left(holder_left);

    // 右端固定区域
    Vecd half_size_right = half_size_left;
    Vecd center_right(full_length-4*resolution_ref, 0.0, 0.0);   // x=full_length 为右端中心

    BodyRegionByParticle holder_right(vessel_boundary, 
    makeShared<TriangleMeshShapeBrick>(half_size_right, 1, center_right));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder_right(holder_right);
    
    SimpleDynamics<VonMisesStress> calculate_stress(vessel_boundary);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

    AlignedBoxPartByCell left_emitter(water_block, AlignedBox(xAxis, Transform(Vec3d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    AlignedBoxPartByCell right_emitter(water_block, AlignedBox(xAxis, Transform(Rotation3d(Pi,Vec3d(0.0,0.0,1.0)),right_bidirectional_translation), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    //AlignedBoxPartByCell right_emitter(water_block, makeShared <AlignedBox> ( xAxis,Transform (Rotation3d(Pi, Vec3d(0.0, 0.0, 1.0)),right_bidirectional_translation),bidirectional_buffer_halfsize));

    //InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_density_by_summation(water_block_inner, water_block_contact);
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(vessel_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(vessel_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);
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
    body_states_recording.addToWrite<Vecd>(vessel_boundary, "Velocity");
    body_states_recording.addToWrite<Matd>(vessel_boundary, "DeformationGradient");
    body_states_recording.addToWrite<Real>(vessel_boundary, "VonMisesStress");
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
    vessel_boundary_normal_direction.exec();//这里有问题 跑到这就跑不出来了 应该是之前生成的是shell 所以这里算不出来
    vessel_boundary_corrected_configuration.exec();//跳过上面的跑这个没问题
    // //----------------------------------------------------------------------
    // //	Setup for time-stepping control
    // //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1;
    Real end_time = 0.1;
    Real output_interval = end_time / 500.0;
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
    Real physical_viscosity_wall = 100.0;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd, FixedDampingRate>>>
    wall_velocity_damping(0.6, vessel_boundary_inner, "Velocity", physical_viscosity_wall);
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
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_correction.exec();
            /** Update normal direction on elastic body.*/
            insert_body_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                // density_relaxation.exec(dt);

                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);
                
                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(vessel_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    //vessel_boundary_stress_relaxation_first_half.exec(dt_s);
                    //constraint_holder_left.exec(dt_s);
                    //constraint_holder_right.exec(dt_s);
                    //wall_velocity_damping.exec(dt_s);
                    //constraint_holder_left.exec(dt_s);
                    //constraint_holder_right.exec(dt_s);
                    //vessel_boundary_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
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
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            vessel_boundary.updateCellLinkedList();
            insert_body_contact.updateConfiguration();
            // interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }

        TickCount t2 = TickCount::now();

        /** write run-time observation into file */
        calculate_stress.exec();
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
