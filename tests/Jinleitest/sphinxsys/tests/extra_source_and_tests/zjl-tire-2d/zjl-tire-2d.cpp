#include "sphinxsys.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
using namespace SPH;

//----------------------------------------------------------------------
// Geometry and resolution (2D)
//----------------------------------------------------------------------
Real resolution_ref = 0.0025;
Real R_outer = 0.316;
Real R_inner = 0.203;
Real water_length = 0.266;
Real water_height = 0.01;
Real road_length = 0.266;
Real road_height = 0.05;
Real gravity_g = 9.81;

Real rho0_s = 1200.0;//轮胎密度
Real Youngs_modulus = 1e5;
Real poisson = 0.45;
Real physical_viscosity = 1e6;
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real rho0_f = 1000.0;
Real c_f = 10.0;
Real mu_f = 0.001;

BoundingBox system_domain_bounds(Vec2d(-0.8, -0.5), Vec2d(1.3, 0.7));

// template <class GravityType>
// class GravityForce : public ForcePrior
// {
//   protected:
//     const GravityType gravity_;
//     Vecd *pos_;
//     Real *mass_;
//     Real *physical_time_;

//   public:
//     GravityForce(SPHBody &sph_body, const GravityType &gravity);
//     virtual ~GravityForce(){};
//     void update(size_t index_i, Real dt = 0.0);
// };
//SPHBody &sph_body,要变Bodypart
//----------------------------------------------------------------------
// Tire geometry (2D circular ring)
//----------------------------------------------------------------------
class Tire : public MultiPolygonShape {
public:
    explicit Tire(const std::string &shape_name) : MultiPolygonShape(shape_name) {
        multi_polygon_.addACircle(Vecd(water_length*0.5,road_height+water_height+R_outer ), R_outer, 50, ShapeBooleanOps::add);
        multi_polygon_.addACircle(Vecd(water_length*0.5, road_height+water_height+R_outer), R_inner, 50, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
// road (solid domain)
//----------------------------------------------------------------------
class Road : public MultiPolygonShape {
public:
    explicit Road(const std::string &shape_name) : MultiPolygonShape(shape_name) {
        std::vector<Vecd> road_shape = {
            Vecd(0.0, 0.0), Vecd(0.0, road_height),
            Vecd(road_length, road_height), Vecd(road_length, 0.0), Vecd(0.0, 0.0)
        };
        multi_polygon_.addAPolygon(road_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
// Water block (fluid domain)
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape {
public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name) {
        std::vector<Vecd> water_shape = {
            Vecd(0.0, road_height), Vecd(0.0, road_height+water_height),
            Vecd(water_length, road_height+water_height), Vecd(water_length, road_height), Vecd(0.0,road_height)
        };
        multi_polygon_.addAPolygon(water_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------
int main(int ac, char *av[]) {
    //--------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //--------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptationRatios(1.15,1.0);
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody tire(sph_system, makeShared<Tire>("Tire"));
    tire.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    tire.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
        (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? tire.generateParticles<BaseParticles, Reload>(tire.getName())
        : tire.generateParticles<BaseParticles, Lattice>();
    
    SolidBody road(sph_system, makeShared<Road>("Road"));
    road.defineMaterial<Solid>();
    road.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation insert_body_inner(tire);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(tire);
        RelaxationStepInner relaxation_step_inner(insert_body_inner);
        BodyStatesRecordingToVtp write_insert_body_to_vtp(tire);    
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_particle_reload_files(tire);
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
    InnerRelation insert_body_inner(tire);
    ContactRelation water_block_contact(water_block, {&road, &tire});
    ContactRelation insert_body_contact(tire, {&water_block});//是不是和两个对象接触就要加上RealBodyVector
    //ContactRelation beam_observer_contact(beam_observer, {&insert_body});
    //ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
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
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    Gravity gravity(Vec2d(0.0,-gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(tire, gravity);//改成部分，轮毂加速度
    
    
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(road);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(tire);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> insert_body_computing_time_step_size(tire);
    //BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createBeamBaseShape()));elasticSolid-shell里没有FixBodyPart
    //SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);
    SimpleDynamics<VonMisesStress> calculate_stress(tire);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_correction(DynamicsArgs(water_block_inner, 0.25), water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    //AlignedBoxPartByCell inflow_buffer(water_block, AlignedBox(xAxis, Transform(Vec2d(buffer_translation)), buffer_halfsize));
    //SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    //PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    //PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);

    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(tire);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(tire);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    //RegressionTestTimeAverage<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid(tire, "ViscousForceFromFluid");
    //RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_beam_tip_displacement("Position", beam_observer_contact);
    //ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    //periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing surface normal direction for the insert body. */
    insert_body_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    insert_body_corrected_configuration.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1;
    Real end_time = 0.1;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    // damping ratio for tire
    Real physical_viscosity_wall = 100.0;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
    wall_velocity_damping(0.6, insert_body_inner, "Velocity", physical_viscosity_wall);
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
            viscous_force.exec();
            transport_correction.exec();

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
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
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt - dt_s_sum);
                    insert_body_stress_relaxation_first_half.exec(dt_s);
                    //constraint_beam_base.exec();
                    wall_velocity_damping.exec(dt_s);
                    insert_body_stress_relaxation_second_half.exec(dt_s);
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
                //write_beam_tip_displacement.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            time_instance = TickCount::now();

            /** Water block configuration and periodic condition. */
            //periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            //periodic_condition.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            /** one need update configuration after periodic condition. */
            tire.updateCellLinkedList();
            insert_body_contact.updateConfiguration();
            boundary_indicator.exec();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        calculate_stress.exec();
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        //write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        //fluid_observer_contact.updateConfiguration();
        //write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

  return 0;
}
