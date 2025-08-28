#include "sphinxsys.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
#include "base_general_dynamics.h"
#include "base_data_package.h"
#include "carweightforce_prior.h"
#include "carweightforce_prior.hpp"
#include "tire_wheel_rim_material_definition.h"

using namespace SPH;
/**
 * @brief 2D Tire-water interaction simulation
 * @author Jinlei Zhou
 * @date 2025-06
 */

//----------------------------------------------------------------------
// Geometry and resolution (2D)
//----------------------------------------------------------------------
Real resolution_ref = 0.0025;
// Real R_outer = 0.300;
// Real R_inner = 0.200;
Real rim_thickness = 0.1;
// Real road_length = 0.400;
// Real road_height = 0.05;
Real gravity_ = 11.052;
Vec2d rim_upper_point(road_length*0.5,road_height+R_outer+0.5*rim_thickness);
//Real rho0_s = 1200.0;//轮胎密度 整体给还是轮毂和橡胶分开给？
Real Youngs_modulus = 1e4;//20-80MPa
//Real poisson = 0.45;
Real physical_viscosity = 1e6;
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real rho0_f = 1000.0;
Real c_f = 10.0;
Real mu_f = 0.001;

BoundingBox system_domain_bounds(Vec2d(-0.8, -0.5), Vec2d(1.3, 0.7));

//----------------------------------------------------------------------
// Tire geometry (2D circular ring)
//----------------------------------------------------------------------
class Tire : public MultiPolygonShape {
public:
    explicit Tire(const std::string &shape_name) : MultiPolygonShape(shape_name) {
        multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_outer, 50, ShapeBooleanOps::add);
        multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
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
// Wheel rim geometry (2D circular ring)
//----------------------------------------------------------------------
/** create the Wheel rim as constrain shape. */
MultiPolygon createRimBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner, 50, ShapeBooleanOps::add);
    multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
    return multi_polygon;
}
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

    //-------------------------------------------------------------------

    SolidBody tire(sph_system, makeShared<Tire>("Tire"));
    tire.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    //tire.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    tire.defineMaterial<TireBodyComposite>();
        (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? tire.generateParticles<BaseParticles, Reload>(tire.getName())
        : tire.generateParticles<BaseParticles, Lattice>();
    
    SolidBody road(sph_system, makeShared<Road>("Road"));
    road.defineMaterial<Solid>();
    road.generateParticles<BaseParticles, Lattice>();
        
    // SolidBody wheel_rim(sph_system, makeShared<crea>("WallBoundary"));
    // wheel_rim.defineMaterial<Solid>();
    // wheel_rim.generateParticles<BaseParticles, Lattice>();

    ObserverBody tire_observer(sph_system, "BallObserver");
    tire_observer.generateParticles<ObserverParticles>(StdVec<Vec2d>{rim_upper_point});

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
    InnerRelation insert_body_inner(tire);
    // ContactRelation tire_road_contact(tire, {&road});  
    // ContactRelation road_tire_contact(road, {&tire});应该是只能用于observe body和SPH body之间
    SurfaceContactRelation tire_road_contact(tire, {&road});
    //SurfaceContactRelation road_tire_contact(road, {&tire});需要反过来在定义一次吗
    //ContactRelation tire_observer_contact(tire_observer, {&tire});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    //车身重力仅作用于轮毂算出来的加速度  a=F/m=5000N/(Area of wheel rim*1*density of wheel) = 5000/（0.377*1200）=11.052
    BodyRegionByParticle beam_base(tire, makeShared<MultiPolygonShape>(createRimBaseShape()));
    Gravity gravity(Vec2d(0.0,-gravity_));
    SimpleDynamics<CarweightForce<Gravity>> constant_gravity(beam_base, gravity);//改成部分，轮毂加速度
    //SimpleDynamics<GravityForce<Gravity>> constant_gravity(tire, gravity);//改成部分，轮毂加速度

    
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(road);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(tire);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);
    SimpleDynamics<TireMaterialInitialization> composite_material_id(tire);
    /** stress relaxation for the tire. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
    //BodyRegionByParticle beam_base(tire, makeShared<MultiPolygonShape>(createRimBaseShape()));
    SimpleDynamics<VonMisesStress> calculate_stress(tire);
    /** Algorithms for tire-road contact. */
    InteractionDynamics<solid_dynamics::ContactFactorSummation> tire_update_contact_density(tire_road_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> tire_compute_solid_contact_forces(tire_road_contact);
    // damping ratio for tire
    Real physical_viscosity_wall = 1.0e8;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
        wall_velocity_damping(1.0, DynamicsArgs(insert_body_inner, "Velocity", physical_viscosity_wall));
    ReduceDynamics<solid_dynamics::AcousticTimeStep> tire_get_time_step_size(tire, 0.45);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    BodyStatesRecordingToVtp write_ball_state(tire);
    //ObservedQuantityRecording<Vecd>write_tire_displacement("Position", tire_observer_contact);
    write_real_body_states.addToWrite<Vecd>(tire, "Position");
    write_real_body_states.addToWrite<Real>(tire, "VonMisesStress");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    insert_body_normal_direction.exec();
    insert_body_corrected_configuration.exec();
    constant_gravity.exec();
    /** initialize material ids for the tire. */
    composite_material_id.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 1.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
    Real Dt = output_interval;
    Real dt = 0.0;
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
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: " << dt << "\n";
                }
                tire_update_contact_density.exec();
                tire_compute_solid_contact_forces.exec();
                insert_body_stress_relaxation_first_half.exec(dt);
                wall_velocity_damping.exec(dt);//速度约束两次前后各一次
                insert_body_stress_relaxation_second_half.exec(dt);
                tire.updateCellLinkedList();
                tire_road_contact.updateConfiguration();

                ite++;
                Real dt_free = tire_get_time_step_size.exec();
                dt = dt_free;
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            //write_tire_displacement.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        calculate_stress.exec();
        write_ball_state.writeToFile(ite);
        write_real_body_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    
  
    return 0;
}
