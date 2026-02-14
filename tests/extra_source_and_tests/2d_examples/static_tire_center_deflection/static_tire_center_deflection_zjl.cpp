/**
 * @file  static_tire_center_deflection_zjl.cpp
 *  */
#include "sphinxsys.h"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "pressure_boundary.h"
#include "base_general_dynamics.h"
#include "base_data_type.h"
#include "base_data_type_package.h"
#include "carweightforce_prior.h"
#include "carweightforce_prior.hpp"
#include "tire_wheel_rim_material_definition.h"
#include <filesystem>
#include <fstream>
using namespace SPH;

// Geometry parameters
Real resolution_ref = 0.0025;
Real rim_thickness = 0.1;
Real acceleration = 200.0;
Vec2d tire_center(road_length * 0.5, road_height + R_outer);
Vec2d tire_top_point = tire_center + Vec2d(0.0, R_outer-10*resolution_ref);  

Real physical_viscosity = 1e6;                                      
Real rho0_f = 1000.0;
BoundingBoxd system_domain_bounds(Vec2d(-0.8, -0.5), Vec2d(1.3, 0.7));

// Tire geometry (2D circular ring)
class Tire : public MultiPolygonShape {
public:
    explicit Tire(const std::string &shape_name) : MultiPolygonShape(shape_name) {
        multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_outer, 50, ShapeBooleanOps::add);
        multi_polygon_.addACircle(Vecd(road_length*0.5,road_height+R_outer), R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
    }
};

// Road geometry
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

// Wheel rim geometry (2D circular ring)
MultiPolygon createRimBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner, 50, ShapeBooleanOps::add);
    multi_polygon.addACircle(Vecd(road_length*0.5,road_height+R_outer ),R_inner-rim_thickness, 50, ShapeBooleanOps::sub);
    return multi_polygon;
}

// ============================
// Main program
// ============================
int main(int ac, char *av[]) {
    //--------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    
    // sph_system.setRunParticleRelaxation(true);
    // sph_system.setReloadParticles(false);
    
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.setGenerateRegressionData(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    // ========= Bodies =========
    
    // Tire particles
    SolidBody tire(sph_system, makeShared<Tire>("Tire"));
    tire.defineBodyLevelSetShape()->writeLevelSet();
    tire.defineMaterial<TireBodyComposite>();
        (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? tire.generateParticles<BaseParticles, Reload>(tire.getName())
        : tire.generateParticles<BaseParticles, Lattice>();
    
    // Road particles
    SolidBody road(sph_system, makeShared<Road>("Road"));
    road.defineMaterial<Solid>();
    road.generateParticles<BaseParticles, Lattice>();

    // Deflection observer
    ObserverBody tire_observer(sph_system, "TireObserver");
    tire_observer.generateParticles<ObserverParticles>(StdVec<Vec2d>{tire_top_point});

    // ====================================================================================
    //	Run particle relaxation for body-fitted distribution if chosen.
    // ====================================================================================
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
        
        // ========= Output results =========
        write_particle_reload_files.writeToFile(0);
        write_real_body_states.writeToFile(0);
        return 0;
    }
    
    // ========= Relations =========
    InnerRelation insert_body_inner(tire);
    SurfaceContactRelation tire_road_contact(tire, {&road});
    BodyRegionByParticle beam_base(tire, makeShared<MultiPolygonShape>(createRimBaseShape()));
    Gravity gravity(Vec2d(0.0,-acceleration));
    SimpleDynamics<CarweightForce<Gravity>> constant_gravity(beam_base, gravity);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(road);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(tire);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);
    SimpleDynamics<TireMaterialInitialization> composite_material_id(tire);
    ContactRelation tire_observer_contact(tire_observer, {&tire});
    // ========= stress relaxation for the tire =========
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
    SimpleDynamics<VonMisesStress> calculate_stress(tire);
    
    // ========= tire-road contact =========
    InteractionDynamics<solid_dynamics::ContactFactorSummation> tire_update_contact_density(tire_road_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> tire_compute_solid_contact_forces(tire_road_contact);
    
    // ========= damping =========
    Real physical_viscosity_tire = 1000;
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
    tire_velocity_damping(0.6, DynamicsArgs(insert_body_inner, "Velocity", physical_viscosity_tire));
    ReduceDynamics<solid_dynamics::AcousticTimeStep> tire_get_time_step_size(tire, 0.45);
    
    // ========= IO =========
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    BodyStatesRecordingToVtp write_tire_state(tire);
    write_real_body_states.addToWrite<Vecd>(tire, "Position");
    write_real_body_states.addToWrite<Real>(tire, "VonMisesStress");
    write_real_body_states.addToWrite<int>(tire, "MaterialID");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
    write_rim_upper_position("Position", tire_observer_contact);

    // ========= Init =========
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    insert_body_normal_direction.exec();
    insert_body_corrected_configuration.exec();
    constant_gravity.exec();
    composite_material_id.exec();
    
    // ========= Time stepping =========
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 0.5;
    Real end_time = T0;
    Real output_interval = T0/50;
    Real Dt = output_interval;
    Real dt = 0.0;
    
    // ========= Statistics for CPU time =========
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    
    // ========= First output before the main loop =========
    write_real_body_states.writeToFile(0);
    tire_observer_contact.updateConfiguration();
    write_rim_upper_position.writeToFile(0);

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
                tire_velocity_damping.exec(dt);
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

        }
        TickCount t2 = TickCount::now();
        calculate_stress.exec();
        write_tire_state.writeToFile();
        write_real_body_states.writeToFile();
        tire_observer_contact.updateConfiguration();
        write_rim_upper_position.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    
    if (sph_system.GenerateRegressionData())
    {
        write_rim_upper_position.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_rim_upper_position.testResult();
    }
  
    return 0;
}
