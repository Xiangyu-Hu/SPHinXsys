/**
 * @file 	Viscous ball collision.cpp
 * @brief 	two balls with viscoplastic and oobleck(shear thickening) material bouncing with a boundary wall. 
 * @author 	Liezhao Wu, Xiaojing Tang and Xiangyu Hu
 */
#include "sphinxsys.h"  
#include "viscous_inelastic_solid.h"                 
using namespace SPH;                      
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------

Real DL = 0.8;                            /* boundary wall length. */
Real DH = 0.4;                            /* height. */
Real resolution_ref = DH  / 200.0;         /* reference resolution. */
Real BW = resolution_ref * 4.0;           /* boundary wall width. */
Real ball_radius = resolution_ref * 25.0; /* ball radius. */
Vec2d ball_1_center(0.3 * DL, 0.75 * DH); /* ball 1 center. */
Vec2d ball_2_center(0.7 * DL, 0.75 * DH); /* ball 2 center. */
BoundingBox system_domain_bounds(Vec2d(0.0, -BW), Vec2d(DL, DH));
// observer location
StdVec<Vecd> observation_location_1 = {ball_1_center};
StdVec<Vecd> observation_location_2 = {ball_2_center};
//----------------------------------------------------------------------
//	Global parameters on material properties.
//----------------------------------------------------------------------
Real gravity_g = 2.0;
Real rho0_s = 1.0e3;                /*  density */
Real Bulk_modulus = 1.09e5;         /*  bulk modulus */
Real Shear_modulus = 1.12e4;        /*  shear modulus */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);     
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); 
Real yield_stress = 0.1;          /* yield stress of plastic material. */                                  
Real viscosity = 10.0;             /* the viscosity. */  
Real Herschel_Bulkley_power_1 = 1.0;         /* viscoplastic material. */    
Real Herschel_Bulkley_power_2 = 2.8;         /* oobleck(shear thickening) material. */
//----------------------------------------------------------------------
//	Geometric shapes.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> wall_boundary_shape;
        wall_boundary_shape.push_back(Vecd(0.0, -BW));
        wall_boundary_shape.push_back(Vecd(0.0, 0.0));
        wall_boundary_shape.push_back(Vecd(DL, 0.0));
        wall_boundary_shape.push_back(Vecd(DL, -BW));
        wall_boundary_shape.push_back(Vecd(0.0, -BW));

        multi_polygon_.addAPolygon(wall_boundary_shape, ShapeBooleanOps::add);
    }
};
class BallOne : public MultiPolygonShape
{
  public:
    explicit BallOne(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(ball_1_center, ball_radius, 100, ShapeBooleanOps::add);
    }
};
class BallTwo : public MultiPolygonShape
{
  public:
    explicit BallTwo(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(ball_2_center, ball_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody ball_1(sph_system, makeShared<BallOne>("BallOne"));
    ball_1.defineBodyLevelSetShape();
    ball_1.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                    yield_stress, viscosity, Herschel_Bulkley_power_1);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_1.generateParticles<ParticleGeneratorReload>(io_environment, ball_1.getName())
        : ball_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody ball_2(sph_system, makeShared<BallTwo>("BallTwo"));
    ball_2.defineBodyLevelSetShape();
    ball_2.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                 yield_stress, viscosity, Herschel_Bulkley_power_2);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_2.generateParticles<ParticleGeneratorReload>(io_environment, ball_2.getName())
        : ball_2.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody ball_1_observer(sph_system, "BallOneObserver");
    ball_1_observer.generateParticles<ObserverParticleGenerator>(observation_location_1);
    ObserverBody ball_2_observer(sph_system, "BallTwoObserver");
    ball_2_observer.generateParticles<ObserverParticleGenerator>(observation_location_2);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_1_inner(ball_1);
        InnerRelation ball_2_inner(ball_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_1_random_particles(ball_1);
        SimpleDynamics<RandomizeParticlePosition> ball_2_random_particles(ball_2);
        relax_dynamics::RelaxationStepInner ball_1_relaxation_step_inner(ball_1_inner);
        relax_dynamics::RelaxationStepInner ball_2_relaxation_step_inner(ball_2_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_ball_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&ball_1, &ball_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_1_random_particles.exec(0.25);
        ball_2_random_particles.exec(0.25);
        write_ball_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_1_relaxation_step_inner.exec();
            ball_2_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_ball_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation ball_1_inner(ball_1);
    SurfaceContactRelation ball_1_contact(ball_1, {&wall_boundary});
    InnerRelation ball_2_inner(ball_2);
    SurfaceContactRelation ball_2_contact(ball_2, {&wall_boundary});
    ContactRelation ball_1_observer_contact(ball_1_observer, {&ball_1});
    ContactRelation ball_2_observer_contact(ball_2_observer, {&ball_2});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> ball_1_initialize_timestep(ball_1, gravity_ptr);
    SimpleDynamics<TimeStepInitialization> ball_2_initialize_timestep(ball_2, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_1_corrected_configuration(ball_1_inner);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_2_corrected_configuration(ball_2_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_1_get_time_step_size(ball_1, 0.2);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_2_get_time_step_size(ball_2, 1.0);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_1_stress_relaxation_first_half(ball_1_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_1_stress_relaxation_second_half(ball_1_inner);
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_2_stress_relaxation_first_half(ball_2_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_2_stress_relaxation_second_half(ball_2_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_1_update_contact_density(ball_1_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_1_compute_solid_contact_forces(ball_1_contact);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_2_update_contact_density(ball_2_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_2_compute_solid_contact_forces(ball_2_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        ball_1_position_recording("Position", io_environment, ball_1_observer_contact);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        ball_2_position_recording("Position", io_environment, ball_2_observer_contact);
    
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_1_corrected_configuration.exec();
    ball_2_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    ball_1_position_recording.writeToFile(0);
    ball_2_position_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    int ite = 0;
    Real End_time = 1.5;
    Real output_interval = 0.01 * End_time;
    Real Dt = 0.1 * output_interval;
    Real dt = 0.0;
    Real dt_1 = 0.0;
    Real dt_2 = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < End_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                ball_1_initialize_timestep.exec();
                ball_2_initialize_timestep.exec();
 
                ball_1_update_contact_density.exec();
                ball_1_compute_solid_contact_forces.exec();
                ball_1_stress_relaxation_first_half.exec(dt);
                ball_1_stress_relaxation_second_half.exec(dt);

                ball_1.updateCellLinkedList();
                ball_1_contact.updateConfiguration();

                ball_2_update_contact_density.exec();
                ball_2_compute_solid_contact_forces.exec();
                ball_2_stress_relaxation_first_half.exec(dt);
                ball_2_stress_relaxation_second_half.exec(dt);

                ball_2.updateCellLinkedList();
                ball_2_contact.updateConfiguration();

                ite++;
                dt_1 = ball_1_get_time_step_size.exec();
                dt_2 = ball_1_get_time_step_size.exec();
                dt = SMIN(dt_1, dt_2);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                if (ite % 1000 == 0)
                { 
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                              
            }
        }

        ball_1_position_recording.writeToFile(ite);
        ball_2_position_recording.writeToFile(ite);
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        ball_1_position_recording.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
        ball_2_position_recording.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
    }
    else
    {
        ball_1_position_recording.testResult();
        ball_2_position_recording.testResult();
    }
    return 0;
}
