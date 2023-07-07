/**
 * @file 	collision.cpp
 * @brief 	two soft balls with and without internal damping bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0;                /**< box length. */
Real DH = 4.0;                /**< box height. */
Real resolution_ref = 0.025;  /**< reference resolution. */
Real BW = resolution_ref * 4; /**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center_1(2.0, 2.0);
Vec2d ball_center_2(6.0, 2.0);
Real ball_radius = 0.5;
// observer location
StdVec<Vecd> observation_location_1 = {ball_center_1};
StdVec<Vecd> observation_location_2 = {ball_center_2};
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real gravity_g = 1.0;
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
Real physical_viscosity = 10000.0;
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, -BW));

        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(0.0, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, DH));
        inner_wall_shape.push_back(Vecd(DL, DH));
        inner_wall_shape.push_back(Vecd(DL, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};
class FreeBall : public MultiPolygonShape
{
  public:
    explicit FreeBall(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
    }
};
class DampingBall : public MultiPolygonShape
{
  public:
    explicit DampingBall(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(ball_center_2, ball_radius, 100, ShapeBooleanOps::add);
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
    SolidBody free_ball(sph_system, makeShared<FreeBall>("FreeBall"));
    free_ball.defineBodyLevelSetShape();
    free_ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? free_ball.generateParticles<ParticleGeneratorReload>(io_environment, free_ball.getName())
        : free_ball.generateParticles<ParticleGeneratorLattice>();

    SolidBody damping_ball(sph_system, makeShared<DampingBall>("DampingBall"));
    damping_ball.defineBodyLevelSetShape();
    damping_ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? damping_ball.generateParticles<ParticleGeneratorReload>(io_environment, damping_ball.getName())
        : damping_ball.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody free_ball_observer(sph_system, "FreeBallObserver");
    free_ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_1);
    ObserverBody damping_ball_observer(sph_system, "DampingBallObserver");
    damping_ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_2);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation free_ball_inner(free_ball);
        InnerRelation damping_ball_inner(damping_ball);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> free_ball_random_particles(free_ball);
        SimpleDynamics<RandomizeParticlePosition> damping_ball_random_particles(damping_ball);
        relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
        relax_dynamics::RelaxationStepInner damping_ball_relaxation_step_inner(damping_ball_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_ball_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&free_ball, &damping_ball});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        free_ball_random_particles.exec(0.25);
        damping_ball_random_particles.exec(0.25);
        write_ball_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            free_ball_relaxation_step_inner.exec();
            damping_ball_relaxation_step_inner.exec();
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
    InnerRelation free_ball_inner(free_ball);
    SurfaceContactRelation free_ball_contact(free_ball, {&wall_boundary});
    InnerRelation damping_ball_inner(damping_ball);
    SurfaceContactRelation damping_ball_contact(damping_ball, {&wall_boundary});
    ContactRelation free_ball_observer_contact(free_ball_observer, {&free_ball});
    ContactRelation damping_all_observer_contact(damping_ball_observer, {&damping_ball});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> free_ball_initialize_timestep(free_ball, gravity_ptr);
    SimpleDynamics<TimeStepInitialization> damping_ball_initialize_timestep(damping_ball, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> free_ball_corrected_configuration(free_ball_inner);
    InteractionWithUpdate<CorrectedConfigurationInner> damping_ball_corrected_configuration(damping_ball_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> free_ball_get_time_step_size(free_ball);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> damping_ball_get_time_step_size(damping_ball);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> free_ball_stress_relaxation_first_half(free_ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> free_ball_stress_relaxation_second_half(free_ball_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> damping_ball_stress_relaxation_first_half(damping_ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> damping_ball_stress_relaxation_second_half(damping_ball_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> free_ball_update_contact_density(free_ball_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> free_ball_compute_solid_contact_forces(free_ball_contact);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> damping_ball_update_contact_density(damping_ball_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> damping_ball_compute_solid_contact_forces(damping_ball_contact);
    /** Damping for one ball */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        damping(0.5, damping_ball_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        free_ball_displacement_recording("Position", io_environment, free_ball_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        damping_ball_displacement_recording("Position", io_environment, damping_all_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    free_ball_corrected_configuration.exec();
    damping_ball_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    free_ball_displacement_recording.writeToFile(0);
    damping_ball_displacement_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 10.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
    Real Dt = 0.1 * output_interval;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                free_ball_initialize_timestep.exec();
                damping_ball_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                free_ball_update_contact_density.exec();
                free_ball_compute_solid_contact_forces.exec();
                free_ball_stress_relaxation_first_half.exec(dt);
                free_ball_stress_relaxation_second_half.exec(dt);

                free_ball.updateCellLinkedList();
                free_ball_contact.updateConfiguration();

                damping_ball_update_contact_density.exec();
                damping_ball_compute_solid_contact_forces.exec();
                damping_ball_stress_relaxation_first_half.exec(dt);
                damping.exec(dt);
                damping_ball_stress_relaxation_second_half.exec(dt);

                damping_ball.updateCellLinkedList();
                damping_ball_contact.updateConfiguration();

                ite++;
                Real dt_free = free_ball_get_time_step_size.exec();
                Real dt_damping = damping_ball_get_time_step_size.exec();
                dt = SMIN(dt_free, dt_damping);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                free_ball_displacement_recording.writeToFile(ite);
                damping_ball_displacement_recording.writeToFile(ite);
            }
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    free_ball_displacement_recording.testResult();
    damping_ball_displacement_recording.testResult();

    return 0;
}
