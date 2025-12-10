/**
 * @file ball_shell_collision.cpp
 * @brief an elastic ball bouncing within a rigid shell boundary
 * @details This is a case to test elasticSolid -> shell impact/collision.
 * @author Massoud Rezavand and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.025;
Real shell_shape_radius = 2.0;
Vec2d shell_shape_center(2.0, 2.0);
Real thickness = resolution_ref * 1.;
Vec2d ball_center(3.0, 1.5);
Real ball_radius = 0.5;
Vec2d gravity(0.0, -1.0);
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 2.0e4;
Real poisson = 0.45;
//----------------------------------------------------------------------
//	Case dependent geometries.
//----------------------------------------------------------------------
class ShellShape : public ComplexShape
{
  public:
    explicit ShellShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(shell_shape_center, shell_shape_radius + thickness);
        subtract<GeometricShapeBall>(shell_shape_center, shell_shape_radius);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem and IO environment.
    //----------------------------------------------------------------------
    Vec2d domain_lower_bound(-thickness, -thickness);
    Real domain_box_size = 2.0 * shell_shape_radius + thickness;
    Vec2d domain_upper_bound(domain_box_size, domain_box_size);
    BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(true);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody ball(sph_system, makeShared<GeometricShapeBall>(ball_center, ball_radius, "BallBody"));
    ball.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        ball.generateParticles<BaseParticles, Reload>(ball.getName());
    }
    else
    {
        ball.defineBodyLevelSetShape()->writeLevelSet(sph_system);
        ball.generateParticles<BaseParticles, Lattice>();
    }

    SolidBody rigid_shell(sph_system, makeShared<ShellShape>("ShellShape"));
    rigid_shell.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    rigid_shell.defineMaterial<Solid>();
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        rigid_shell.generateParticles<SurfaceParticles, Reload>(rigid_shell.getName());
    }
    else if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }
    else
    {
        Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
        rigid_shell.defineBodyLevelSetShape(level_set_refinement_ratio, UsageType::Surface)
            ->writeLevelSet(sph_system);
        rigid_shell.generateParticles<SurfaceParticles, Lattice>(thickness);
    }

    ObserverBody ball_observer(sph_system, "BallObserver");
    ball_observer.generateParticles<ObserverParticles>(StdVec<Vecd>{ball_center});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner(ball);
        InnerRelation rigid_shell_inner(rigid_shell);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for the ball.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
        RelaxationStepInner ball_relaxation_step(ball_inner);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for the rigid shell.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> rigid_shell_random_particles(rigid_shell);
        ShellRelaxationStep rigid_shell_relaxation_step(rigid_shell_inner);
        ShellNormalDirectionPrediction shell_normal_prediction(rigid_shell_inner, thickness, cos(Pi / 3.75));
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(sph_system);
        write_relaxed_particles.addToWrite<int>(rigid_shell, "UpdatedIndicator");
        ReloadParticleIO write_particle_reload({&ball, &rigid_shell});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles.exec(0.25);
        ball_relaxation_step.SurfaceBounding().exec();
        rigid_shell_random_particles.exec(0.25);
        rigid_shell_relaxation_step.MidSurfaceBounding().exec();
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step.exec();
            for (int k = 0; k < 2; ++k)
                rigid_shell_relaxation_step.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        shell_normal_prediction.exec();
        write_relaxed_particles.writeToFile(ite);
        write_particle_reload.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation ball_inner(ball);
    ShellSurfaceContactRelation ball_contact(ball, {&rigid_shell});
    ContactRelation ball_observer_contact(ball_observer, {&ball});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    Gravity constant_gravity(gravity);
    SimpleDynamics<GravityForce<Gravity>> ball_constant_gravity(ball, constant_gravity);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> ball_corrected_configuration(ball_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> ball_stress_relaxation_first_half(ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half(ball_inner);
    InteractionDynamics<solid_dynamics::ShellContactFactor> ball_update_contact_density(ball_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces(ball_contact);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> ball_get_time_step_size(ball);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_ball_center_displacement("Position", ball_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_corrected_configuration.exec();
    ball_constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
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
                ball_update_contact_density.exec();
                ball_compute_solid_contact_forces.exec();
                ball_stress_relaxation_first_half.exec(dt);
                ball_stress_relaxation_second_half.exec(dt);

                ball.updateCellLinkedList();
                ball_contact.updateConfiguration();

                ite++;
                Real dt_ball = ball_get_time_step_size.exec();
                dt = dt_ball;
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            write_ball_center_displacement.writeToFile(ite);
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

    if (sph_system.GenerateRegressionData())
    {
        write_ball_center_displacement.generateDataBase(1.0e-2);
    }
    else
    {
        write_ball_center_displacement.testResult();
    }

    return 0;
}
