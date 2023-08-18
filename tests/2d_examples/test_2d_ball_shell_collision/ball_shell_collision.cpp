/**
 * @file 	ball_shell_collision.cpp
 * @brief 	an elastic ball bouncing within a confined shell boundary
 * @details This is a case to test elasticSolid -> shell impact/collision.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.025; /**< reference resolution. */
Real circle_radius = 2.0;
Vec2d circle_center(2.0, 2.0);
Real thickness = resolution_ref * 1.; /**< shell thickness. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
BoundingBox system_domain_bounds(Vec2d(-thickness, -thickness), Vec2d(2.0 * circle_radius + thickness, 2.0 * circle_radius + thickness));
Vec2d ball_center(3.0, 1.5);
Real ball_radius = 0.5;
StdVec<Vecd> beam_observation_location = {ball_center};
Real gravity_g = 1.0;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 2.0e4;
Real poisson = 0.45;
Real physical_viscosity = 1.0e6;
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(circle_center, circle_radius + resolution_ref);
        subtract<GeometricShapeBall>(circle_center, circle_radius);
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
    SolidBody ball(sph_system, makeShared<GeometricShapeBall>(ball_center, ball_radius, "BallBody"));
    ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName());
    }
    else
    {
        ball.defineBodyLevelSetShape()->writeLevelSet(io_environment);
        ball.generateParticles<ParticleGeneratorLattice>();
    }

    // Note the wall boundary here has sharp corner, and is a numerical invalid elastic shell structure,
    // and its dynamics is not able to be modeled by the shell dynamics in SPHinXsys in the current version.
    // Here, we use it simply as a rigid shell.
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    // here dummy linear elastic solid is use because no solid dynamics in particle relaxation
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName());
    }
    else
    {
        wall_boundary.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(io_environment);
        wall_boundary.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness);
    }

    if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }

    ObserverBody ball_observer(sph_system, "BallObserver");
    ball_observer.generateParticles<ObserverParticleGenerator>(beam_observation_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner(ball);
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for wall boundary.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> wall_boundary_random_particles(wall_boundary);
        relax_dynamics::ShellRelaxationStepInner
            relaxation_step_wall_boundary_inner(wall_boundary_inner, thickness, level_set_refinement_ratio);
        relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(wall_boundary_inner, thickness, cos(Pi / 3.75));
        wall_boundary.addBodyStateForRecording<int>("UpdatedIndicator");
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        MeshRecordingToPlt write_mesh_cell_linked_list(io_environment, wall_boundary.getCellLinkedList());
        ReloadParticleIO write_particle_reload(io_environment, {&ball, &wall_boundary});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles.exec(0.25);
        wall_boundary_random_particles.exec(0.25);

        relaxation_step_wall_boundary_inner.mid_surface_bounding_.exec();
        write_relaxed_particles.writeToFile(0);
        wall_boundary.updateCellLinkedList();
        write_mesh_cell_linked_list.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step_inner.exec();
            for (int k = 0; k < 2; ++k)
                relaxation_step_wall_boundary_inner.exec();
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
    //----------------------------------------------------------------------
    InnerRelation ball_inner(ball);
    SurfaceContactRelation ball_contact(ball, {&wall_boundary});
    ContactRelation ball_observer_contact(ball_observer, {&ball});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Define external force.*/
    SimpleDynamics<TimeStepInitialization> ball_initialize_timestep(ball, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration(ball_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size(ball);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> ball_stress_relaxation_first_half(ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half(ball_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ShellContactDensity> ball_update_contact_density(ball_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces(ball_contact);
    // DampingWithRandomChoice<InteractionSplit<solid_dynamics::PairwiseFrictionFromWall>>
    // 	ball_friction(0.1, ball_contact, physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_ball_center_displacement("Position", io_environment, ball_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_corrected_configuration.exec();

    /** Initial states output. */
    body_states_recording.writeToFile(0);
    /** Main loop. */
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
                ball_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                ball_update_contact_density.exec();
                ball_compute_solid_contact_forces.exec();
                ball_stress_relaxation_first_half.exec(dt);
                // ball_friction.exec(dt);
                ball_stress_relaxation_second_half.exec(dt);

                ball.updateCellLinkedList();
                ball_contact.updateConfiguration();

                ite++;
                Real dt_ball = ball_get_time_step_size.exec();
                dt = dt_ball;
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
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

    if (sph_system.generate_regression_data_)
    {
        write_ball_center_displacement.generateDataBase(1.0e-2);
    }
    else
    {
        write_ball_center_displacement.testResult();
    }

    return 0;
}
