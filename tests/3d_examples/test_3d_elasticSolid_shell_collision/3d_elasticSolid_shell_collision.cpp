/**
 * @file 	3d_elasticSolid_shell_collision.cpp
 * @brief 	This is a benchmark test of the 3D elastic solid->shell contact/impact formulations.
 * @details  We consider the collision of an elastic ball bouncing in a spherical shell box.
 * @author 	Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.05;                         /**< reference resolution. */
Real thickness = resolution_ref * 1.;               /**< shell thickness. */
Real radius = 2.0;                                  /**< cylinder radius. */
Real half_height = 1.0;                             /** Height of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0; /** Radius of the mid surface. */
Vec3d ball_center(radius / 2.0, 0.0, 0.0);
int particle_number_mid_surface = int(2.0 * radius_mid_surface * Pi * 215.0 / 360.0 / resolution_ref);
int particle_number_height = 2 * int(half_height / resolution_ref);
int BWD = 1; /** Width of the boundary layer measured by number of particles. */
BoundingBox system_domain_bounds(Vec3d(-radius - thickness, -half_height - thickness, -radius - thickness),
                                 Vec3d(radius + thickness, half_height + thickness, radius + thickness));
StdVec<Vecd> ball_observation_location = {ball_center};
Real ball_radius = 0.5;
Real gravity_g = 1.0;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 2.0e4;
Real poisson = 0.45;
Real physical_viscosity = 1.0e6;
//----------------------------------------------------------------------
/** Define application dependent particle generator for thin structure. */
class CylinderParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit CylinderParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        // the cylinder and boundary
        for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
        {
            for (int j = 0; j < particle_number_height; j++)
            {
                Real x = radius_mid_surface * cos(162.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
                Real y = (j - particle_number_height / 2) * resolution_ref + resolution_ref * 0.5;
                Real z = radius_mid_surface * sin(162.5 / 180.0 * Pi + (i - BWD + 0.5) * 215.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
                initializePositionAndVolumetricMeasure(Vecd(x, y, z), resolution_ref * resolution_ref);
                Vec3d n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
                initializeSurfaceProperties(n_0, thickness);
            }
        }
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
    // sph_system.generate_regression_data_ = true;
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    /** create a shell body. */
    SolidBody shell(sph_system, makeShared<DefaultShape>("shell"));
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    shell.generateParticles<CylinderParticleGenerator>();

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

    ObserverBody ball_observer(sph_system, "BallObserver");
    ball_observer.generateParticles<ObserverParticleGenerator>(ball_observation_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner(ball);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, ball);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles.exec(0.25);
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation ball_inner(ball);
    SurfaceContactRelation ball_contact(ball, {&shell});
    ContactRelation ball_observer_contact(ball_observer, {&ball});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> ball_initialize_timestep(ball, makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g)));
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration(ball_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size(ball, 0.45);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> ball_stress_relaxation_first_half(ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half(ball_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ShellContactDensity> ball_update_contact_density(ball_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces(ball_contact);
    DampingWithRandomChoice<InteractionSplit<solid_dynamics::PairwiseFrictionFromWall>>
        ball_friction(0.1, ball_contact, physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    BodyStatesRecordingToVtp write_ball_state(io_environment, {&ball});
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
                ball_friction.exec(dt);
                ball_stress_relaxation_second_half.exec(dt);

                ball.updateCellLinkedList();
                ball_contact.updateConfiguration();

                ite++;
                Real dt_free = ball_get_time_step_size.exec();
                dt = dt_free;
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            write_ball_center_displacement.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        write_ball_state.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.generate_regression_data_)
    {
        write_ball_center_displacement.generateDataBase(0.005);
    }
    else
    {
        write_ball_center_displacement.testResult();
    }
    return 0;
}
