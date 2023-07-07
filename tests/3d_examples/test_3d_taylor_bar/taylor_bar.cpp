/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar.
 * @author Chi Zhang  and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 */
#include "sphinxsys.h"

#include "taylor_bar.h" /**< Case setup for this example. */

using namespace SPH;

int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    system.setRunParticleRelaxation(false);
    system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(system);

    /** create a body with corresponding material, particles and reaction model. */
    SolidBody column(system, makeShared<Column>("Column"));
    column.defineAdaptationRatios(1.3, 1.0);
    column.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    column.defineParticlesAndMaterial<ElasticSolidParticles, HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? column.generateParticles<ParticleGeneratorReload>(io_environment, column.getName())
        : column.generateParticles<ParticleGeneratorLattice>();
    column.addBodyStateForRecording<Vecd>("NormalDirection");

    SolidBody wall(system, makeShared<Wall>("Wall"));
    wall.defineParticlesAndMaterial<SolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall.generateParticles<ParticleGeneratorLattice>();

    /** Define Observer. */
    ObserverBody my_observer(system, "MyObserver");
    my_observer.generateParticles<ColumnObserverParticleGenerator>();

    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    SurfaceContactRelation column_wall_contact(column, {&wall});
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);

    if (system.RunParticleRelaxation())
    {
        /**
         * @brief 	Methods used for particle relaxation.
         */
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_column_particles(column);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_column_to_vtp(io_environment, column);
        /** Write the particle reload files. */

        ReloadParticleIO write_particle_reload_files(io_environment, column);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(column_inner);
        /**
         * @brief 	Particle relaxation starts here.
         */
        random_column_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_states.writeToFile(0.0);

        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
                write_column_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0.0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_relaxation_first_half(column_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(column_inner);
    InteractionDynamics<solid_dynamics::DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    SimpleDynamics<InitialCondition> initial_condition(column);
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(column_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(column, 0.3);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_velocity("Velocity", io_environment, my_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", io_environment, my_observer_contact);

    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    corrected_configuration.exec();
    initial_condition.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    int ite = 0;
    Real end_time = 1.0e-4;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real output_period = 1.0e-6; // anyway 50 write_states files in total
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_displacement.writeToFile(0);
    write_velocity.writeToFile(0);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % screen_output_interval == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";

                if (ite != 0 && ite % observation_sample_interval == 0)
                {
                    write_displacement.writeToFile(ite);
                    write_velocity.writeToFile(ite);
                }
            }
            column_wall_contact_force.exec(dt);
            stress_relaxation_first_half.exec(dt);
            stress_relaxation_second_half.exec(dt);

            column.updateCellLinkedList();
            column_wall_contact.updateConfiguration();

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    write_displacement.testResult();
    write_velocity.testResult();

    return 0;
}
