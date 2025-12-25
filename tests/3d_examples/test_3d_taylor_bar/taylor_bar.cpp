/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar.
 * @author Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 * //TODO: Seems that the wall contact force should be improved.
 */
#include "taylor_bar.h" /**< Case setup for this example. */
#include "sphinxsys.h"

using namespace SPH;

int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.setGenerateRegressionData(false);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    /** create a body with corresponding material, particles and reaction model. */
    SolidBody column(sph_system, makeShared<Column>("Column"));
    column.defineAdaptationRatios(1.3, 1.0);
    column.defineBodyLevelSetShape(2.0)->writeLevelSet();
    column.defineMaterial<HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? column.generateParticles<BaseParticles, Reload>(column.getName())
        : column.generateParticles<BaseParticles, Lattice>();

    SolidBody wall(sph_system, makeShared<WallShape>("Wall"));
    wall.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall.generateParticles<BaseParticles, Lattice>();

    /** Define Observer. */
    ObserverBody my_observer(sph_system, "MyObserver");
    StdVec<Vecd> observation_location = {Vecd(0.0, 0.0, PW)};
    my_observer.generateParticles<ObserverParticles>(observation_location);

    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    SurfaceContactRelation column_wall_contact(column, {&wall});
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(sph_system);

    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_column_particles(column);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_column_to_vtp(column);
        /** Write the particle reload files. */

        ReloadParticleIO write_particle_reload_files(column);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(column_inner);
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
    SimpleDynamics<InitialCondition> initial_condition(column);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(column_inner);

    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf> stress_relaxation_first_half(column_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(column_inner);
    InteractionDynamics<DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(column, 0.2);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", my_observer_contact);

    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    corrected_configuration.exec();
    initial_condition.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % screen_output_interval == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";

                if (ite != 0 && ite % observation_sample_interval == 0)
                {
                    write_displacement.writeToFile(ite);
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
            physical_time += dt;
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

    if (sph_system.GenerateRegressionData())
    {
        write_displacement.generateDataBase(0.1);
    }
    else
    {
        write_displacement.testResult();
    }

    return 0;
}
