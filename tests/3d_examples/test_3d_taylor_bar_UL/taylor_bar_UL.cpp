/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar using updated Lagragian SPH.
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu
 */
#include "taylor_bar_UL.h"
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    system.setRunParticleRelaxation(false);
    system.setReloadParticles(true);

    RealBody column(system, makeShared<Column>("Column"));
    column.defineBodyLevelSetShape()->writeLevelSet(system);
    column.defineMaterial<J2Plasticity>(rho0_s, c0, Youngs_modulus, poisson, yield_stress);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? column.generateParticles<BaseParticles, Reload>(column.getName())
        : column.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    
    ObserverBody my_observer(system, "MyObserver");
    StdVec<Vecd> observation_location = {Vecd(0.0, 0.0, PW)};
    my_observer.generateParticles<ObserverParticles>(observation_location);

    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, {&column});
    SurfaceContactRelation column_wall_contact(column, {&wall_boundary});
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(system);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
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
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(column_inner);
    Dynamics1Level<continuum_dynamics::ShearStressRelaxationHourglassControlJ2Plasticity> column_shear_stress_relaxation(column_inner);
    Dynamics1Level<continuum_dynamics::Integration1stHalf> column_pressure_relaxation(column_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerDissipativeRiemann> column_density_relaxation(column_inner);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> beam_volume_update(column);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> advection_time_step(column, U_max, 0.2);
    ReduceDynamics<continuum_dynamics::AcousticTimeStepSize> acoustic_time_step(column, 0.4);
    InteractionDynamics<DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    write_states.addToWrite<Real>(column, "VonMisesStress");
    write_states.addToWrite<Real>(column, "VonMisesStrain");
    write_states.addToWrite<Real>(column, "Pressure");
    write_states.addToWrite<Real>(column, "Density");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_displacement("Position", my_observer_contact);
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
    Real end_time = 6.0e-5;
    int screen_output_interval = 100;
    Real output_period = end_time / 60;
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            Real relaxation_time = 0.0;
            Real advection_dt = advection_time_step.exec();
            beam_volume_update.exec();
            while (relaxation_time < advection_dt)
            {
                Real acoustic_dt = acoustic_time_step.exec();
                if (ite % screen_output_interval == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	advection_dt: "
                              << advection_dt << "	acoustic_dt: "
                              << acoustic_dt << "\n";
                }
                column_wall_contact_force.exec(acoustic_dt);

                column_pressure_relaxation.exec(acoustic_dt);
                column_shear_stress_relaxation.exec(acoustic_dt);
                column_density_relaxation.exec(acoustic_dt);

                ite++;
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                GlobalStaticVariables::physical_time_ += acoustic_dt;
            }
            column.updateCellLinkedList();
            column_inner.updateConfiguration();
            column_wall_contact.updateConfiguration();
            corrected_configuration.exec();
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile(ite);
        write_displacement.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall_boundary time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (system.GenerateRegressionData())
    {
        write_displacement.generateDataBase(0.1);
    }
    else
    {
        write_displacement.testResult();
    }

    return 0;
}
