/**
 * @file taylor_bar_sycl.cpp
 * @brief This test offloads the levelset computation to the GPU device.
 *        All other components remain unchanged from the original Taylor bar test.
 *        The regression test is also adapted to validate the results.
 * @author Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 * //TODO: Seems that the wall contact force should be improved.
 */
#include "taylor_bar_sycl.h" /**< Case setup for this example. */
// #include "sphinxsys.h"

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

    SolidBody wall(sph_system, makeShared<WallShape>("Wall"));

    if (sph_system.RunParticleRelaxation())
    {
        LevelSetShape *level_set_shape = column.defineBodyLevelSetShape(par_ck, 2.0)->writeLevelSet();
        column.generateParticles<BaseParticles, Lattice>();
        wall.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_body_surface(column);
        Inner<> column_inner(column);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(sph_system);
        auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
        auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

        auto &input_body_cell_linked_list = main_methods.addCellLinkedListDynamics(column);
        auto &input_body_update_inner_relation = main_methods.addRelationDynamics(column_inner);
        auto &random_input_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(column);
        auto &relaxation_residual =
            main_methods.addInteractionDynamics<RelaxationResidualCK, NoKernelCorrectionCK>(column_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(column, *level_set_shape);
        auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(column);
        auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(column);
        auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
        //----------------------------------------------------------------------
        //	Run on CPU after relaxation finished and output results.
        //----------------------------------------------------------------------
        auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(column);
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&column, &wall});
        write_particle_reload_files.addToReload<Vecd>(wall, "NormalDirection");
        //----------------------------------------------------------------------
        //	Prepare the simulation with cell linked list, configuration
        //	and case specified initial condition if necessary.
        //----------------------------------------------------------------------
        random_input_body_particles.exec();

        //----------------------------------------------------------------------
        //	First output before the simulation.
        //----------------------------------------------------------------------
        body_state_recorder.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            input_body_cell_linked_list.exec();
            input_body_update_inner_relation.exec();

            relaxation_residual.exec();
            Real relaxation_step = relaxation_scaling.exec();
            update_particle_position.exec(relaxation_step);
            level_set_bounding.exec();

            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;
        /** Output results. */
        wall_boundary_normal_direction.exec();
        write_particle_reload_files.writeToFile();
        return 0;
    }

    column.defineMaterial<HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
    column.generateParticles<BaseParticles, Reload>(column.getName());

    wall.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall.generateParticles<BaseParticles, Reload>(wall.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");

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
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<::InitialCondition> initial_condition(column);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(column_inner);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
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
