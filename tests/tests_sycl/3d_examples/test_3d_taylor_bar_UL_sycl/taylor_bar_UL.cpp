/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar using updated Lagragian SPH.
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu
 */
#include "taylor_bar_UL.h"
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(true);
    sph_system.handleCommandlineOptions(ac, av);

    auto &column_shape = sph_system.addShape<TriangleMeshShapeCylinder>(
        Vec3d(0, 0, 1.0), column_radius, 0.5 * PW, resolution, translation_column, "Column");
    auto &wall_shape = sph_system.addShape<TriangleMeshShapeBrick>(
        halfsize_holder, resolution, translation_holder, "Wall");

    auto &column = sph_system.addBody<RealBody>(column_shape);
    auto &wall_boundary = sph_system.addBody<SolidBody>(wall_shape);
    if (sph_system.RunParticleRelaxation())
    {
        LevelSetShape *level_set_shape = column.defineBodyLevelSetShape(par_ck, 2.0)->writeLevelSet();
        column.generateParticles<BaseParticles, Lattice>();
        wall_boundary.generateParticles<BaseParticles, Lattice>();
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
            main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(column_inner)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(column, *level_set_shape);
        auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(column);
        auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(column);
        auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
        //----------------------------------------------------------------------
        //	Run on CPU after relaxation finished and output results.
        //----------------------------------------------------------------------
        auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(column);
        auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&column, &wall_boundary});
        write_particle_reload_files.addToReload<Vecd>(wall_boundary, "NormalDirection");
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
    column.defineMaterial<J2Plasticity>(rho0_s, c0, Youngs_modulus, poisson, yield_stress);
    column.generateParticles<BaseParticles, Reload>(column.getName());

    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");

    auto &column_observer = sph_system.addBody<ObserverBody>("ColumnObserver");
    column_observer.generateParticles<ObserverParticles>(observation_location);

    /**body relation topology */
    auto &column_inner = sph_system.addInnerRelation(column);
    auto &column_wall_contact = sph_system.addContactRelation(column, wall_boundary);
    auto &column_observer_contact = sph_system.addContactRelation(column_observer, column);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(column, "Vecolity", Vec3d(0, 0, -vel_0)).exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup update_column_configuration;
    update_column_configuration.add(&main_methods.addCellLinkedListDynamics(column));
    update_column_configuration.add(&main_methods.addRelationDynamics(column_inner, column_wall_contact));
    auto &wall_boundary_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);

    auto &soil_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(column);
    auto &soil_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(column);
    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);

    auto &soil_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel< // to check why not use Riemann solver
            fluid_dynamics::AcousticStep1stHalf, NoRiemannSolverCK, NoKernelCorrectionCK>(column_inner);
    auto &soil_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf, DissipativeRiemannSolverCK, NoKernelCorrectionCK>(column_inner);

    InteractionWithUpdate<continuum_dynamics::ShearStressRelaxationHourglassControl1stHalfJ2Plasticity> column_shear_stress(column_inner);
    InteractionDynamics<continuum_dynamics::ShearStressRelaxationHourglassControl2ndHalf> column_shear_acceleration(column_inner);
    auto &column_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(column, U_max, 0.2);
    auto &column_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(column, 0.4);
    InteractionDynamics<DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(column, "Pressure");
    write_states.addToWrite<Real>(column, "Density");
    SimpleDynamics<continuum_dynamics::VonMisesStress> column_von_mises_stress(column);
    write_states.addToWrite<Real>(column, "VonMisesStress");
    ObservedQuantityRecording<Vecd> write_displacement("Position", column_observer_contact);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_kinetic_energy(column);
    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
    write_kinetic_energy.writeToFile(0);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            Real relaxation_time = 0.0;
            Real advection_dt = advection_time_step.exec();
            column_volume_update.exec();
            while (relaxation_time < advection_dt)
            {
                Real acoustic_dt = acoustic_time_step.exec();
                if (ite % screen_output_interval == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	advection_dt: "
                              << advection_dt << "	acoustic_dt: "
                              << acoustic_dt << "\n";
                }
                column_wall_contact_force.exec(acoustic_dt);

                column_pressure_relaxation.exec(acoustic_dt);
                column_shear_stress.exec(acoustic_dt);
                column_shear_acceleration.exec(acoustic_dt);
                column_density_relaxation.exec(acoustic_dt);

                ite++;
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                physical_time += acoustic_dt;
            }
            column.updateCellLinkedList();
            column_inner.updateConfiguration();
            column_wall_contact.updateConfiguration();
            corrected_configuration.exec();
        }
        TickCount t2 = TickCount::now();
        column_von_mises_stress.exec();
        write_states.writeToFile(ite);
        write_displacement.writeToFile(ite);
        write_kinetic_energy.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall_boundary time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_kinetic_energy.generateDataBase(5.0e-2);
    }
    else
    {
        write_kinetic_energy.testResult();
    }

    return 0;
}
