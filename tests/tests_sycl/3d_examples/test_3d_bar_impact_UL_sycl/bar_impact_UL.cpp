/**
 * @file bar_impact_UL.cpp
 * @brief This is the case setup for a plastic bar impacting on rigid wall
 * and feedbacked by contact force using updated Lagrangian SPH.
 * @author Shuaihao Zhang, Dong Wu and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real total_physical_time = 6.0e-5; /**< TOTAL SIMULATION TIME*/
Real PL = 0.00391;                 /**< X-direction domain. */
Real PW = 0.02346;                 /**< Z-direction domain. */
Real particle_spacing_ref = PL / 12.0;
Real column_radius = PL;
Vecd translation_column(0.0, 0.0, 0.5 * PW + particle_spacing_ref);
Real SL = particle_spacing_ref * 4.0;
Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
Vecd translation_holder(0.0, 0.0, -0.5 * SL);
int resolution(20);
StdVec<Vecd> observation_location = {Vecd(0.0, 0.0, PW)};
Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 2700.0; /**< Reference density. */
Real poisson = 0.3;   /**< Poisson ratio. */
Real Youngs_modulus = 78.2e9;
Real yield_stress = 0.29e9;
Real vel_0 = 373.0;
Real U_max = vel_0;
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
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
        if (!sph_system.ReloadParticles())
        {
            return 0;
        }
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
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(column, "Velocity", Vec3d(0, 0, -vel_0)).exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    ParticleDynamicsGroup update_column_configuration;
    update_column_configuration.add(&main_methods.addCellLinkedListDynamics(column));
    update_column_configuration.add(&main_methods.addRelationDynamics(column_inner, column_wall_contact));
    auto &wall_boundary_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &column_observer_contact_relation = main_methods.addRelationDynamics(column_observer_contact);

    auto &column_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(column);
    auto &column_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(column);
    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);

    ParticleDynamicsGroup column_shear_force;
    column_shear_force.add(&main_methods.addInteractionDynamics<LinearGradient, Vecd>(column_inner, "Velocity"));
    column_shear_force.add(&main_methods.addInteractionDynamicsOneLevel<continuum_dynamics::ShearIntegration, J2Plasticity>(column_inner));

    auto &column_acoustic_step_1st_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep1stHalf, DissipativeRiemannSolverCK, NoKernelCorrectionCK>(column_inner);
    auto &column_acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf, DissipativeRiemannSolverCK, NoKernelCorrectionCK>(column_inner);

    auto &column_wall_contact_factor = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(column_wall_contact);
    auto &column_wall_contact_force = main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::RepulsionForceCK, Wall>(column_wall_contact);

    auto &column_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(column, U_max, 0.2);
    auto &column_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(column, 0.4);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(column, "Pressure");
    body_state_recorder.addToWrite<Real>(column, "Density");
    auto &record_column_mechanical_energy = main_methods.addReduceRegression<
        RegressionTestDynamicTimeWarping, TotalKineticEnergyCK>(column);
    auto &column_observer_position = main_methods.addObserveRecorder<Vecd>("Position", column_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(column_advection_time_step.exec());
    size_t advection_steps = sph_system.RestartStep() + 1;
    int screening_interval = 100;
    int observation_interval = screening_interval * 2;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 60.0);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_boundary_cell_linked_list.exec();
    update_column_configuration.exec();

    column_advection_step_setup.exec();
    column_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    record_column_mechanical_energy.writeToFile(advection_steps);
    column_observer_position.writeToFile(advection_steps);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;
    //----------------------------------------------------------------------
    //	Single time stepping loop is used for multi-time stepping.
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(column_acoustic_time_step);
        column_shear_force.exec(acoustic_dt);
        column_wall_contact_force.exec();
        column_acoustic_step_1st_half.exec(acoustic_dt);
        column_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(column_advection_time_step))
        {
            advection_steps++;
            column_update_particle_position.exec();

            /** Output body state during the simulation according output_interval. */
            time_instance = TickCount::now();
            /** screen output, write body observables and restart files  */
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                record_column_mechanical_energy.writeToFile(advection_steps);
                column_observer_contact_relation.exec();
                column_observer_position.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
            update_column_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            column_advection_step_setup.exec();
            column_wall_contact_factor.exec();
            column_linear_correction_matrix.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_advection_step ="
              << interval_advection_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        record_column_mechanical_energy.generateDataBase(5.0e-2);
    }
    else
    {
        record_column_mechanical_energy.testResult();
    }

    return 0;
}
