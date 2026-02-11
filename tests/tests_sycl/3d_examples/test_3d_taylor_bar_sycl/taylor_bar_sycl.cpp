/**
 * @file taylor_bar_sycl.cpp
 * @brief This test offloads the levelset computation to the GPU device.
 *        All other components remain unchanged from the original Taylor bar test.
 *        The regression test is also adapted to validate the results.
 * @author Xiaojing Tang, Dong Wu and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 */
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real total_physical_time = 1.0e-4; /**< TOTAL SIMULATION TIME*/
Real PL = 0.0032;                  /**< X-direction domain. */
Real PW = 0.0324;                  /**< Z-direction domain. */
Real particle_spacing_ref = PL / 5.0;
/** YOU can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real SL = particle_spacing_ref * 4.0; /**< Length of the holder is one layer particle. */
Real column_radius = PL;
Vec3d translation_column(0.0, 0.0, 0.6 * PW);
Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
Vecd translation_holder(0.0, 0.0, -0.5 * SL);
Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
StdVec<Vecd> observation_location{Vecd(0.0, 0.0, PW)};
Real impact_speed = 227.0;
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; /**< Reference density. */
Real poisson = 0.35;  /**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Setup geometry first.
    //----------------------------------------------------------------------
    auto &column_shape = sph_system.addShape<TriangleMeshShapeCylinder>(
        Vec3d(0, 0, 1.0), column_radius, 0.5 * PW, resolution, translation_column, "Column");
    auto &wall_shape = sph_system.addShape<TriangleMeshShapeBrick>(
        halfsize_holder, resolution, translation_holder, "Wall");
    //----------------------------------------------------------------------
    //	Run particle relaxation first if needed.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        // setup a sub-system for particle relaxation and delete it after particle relaxation.
        RelaxationSystem relaxation_system(system_domain_bounds, particle_spacing_ref);
        auto &column = relaxation_system.addBody<RealBody>(column_shape);
        auto &wall = relaxation_system.addBody<SolidBody>(wall_shape);

        LevelSetShape *level_set_shape = column.defineBodyLevelSetShape(par_ck, 2.0)->writeLevelSet();
        column.generateParticles<BaseParticles, Lattice>();
        wall.generateParticles<BaseParticles, Lattice>();
        NearShapeSurface near_body_surface(column);
        Inner<> column_inner(column);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SPHSolver sph_solver(relaxation_system);
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
        auto &wall_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall);
        //----------------------------------------------------------------------
        //	Define simple file input and outputs functions.
        //----------------------------------------------------------------------
        auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(relaxation_system);
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
        wall_normal_direction.exec();
        write_particle_reload_files.writeToFile();

        if (!sph_system.ReloadParticles())
        {
            return 0;
        }
        else
        {
            std::cout << "To reload particles and start the main simulation." << std::endl;
        }
    }
    //----------------------------------------------------------------------
    //	Simulation setup continues to define bodies.
    //----------------------------------------------------------------------
    auto &column = sph_system.addBody<RealBody>(column_shape);
    column.defineMaterial<HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
    column.generateParticles<BaseParticles, Reload>(column.getName());

    auto &wall = sph_system.addBody<RealBody>(wall_shape);
    wall.defineMaterial<Solid>();
    wall.generateParticles<BaseParticles, Reload>(wall.getName())
        ->reloadExtraVariable<Vecd>("NormalDirection");

    auto &my_observer = sph_system.addBody<ObserverBody>("MyObserver");
    my_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	All body relations.
    //----------------------------------------------------------------------
    auto &column_inner = sph_system.addInnerRelation(column, ConfigType::Lagrangian);
    auto &column_wall_contact = sph_system.addContactRelation(column, wall);
    auto &my_observer_contact = sph_system.addContactRelation(my_observer, column, ConfigType::Lagrangian);
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
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(column, "Velocity", Vec3d(0, 0, -impact_speed)).exec();

    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall);
    ParticleDynamicsGroup update_conact_configuration;
    update_conact_configuration.add(&main_methods.addCellLinkedListDynamics(column));
    update_conact_configuration.add(&main_methods.addRelationDynamics(column_wall_contact));
    auto &column_inner_configuration = main_methods.addRelationDynamics(column_inner);
    auto &my_observer_contact_relation = main_methods.addRelationDynamics(my_observer_contact);

    auto &column_linear_correction_matrix = main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(column_inner);
    auto &column_wall_contact_factor = main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(column_wall_contact);
    auto &column_wall_contact_force = main_methods.addInteractionDynamicsWithUpdate<solid_dynamics::RepulsionForceCK, Wall>(column_wall_contact);

    auto &column_acoustic_step_1st_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration1stHalf, HardeningPlasticSolid>(column_inner);
    auto &column_acoustic_step_2nd_half = main_methods.addInteractionDynamicsOneLevel<
        solid_dynamics::StructureIntegration2ndHalf>(column_inner);

    auto &column_acoustic_time_step = main_methods.addReduceDynamics<solid_dynamics::AcousticTimeStepCK>(column, 0.2);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_displacement = main_methods.addObserveRegression<
        RegressionTestDynamicTimeWarping, Vecd>("Position", my_observer_contact);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    size_t acoustic_steps = 1;
    int screening_interval = 100;
    int observation_interval = screening_interval * 2;
    auto &state_recording = time_stepper.addTriggerByInterval(total_physical_time / 60.0);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_cell_linked_list.exec();
    update_conact_configuration.exec();
    column_inner_configuration.exec();
    my_observer_contact_relation.exec();
    column_linear_correction_matrix.exec();
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    write_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Statistics for the computing time information
    //----------------------------------------------------------------------
    TimeInterval interval_output;
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
        column_wall_contact_factor.exec();
        column_wall_contact_force.exec();
        column_acoustic_step_1st_half.exec(acoustic_dt);
        column_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        /** Output body state during the simulation according output_interval. */
        time_instance = TickCount::now();
        /** screen output, write body observables and restart files  */
        if (acoustic_steps % screening_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << acoustic_steps
                      << "	Time = " << time_stepper.getPhysicalTime() << "	"
                      << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (acoustic_steps % observation_interval == 0)
        {
            write_displacement.writeToFile(acoustic_steps);
        }

        if (state_recording())
        {
            body_state_recorder.writeToFile();
        }
        interval_output += TickCount::now() - time_instance;

        /** Particle sort, update cell linked list and configuration. */
        time_instance = TickCount::now();
        update_conact_configuration.exec();
        interval_updating_configuration += TickCount::now() - time_instance;

        acoustic_steps++;
    }
    //----------------------------------------------------------------------
    // Summary for wall time used for the simulation.
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_acoustic_step = "
              << interval_acoustic_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

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
