/**
 * @file dambreak_decomp.cpp
 * @brief 2D dambreak used to bring up the domain decomposition on the host path.
 * @details Derived from test_2d_dambreak_sycl. The observer, the regression tests
 *          and the restart output are removed; the total mechanical energy is
 *          recorded more often so that runs can be compared line by line.
 *          With SPHINXSYS_MULTI_SUBDOMAIN_HOST=OFF this file must reproduce the
 *          original case bit for bit.
 */
#include "sphinxsys.h"

#include <cstdlib>
#include <string>

using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Run time controls read from the environment, so that the same binary can be
//	compared in different configurations:
//	  SPHINXSYS_THREADS     TBB worker threads; 1 makes a run bit reproducible
//	  SPHINXSYS_SUBDOMAINS  number of host subdomains (decomposed build only)
//----------------------------------------------------------------------
int environmentInt(const char *name, int default_value)
{
    const char *value = std::getenv(name);
    return value != nullptr ? std::atoi(value) : default_value;
}
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
Real end_time = 2.0;                /**< Shorter than the reference case; enough for the front to cross the tank. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Artificial sound speed. */
//----------------------------------------------------------------------
//	Parameters for the shapes.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    const int number_of_threads = environmentInt("SPHINXSYS_THREADS", int(std::thread::hardware_concurrency()));
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    // The runner fixes the number of replicas every variable allocates, so it must be
    // initialized before any particle data exists.
    const int number_of_subdomains = environmentInt("SPHINXSYS_SUBDOMAINS", 1);
    // SPHINXSYS_RUNNER_THREADED=1 drives each subdomain from its own host thread, with
    // the barrier structure of the multi-device path; the default visits the subdomains
    // one after another on the calling thread, which is the configuration to debug in.
    execution::subdomain_runner.initialize(
        number_of_subdomains, environmentInt("SPHINXSYS_RUNNER_THREADED", 0) != 0
                                  ? execution::SubdomainRunner::Mode::Threaded
                                  : execution::SubdomainRunner::Mode::Sequential);
#endif
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref, number_of_threads);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with inital shape, materials and particles.
    //----------------------------------------------------------------------
    GeometricShapeBox initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    ComplexShape wall_complex_shape("WallBoundary");
    wall_complex_shape.add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
    wall_complex_shape.subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    SolidBody wall_boundary(sph_system, wall_complex_shape);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, wall_boundary);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    auto &host_methods = sph_solver.getHostMethodContainer();
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    //----------------------------------------------------------------------
    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &water_block_update_complex_relation = main_methods.addRelationDynamics(water_block_inner, water_wall_contact);
    auto &particle_sort = main_methods.addSortDynamics(water_block);

    Gravity gravity(Vecd(0.0, -gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(water_block, gravity);
    auto &wall_boundary_normal_direction = host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary); // run on CPU
    auto &water_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);

    auto &fluid_linear_correction_matrix =
        main_methods.addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(water_block_inner, 0.5)
            .addPostContactInteraction(water_wall_contact);
    auto &fluid_acoustic_step_1st_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep1stHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact);
    auto &fluid_acoustic_step_2nd_half =
        main_methods.addInteractionDynamics<
                        fluid_dynamics::AcousticStep2ndHalf, OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(water_wall_contact);
    auto &fluid_density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::CompressionSummation>(water_block_inner)
            .addPostContactInteraction(water_wall_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, WeaklyCompressibleFluid, FreeSurface>(water_block);

    auto &fluid_advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_ref);
    auto &fluid_acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<WeaklyCompressibleFluid>>(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_state_recorder.addToWrite<Real>(water_block, "Density");
    auto &record_water_mechanical_energy =
        main_methods.addIODynamics<ReducedQuantityRecording, TotalMechanicalEnergyCK>(water_block, gravity);
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    //----------------------------------------------------------------------
    //	Domain decomposition of the water body. Defined after every dynamics, so that
    //	all state variables are registered before the exchange set is fixed.
    //----------------------------------------------------------------------
    BaseParticles &water_particles = water_block.getBaseParticles();
    // The exchange set is the set of evolving variables: what a particle carries along
    // when it migrates. Two variables read at the neighbors by the interaction kernels
    // are not evolving by default and are added here, so that they can be refreshed
    // on the halo as well.
    water_particles.addEvolvingVariable<Real>("Pressure");
    water_particles.addEvolvingVariable<Matd>("LinearCorrectionMatrix");
    water_particles.addEvolvingVariable<Real>("Density"); // written to the vtp output

    const Real halo_width = water_block.getSPHAdaptation().getKernel()->CutOffRadius();
    SlabDecomposition decomposition(sph_system.getSystemDomainBounds(), halo_width, number_of_subdomains, 0);
    {   // The water column occupies only part of the tank, so an even split of the tank
        // would leave a subdomain empty. Balance the cut planes on the initial particle
        // positions before scattering.
        Vecd *position = water_particles.dvParticlePosition()->Data();
        const UnsignedInt total_particles = water_particles.TotalRealParticles();
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            StdVec<UnsignedInt> counts(number_of_subdomains, 0);
            for (UnsignedInt i = 0; i < total_particles; ++i)
                counts[decomposition.getSubdomainMap().subdomainOf(position[i])]++;
            if (!decomposition.rebalance(counts, 1.0))
                break;
        }
        std::cout << decomposition.describe();
    }
    SubdomainExchange<MainExecutionPolicy> exchange(decomposition, water_particles, water_particles.EvolvingVariables());
    UpdateHaloCK<MainExecutionPolicy> update_halo(exchange);
    MigrateParticlesCK<MainExecutionPolicy> migrate_particles(exchange);
    // Refreshed in the advection step, each right after the stage that writes it: the
    // volume is read at the neighbors by the kernel correction and by both acoustic
    // halves, the correction matrix by the first acoustic half. Both change once per
    // advection step only, which is why they are not interact variables of the acoustic
    // steps (those would be re-sent every acoustic step).
    SyncHaloStateCK<MainExecutionPolicy> sync_volume(exchange, water_particles);
    sync_volume.addVariable<Real>("VolumetricMeasure");
    SyncHaloStateCK<MainExecutionPolicy> sync_correction(exchange, water_particles);
    sync_correction.addVariable<Matd>("LinearCorrectionMatrix");
    // The pressure (first half) and the velocity (second half) are refreshed by the
    // acoustic steps themselves: each interaction algorithm refreshes the halo of its
    // interact variables right before its interaction step, through the HaloRefresher
    // installed on the particles by the exchange.

    auto check_consistency = [&](size_t step)
    {
        std::string report = exchange.checkConsistency();
        if (!report.empty())
        {
            std::cout << "Decomposition inconsistent at advection step " << step << ": " << report << std::endl;
            exit(1);
        }
    };
    // Output only: which subdomain owns each particle. gatherToHost() lays the owned
    // particles out subdomain by subdomain in the host arrays, so the tag can be filled
    // on the host right before writing, without any exchange.
    DiscreteVariable<int> *dv_subdomain_id = water_particles.registerStateVariable<int>("SubdomainID");
    body_state_recorder.addToWrite<int>(water_block, "SubdomainID");
    auto tag_subdomains_for_output = [&]()
    {
        int *subdomain_id = dv_subdomain_id->Data();
        UnsignedInt offset = 0;
        StdVec<UnsignedInt> owned = exchange.OwnedParticlesPerSubdomain();
        for (int s = 0; s < number_of_subdomains; ++s)
        {
            for (UnsignedInt i = offset; i < offset + owned[s]; ++i)
                subdomain_id[i] = s;
            offset += owned[s];
        }
    };
    auto report_subdomains = [&]()
    {
        StdVec<UnsignedInt> owned = exchange.OwnedParticlesPerSubdomain();
        std::cout << "    owned particles per subdomain:";
        for (UnsignedInt count : owned)
            std::cout << " " << count;
        std::cout << "  total " << exchange.TotalOwnedParticles()
                  << "  halo load factor " << exchange.HaloLoadFactor() << "\n";
    };
#endif
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    //----------------------------------------------------------------------
    //	Setup for advection-step based time-stepping control
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t advection_steps = 1;
    int screening_interval = 100;
    int observation_interval = 20;
    // SPHINXSYS_VTP_INTERVAL (seconds, default 0.1) sets the state recording interval; a
    // very small value records every advection step, which is what a step-by-step
    // comparison between two configurations needs.
    const char *vtp_interval = std::getenv("SPHINXSYS_VTP_INTERVAL");
    auto &state_recording = time_stepper.addTriggerByInterval(
        vtp_interval != nullptr ? std::atof(vtp_interval) : 0.1);
    //----------------------------------------------------------------------
    //	Prepare for the time integration loop.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec(); // run particle dynamics with host kernels first
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    exchange.scatterFromHost(); // before any dynamics touches the water replicas
#endif
    constant_gravity.exec();

#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    update_halo.exec(); // publishes n_local, which the cell linked list is built over
#endif
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();

    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    sync_volume.exec();
#endif
    fluid_linear_correction_matrix.exec();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    sync_correction.exec();
    check_consistency(advection_steps);
    report_subdomains();
#endif
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    exchange.gatherToHost();
    tag_subdomains_for_output();
#endif
    body_state_recorder.writeToFile();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
    exchange.finishHostAccess();
#endif
    record_water_mechanical_energy.writeToFile(advection_steps);
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
    // SPHINXSYS_DEBUG_STEPS = N prints the time step of the first N acoustic steps and
    // of every advection step among them, to compare two configurations step by step.
    const int debug_steps = environmentInt("SPHINXSYS_DEBUG_STEPS", 0);
    size_t acoustic_steps = 0;
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime(end_time))
    {
        //----------------------------------------------------------------------
        //	the fastest and most frequent acostic time stepping.
        //----------------------------------------------------------------------
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        if (debug_steps > 0 && acoustic_steps < debug_steps)
        {
            std::cout << std::setprecision(17) << "acoustic step " << acoustic_steps
                      << " dt = " << acoustic_dt << " time = " << time_stepper.getPhysicalTime() << "\n";
        }
        acoustic_steps++;
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;
        //----------------------------------------------------------------------
        //	the following are slower and less frequent time stepping.
        //----------------------------------------------------------------------
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            if (debug_steps > 0 && acoustic_steps < debug_steps)
            {
                std::cout << std::setprecision(17) << "advection step " << advection_steps
                          << " advection_dt = " << advection_step.getInterval() << "\n";
            }
            water_update_particle_position.exec();

            /** Output body state during the simulation according output_interval. */
            time_instance = TickCount::now();
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << advection_steps
                          << "	Time = " << time_stepper.getPhysicalTime() << "	"
                          << "	advection_dt = " << advection_step.getInterval()
                          << "	acoustic_dt = " << time_stepper.getGlobalTimeStepSize() << "\n";
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
                report_subdomains();
#endif
            }

            if (advection_steps % observation_interval == 0)
            {
                record_water_mechanical_energy.writeToFile(advection_steps);
            }

            if (state_recording())
            {
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
                exchange.gatherToHost();
                tag_subdomains_for_output();
#endif
                body_state_recorder.writeToFile();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
                exchange.finishHostAccess();
#endif
            }
            interval_output += TickCount::now() - time_instance;

            /** Particle sort, update cell linked list and configuration. */
            time_instance = TickCount::now();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
            migrate_particles.exec(); // ownership follows the new positions
#endif
            if (advection_steps % 100)
            {
                particle_sort.exec();
            }
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
            update_halo.exec(); // new halo plan and full refresh, before the cell linked list
#endif
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
            sync_volume.exec();
#endif
            fluid_linear_correction_matrix.exec();
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
            sync_correction.exec();
            check_consistency(advection_steps);
#endif
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

    return 0;
};
