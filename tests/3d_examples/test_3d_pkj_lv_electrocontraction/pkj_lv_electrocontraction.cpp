/**
 * @file 	pkj_lv_electrocontraction.cpp
 * @brief 	This is the case studying the electromechanics on a left ventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */
#include "pkj_lv_electrocontraction.h"
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/**
 * The main program.
 */
int main(int ac, char *av[])
{
    /**
     * Build up context -- a SPHSystem.
     */
    SPHSystem sph_system(system_domain_bounds, dp_0);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for reload initially relaxed particles. */
    sph_system.setReloadParticles(false);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    Heart triangle_mesh_heart_model("HeartModel");
    SharedPtr<LevelSetShape> level_set_heart_model =
        makeShared<LevelSetShape>(par_ck, sph_system, SPHAdaptation(dp_0), triangle_mesh_heart_model);
    level_set_heart_model->correctLevelSetSign()->writeLevelSet();
    std::cout << "I am here 0! " << "\n";
    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    /** check whether run particle relaxation for body fitted particle distribution. */
    if (sph_system.RunParticleRelaxation())
    {
        SolidBody heart_model(sph_system, level_set_heart_model);
        std::cout << "I am here 1! " << "\n";
        heart_model.defineClosure<LocallyOrthotropicMuscle, IsotropicDiffusion>(
            ConstructArgs(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0),
            ConstructArgs(diffusion_species_name, diffusion_coeff));
        heart_model.generateParticles<BaseParticles, Lattice>();
        /** topology */
        InnerRelation heart_model_inner(heart_model);
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_particles(heart_model);
        RelaxationStepInner relaxation_step_inner(heart_model_inner);
        BodyStatesRecordingToVtp write_heart_model_state_to_vtp({heart_model});
        ReloadParticleIO write_particle_reload_files(heart_model);
        write_particle_reload_files.addToReload<Vecd>(heart_model, "Fiber");
        write_particle_reload_files.addToReload<Vecd>(heart_model, "Sheet");
        //----------------------------------------------------------------------
        //	Physics relaxation starts here.
        //----------------------------------------------------------------------
        random_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_heart_model_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        // From here the time stepping begins.
        //----------------------------------------------------------------------
        std::cout << "I am here 2! " << "\n";
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            relaxation_step_inner.exec();
            ite++;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_heart_model_state_to_vtp.writeToFile(ite);
            }
        }
        std::cout << "I am here 3! " << "\n";
        //----------------------------------------------------------------------
        //	Diffusion process to initialize fiber direction
        //----------------------------------------------------------------------
        GetDiffusionTimeStepSize get_time_step_size(heart_model);
        FiberDirectionDiffusionRelaxation diffusion_relaxation(heart_model_inner);
        SimpleDynamics<ComputeFiberAndSheetDirections> compute_fiber_sheet(heart_model, diffusion_species_name);
        BodySurface surface_part(heart_model);
        SimpleDynamics<DiffusionBCs> impose_diffusion_bc(surface_part, diffusion_species_name);
        impose_diffusion_bc.exec();
        write_heart_model_state_to_vtp.addToWrite<Real>(heart_model, diffusion_species_name);
        write_heart_model_state_to_vtp.writeToFile(ite);

        int diffusion_step = 100;
        Real dt = get_time_step_size.exec();
        while (ite <= diffusion_step + relax_step)
        {
            diffusion_relaxation.exec(dt);
            impose_diffusion_bc.exec();
            if (ite % 10 == 0)
            {
                std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
                write_heart_model_state_to_vtp.writeToFile(ite);
            }
            ite++;
        }
        compute_fiber_sheet.exec();
        ite++;
        write_heart_model_state_to_vtp.writeToFile(ite);
        compute_fiber_sheet.exec();
        write_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	SPH simulation section
    //----------------------------------------------------------------------
    /** create a SPH body, material and particles */
    SolidBody physiology_heart(sph_system, level_set_heart_model, "PhysiologyHeart");
    AlievPanfilowModel aliev_panfilow_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
    physiology_heart.defineClosure<Solid, MonoFieldElectroPhysiology<LocalDirectionalDiffusion>>(
        Solid(), ConstructArgs(&aliev_panfilow_model, ConstructArgs(diffusion_coeff, bias_coeff, fiber_direction)));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? physiology_heart.generateParticles<BaseParticles, Reload>("HeartModel")
        : physiology_heart.generateParticles<BaseParticles, Lattice>();

    /** create a SPH body, material and particles */
    SolidBody mechanics_heart(sph_system, level_set_heart_model, "MechanicalHeart");
    mechanics_heart.defineMaterial<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? mechanics_heart.generateParticles<BaseParticles, Reload>("HeartModel")
        : mechanics_heart.generateParticles<BaseParticles, Lattice>();

    /** Creat a Purkinje network for fast diffusion, material and particles */
    TreeBody pkj_body(sph_system, level_set_heart_model, "Purkinje");
    pkj_body.defineClosure<Solid, MonoFieldElectroPhysiology<IsotropicDiffusion>>(
        Solid(), ConstructArgs(&aliev_panfilow_model, ConstructArgs(diffusion_coeff * acceleration_factor)));
    pkj_body.generateParticles<BaseParticles, NetworkWithExtraCheck>(starting_point, second_point, 50, 1.0);
    TreeTerminates pkj_leaves(pkj_body);
    //----------------------------------------------------------------------
    //	SPH Observation section
    //----------------------------------------------------------------------
    ObserverBody voltage_observer(sph_system, "VoltageObserver");
    voltage_observer.generateParticles<ObserverParticles>(createObservationPoints());

    ObserverBody myocardium_observer(sph_system, "MyocardiumObserver");
    myocardium_observer.generateParticles<ObserverParticles>(createObservationPoints());

    /** topology */
    InnerRelation physiology_heart_inner(physiology_heart);
    InnerRelation mechanics_heart_inner(mechanics_heart);
    ContactRelation physiology_heart_contact(physiology_heart, {&mechanics_heart});
    ContactRelation mechanics_heart_contact(mechanics_heart, {&physiology_heart});
    ContactRelation voltage_observer_contact(voltage_observer, {&physiology_heart});
    ContactRelation myocardium_observer_contact(myocardium_observer, {&mechanics_heart});
    // ComplexRelation physiology_heart_complex(physiology_heart, {&pkj_leaves});
    ContactRelationToBodyPart physiology_heart_contact_with_pkj_leaves(physiology_heart, {&pkj_leaves});
    TreeInnerRelation pkj_inner(pkj_body);

    /** Corrected configuration. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_excitation(physiology_heart_inner);
    /** Time step size calculation. */
    GetDiffusionTimeStepSize get_myocardium_physiology_time_step(physiology_heart);
    /** Diffusion process for diffusion body. */
    electro_physiology::ElectroPhysiologyDiffusionRelaxationComplex<LocalDirectionalDiffusion, Dirichlet> myocardium_diffusion_relaxation(
        physiology_heart_inner, physiology_heart_contact_with_pkj_leaves);
    /** Solvers for ODE system */
    electro_physiology::ElectroPhysiologyReactionRelaxationForward myocardium_reaction_relaxation_forward(physiology_heart, aliev_panfilow_model);
    electro_physiology::ElectroPhysiologyReactionRelaxationBackward myocardium_reaction_relaxation_backward(physiology_heart, aliev_panfilow_model);
    /** Physiology for PKJ*/
    /** Time step size calculation. */
    GetDiffusionTimeStepSize get_pkj_physiology_time_step(pkj_body);
    electro_physiology::ElectroPhysiologyDiffusionNetworkRK2 pkj_diffusion_relaxation(pkj_inner);
    /** Solvers for ODE system */
    electro_physiology::ElectroPhysiologyReactionRelaxationForward pkj_reaction_relaxation_forward(pkj_body, aliev_panfilow_model);
    electro_physiology::ElectroPhysiologyReactionRelaxationBackward pkj_reaction_relaxation_backward(pkj_body, aliev_panfilow_model);
    /**Apply the ion stimulus.*/
    SimpleDynamics<ApplyStimulusCurrentToMyocardium> apply_stimulus_myocardium(physiology_heart);
    SimpleDynamics<ApplyStimulusCurrentToPKJ> apply_stimulus_pkj(pkj_body);
    /** Active mechanics. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_contraction(mechanics_heart_inner);
    /** Observer Dynamics */
    InteractionDynamics<CorrectInterpolationKernelWeights>
        correct_kernel_weights_for_interpolation(mechanics_heart_contact);
    /** Interpolate the active contract stress from electrophysiology body. */
    InteractionDynamics<InterpolatingAQuantity<Real>>
        active_stress_interpolation(mechanics_heart_contact, "ActiveContractionStress", "ActiveContractionStress");
    /** Interpolate the particle position in physiology_heart  from mechanics_heart. */
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_particle_position(physiology_heart_contact, "Position", "Position");
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(mechanics_heart_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(mechanics_heart_inner);

    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> get_mechanics_time_step(mechanics_heart);
    /** Constrain region of the inserted body. */
    MuscleBaseShapeParameters muscle_base_parameters;
    BodyRegionByParticle muscle_base(mechanics_heart, makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(muscle_base);

    /**IO for observer.*/
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(physiology_heart, "Voltage");
    write_states.addToWrite<Real>(physiology_heart, "GateVariable");
    write_states.addToWrite<Real>(physiology_heart, "ActiveContractionStress");
    write_states.addToWrite<Real>(mechanics_heart, "ActiveContractionStress");
    ObservedQuantityRecording<Real> write_voltage("Voltage", voltage_observer_contact);
    ObservedQuantityRecording<Vecd> write_displacement("Position", myocardium_observer_contact);
    /**
     * Pre-simulation.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correct_configuration_excitation.exec();
    correct_configuration_contraction.exec();
    correct_kernel_weights_for_interpolation.exec();
    /**Output global basic parameters. */
    write_states.writeToFile(0);
    write_voltage.writeToFile(0);
    write_displacement.writeToFile(0);
    write_states.writeToFile(0);
    /**
     * main loop.
     */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int screen_output_interval = 10;
    int ite = 0;
    int reaction_step = 2;
    Real end_time = 80;
    Real Ouput_T = end_time / 200.0;
    Real Observer_time = 0.01 * Ouput_T;
    Real dt_myocardium = 0.0;
    Real dt_pkj = 0.0;
    Real dt_muscle = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    std::cout << "Main Loop Starts Here : "
              << "\n";
    /** Main loop starts here. */
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Ouput_T)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observer_time)
            {
                if (ite % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
                              << physical_time
                              << "	dt_pkj = " << dt_pkj
                              << "	dt_myocardium = " << dt_myocardium
                              << "	dt_muscle = " << dt_muscle << "\n";
                }
                /** Apply stimulus excitation. */
                // if( 0 <= physical_time
                // 	&&  physical_time <= 0.5)
                // {
                // 	apply_stimulus_myocardium.exec(dt_myocardium);
                // }

                Real dt_pkj_sum = 0.0;
                while (dt_pkj_sum < dt_myocardium)
                {
                    /**
                     * When network generates particles, the final particle spacing, which is after particle projected in to
                     * complex geometry, may small than the reference one, therefore, a smaller time step size is required.
                     */
                    dt_pkj = 0.5 * get_pkj_physiology_time_step.exec();
                    if (dt_myocardium - dt_pkj_sum < dt_pkj)
                        dt_pkj = dt_myocardium - dt_pkj_sum;

                    if (0 <= physical_time && physical_time <= 0.5)
                    {
                        apply_stimulus_pkj.exec(dt_pkj);
                    }
                    /**Strang splitting method. */
                    int ite_pkj_forward = 0;
                    while (ite_pkj_forward < reaction_step)
                    {
                        pkj_reaction_relaxation_forward.exec(0.5 * dt_pkj / Real(reaction_step));
                        ite_pkj_forward++;
                    }
                    /** 2nd Runge-Kutta scheme for diffusion. */
                    pkj_diffusion_relaxation.exec(dt_pkj);
                    // backward reaction
                    int ite_pkj_backward = 0;
                    while (ite_pkj_backward < reaction_step)
                    {
                        pkj_reaction_relaxation_backward.exec(0.5 * dt_pkj / Real(reaction_step));
                        ite_pkj_backward++;
                    }

                    dt_pkj_sum += dt_pkj;
                }

                /**Strang splitting method. */
                int ite_forward = 0;
                while (ite_forward < reaction_step)
                {
                    myocardium_reaction_relaxation_forward.exec(0.5 * dt_myocardium / Real(reaction_step));
                    ite_forward++;
                }
                /** 2nd Runge-Kutta scheme for diffusion. */
                myocardium_diffusion_relaxation.exec(dt_myocardium);

                // backward reaction
                int ite_backward = 0;
                while (ite_backward < reaction_step)
                {
                    myocardium_reaction_relaxation_backward.exec(0.5 * dt_myocardium / Real(reaction_step));
                    ite_backward++;
                }

                active_stress_interpolation.exec();
                Real dt_muscle_sum = 0.0;
                while (dt_muscle_sum < dt_myocardium)
                {
                    dt_muscle = get_mechanics_time_step.exec();
                    if (dt_myocardium - dt_muscle_sum < dt_muscle)
                        dt_muscle = dt_myocardium - dt_muscle_sum;
                    stress_relaxation_first_half.exec(dt_muscle);
                    constraint_holder.exec(dt_muscle);
                    stress_relaxation_second_half.exec(dt_muscle);
                    dt_muscle_sum += dt_muscle;
                }

                ite++;
                dt_myocardium = get_myocardium_physiology_time_step.exec();

                relaxation_time += dt_myocardium;
                integration_time += dt_myocardium;
                physical_time += dt_myocardium;
            }
            write_voltage.writeToFile(ite);
            write_displacement.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        interpolation_particle_position.exec();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
