/**
 * @file 	heart_volume_change.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D, including the volume change of the ventricles.
 * @author 	John Benjamin, Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */
#include "excitation_contraction.h"
using namespace SPH; // Namespace cite here.

/**
 * The main program.
 */
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	SPHSystem section
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, dp_0);
    GlobalStaticVariables::physical_time_ = 0.0;
    Real mechanical_time_ = 0.0;
    system.setRunParticleRelaxation(true); // Tag for run particle relaxation for body-fitted distribution
    system.setReloadParticles(true);       // Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    /** check whether run particle relaxation for body fitted particle distribution. */
    if (system.RunParticleRelaxation())
    {
        SolidBody herat_model(system, makeShared<Heart>("HeartModel"));
        herat_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment);
        herat_model.defineParticlesAndMaterial<FiberDirectionDiffusionParticles, FiberDirectionDiffusion>();
        herat_model.generateParticles<ParticleGeneratorLattice>();
        /** topology */
        InnerRelation herat_model_inner(herat_model);
        /** Random reset the relax solid particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_particles(herat_model);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(herat_model_inner);
        /** Time step for diffusion. */
        GetDiffusionTimeStepSize<FiberDirectionDiffusionParticles> get_time_step_size(herat_model);
        /** Diffusion process for diffusion body. */
        DiffusionRelaxation diffusion_relaxation(herat_model_inner);
        /** Compute the fiber and sheet after diffusion. */
        SimpleDynamics<ComputeFiberAndSheetDirections> compute_fiber_sheet(herat_model);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_herat_model_state_to_vtp(io_environment, {herat_model});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, herat_model);
        //----------------------------------------------------------------------
        //	Physics relaxation starts here.
        //----------------------------------------------------------------------
        random_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_herat_model_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        // From here the time stepping begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        int diffusion_step = 100;
        while (ite < relax_step)
        {
            relaxation_step_inner.exec();
            ite++;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_herat_model_state_to_vtp.writeToFile(ite);
            }
        }

        BodySurface surface_part(herat_model);
        /** constraint boundary condition for diffusion. */
        SimpleDynamics<DiffusionBCs> impose_diffusion_bc(surface_part, "Phi");
        impose_diffusion_bc.exec();

        write_herat_model_state_to_vtp.writeToFile(ite);

        Real dt = get_time_step_size.exec();
        while (ite <= diffusion_step + relax_step)
        {
            diffusion_relaxation.exec(dt);
            impose_diffusion_bc.exec();
            if (ite % 10 == 0)
            {
                std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
                write_herat_model_state_to_vtp.writeToFile(ite);
            }
            ite++;
        }
        compute_fiber_sheet.exec();
        ite++;
        write_herat_model_state_to_vtp.writeToFile(ite);
        compute_fiber_sheet.exec();
        write_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	SPH simulation section
    //----------------------------------------------------------------------
    /** create a SPH body, material and particles */
    SolidBody physiology_heart(system, makeShared<Heart>("PhysiologyHeart"));
    SharedPtr<AlievPanfilowModel> muscle_reaction_model_ptr = makeShared<AlievPanfilowModel>(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
    physiology_heart.defineParticlesAndMaterial<
        ElectroPhysiologyParticles, MonoFieldElectroPhysiology>(
        muscle_reaction_model_ptr, TypeIdentity<LocalDirectionalDiffusion>(), diffusion_coff, bias_coff, fiber_direction);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? physiology_heart.generateParticles<ParticleGeneratorReload>(io_environment, "HeartModel")
        : physiology_heart.generateParticles<ParticleGeneratorLattice>();

    /** create a SPH body, material and particles */
    SolidBody mechanics_heart(system, makeShared<Heart>("MechanicalHeart"));
    mechanics_heart.defineParticlesAndMaterial<
        ElasticSolidParticles, ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? mechanics_heart.generateParticles<ParticleGeneratorReload>(io_environment, "HeartModel")
        : mechanics_heart.generateParticles<ParticleGeneratorLattice>();
    auto myocardium_particles = dynamic_cast<ElasticSolidParticles *>(&mechanics_heart.getBaseParticles());

    //----------------------------------------------------------------------
    //	SPH Observation section
    //----------------------------------------------------------------------
    ObserverBody voltage_observer(system, "VoltageObserver");
    voltage_observer.generateParticles<HeartObserverParticleGenerator>();
    ObserverBody myocardium_observer(system, "MyocardiumObserver");
    myocardium_observer.generateParticles<HeartObserverParticleGenerator>();
    //----------------------------------------------------------------------
    //	SPHBody relation (topology) section
    //----------------------------------------------------------------------
    InnerRelation physiology_heart_inner(physiology_heart);
    InnerRelation mechanics_body_inner(mechanics_heart);
    ContactRelation physiology_heart_contact(physiology_heart, {&mechanics_heart});
    ContactRelation mechanics_body_contact(mechanics_heart, {&physiology_heart});
    ContactRelation voltage_observer_contact(voltage_observer, {&physiology_heart});
    ContactRelation myocardium_observer_contact(myocardium_observer, {&mechanics_heart});
    //----------------------------------------------------------------------
    //	SPH Method section
    //----------------------------------------------------------------------
    // Corrected configuration.
    InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration_excitation(physiology_heart_inner);
    // Time step size calculation.
    electro_physiology::GetElectroPhysiologyTimeStepSize get_physiology_time_step(physiology_heart);
    // Diffusion process for diffusion body.
    electro_physiology::ElectroPhysiologyDiffusionInnerRK2 diffusion_relaxation(physiology_heart_inner);
    // Solvers for ODE system.
    electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_heart);
    electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_heart);
    //	Apply the Iron stimulus.
    SimpleDynamics<ApplyStimulusCurrentSI> apply_stimulus_s1(physiology_heart);
    SimpleDynamics<ApplyStimulusCurrentSII> apply_stimulus_s2(physiology_heart);
    // Active mechanics.
    InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration_contraction(mechanics_body_inner);
    InteractionDynamics<CorrectInterpolationKernelWeights> correct_kernel_weights_for_interpolation(mechanics_body_contact);
    /** Interpolate the active contract stress from electrophysiology body. */
    InteractionDynamics<InterpolatingAQuantity<Real>>
        active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
    /** Interpolate the particle position in physiology_heart  from mechanics_heart. */
    // TODO: this is a bug, we should interpolate displacement other than position.
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_particle_position(physiology_heart_contact, "Position", "Position");
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> get_mechanics_time_step(mechanics_heart);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(mechanics_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(mechanics_body_inner);
    // initialize and update of normal direction
    SimpleDynamics<NormalDirectionFromBodyShape>(mechanics_heart).exec();
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> body_update_normal(mechanics_heart);
    /** Constrain region of the inserted body. */
    MuscleBaseShapeParameters muscle_base_parameters;
    BodyRegionByParticle muscle_base(mechanics_heart, makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(muscle_base);
    //----------------------------------------------------------------------
    //	SPH Output section
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_voltage("Voltage", io_environment, voltage_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", io_environment, myocardium_observer_contact);
    //----------------------------------------------------------------------
    //	 Pre-simulation.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    correct_configuration_excitation.exec();
    correct_configuration_contraction.exec();
    correct_kernel_weights_for_interpolation.exec();

    //----------------------------------------------------------------------
    //	 Surfaces and operations - must be after system initialized
    //----------------------------------------------------------------------
    MeshData myo_mesh;
    myo_mesh.load(full_path_to_myocardium, length_scale);
    myo_mesh.translate(translation);
    myo_mesh.initialize();
    MeshData lv_mesh;
    lv_mesh.load(full_path_to_lv, length_scale);
    lv_mesh.translate(translation);
    lv_mesh.initialize();
    MeshData rv_mesh;
    rv_mesh.load(full_path_to_rv, length_scale);
    rv_mesh.translate(translation);
    rv_mesh.initialize();
    MyocardiumSurfaces myo_srf(mechanics_heart, *myocardium_particles);
    Real smoothing_length = dp_0 * 1.15;
    myo_srf.init_surfaces(myo_mesh, lv_mesh, rv_mesh, 1.0, 1.0, smoothing_length);
    myo_srf.write_all_surfaces_as_obj("output/");
    SurfaceOperationsVentricle surface_ops_LV(*myocardium_particles, myo_srf.get_lv_ids(), mechanics_body_inner);
    SurfaceOperationsVentricle surface_ops_RV(*myocardium_particles, myo_srf.get_rv_ids(), mechanics_body_inner);
    // recording
    StdVec<Real> simulation_time,
        flow_rate_LV, delta_volume_LV,
        flow_rate_RV, delta_volume_RV;
    Real delta_V_total_LV = 0;
    Real delta_V_total_RV = 0;
    /** Output initial states and observations */
    write_states.writeToFile(0);
    write_voltage.writeToFile(0);
    write_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	 Physical parameters for main loop.
    //----------------------------------------------------------------------
    int screen_output_interval = 10;
    int ite = 0;
    int reaction_step = 2;
    Real end_time = 100;
    Real Ouput_T = end_time / 200.0;
    Real Observer_time = 0.01 * Ouput_T;
    Real dt = 0.0;   /**< Default acoustic time step sizes for physiology. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for mechanics. */
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    std::cout << "Main Loop Starts Here : "
              << "\n";
    /** Main loop starts here. */
    while (GlobalStaticVariables::physical_time_ < end_time)
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
                              << GlobalStaticVariables::physical_time_
                              << "	dt = " << dt
                              << "	dt_s = " << dt_s << "\n";
                }
                /** Apply stimulus excitation. */
                if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= 0.5)
                {
                    apply_stimulus_s1.exec(dt);
                }
                /** Single spiral wave. */
                // if( 60 <= GlobalStaticVariables::physical_time_
                // 	&&  GlobalStaticVariables::physical_time_ <= 65)
                // {
                // 	apply_stimulus_s2.exec(dt);
                // }
                /**Strong splitting method. */
                // forward reaction
                int ite_forward = 0;
                while (ite_forward < reaction_step)
                {
                    reaction_relaxation_forward.exec(0.5 * dt / Real(reaction_step));
                    ite_forward++;
                }
                /** 2nd Runge-Kutta scheme for diffusion. */
                diffusion_relaxation.exec(dt);

                // backward reaction
                int ite_backward = 0;
                while (ite_backward < reaction_step)
                {
                    reaction_relaxation_backward.exec(0.5 * dt / Real(reaction_step));
                    ite_backward++;
                }

                active_stress_interpolation.exec();

                Real dt_s_sum = 0.0;
                while (dt_s_sum < dt)
                {
                    dt_s = get_mechanics_time_step.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    stress_relaxation_first_half.exec(dt_s);
                    constraint_holder.exec(dt_s);
                    stress_relaxation_second_half.exec(dt_s);

                    // update normal and surface area
                    body_update_normal.exec();

                    { // updates for flow rate calculations for Windkessel
                        // surface areas
                        surface_ops_LV.update_srf_area();
                        surface_ops_RV.update_srf_area();
                        // update volume and flow rates
                        surface_ops_LV.update_flow_rate(dt_s);
                        surface_ops_LV.update_flow_acc(dt_s);
                        surface_ops_RV.update_flow_rate(dt_s);
                        surface_ops_RV.update_flow_acc(dt_s);
                    }
                    dt_s_sum += dt_s;
                    mechanical_time_ += dt_s;
                    { // Recording - can be plotted via writing to file
                        // time
                        simulation_time.push_back(mechanical_time_);
                        // LV
                        flow_rate_LV.push_back(surface_ops_LV.get_Q_current() * 1e-3);
                        delta_V_total_LV += surface_ops_LV.get_delta_V();
                        delta_volume_LV.push_back(delta_V_total_LV * 1e-3);
                        // RV
                        flow_rate_RV.push_back(surface_ops_RV.get_Q_current() * 1e-3);
                        delta_V_total_RV += surface_ops_RV.get_delta_V();
                        delta_volume_RV.push_back(delta_V_total_RV * 1e-3);
                    }
                }

                ite++;
                dt = get_physiology_time_step.exec();

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
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
    write_csv_files("LV_results.csv", "Time", "flow_rate_LV", "delta_volume_LV", simulation_time, flow_rate_LV, delta_volume_LV);
    write_csv_files("RV_results.csv", "Time", "flow_rate_RV", "delta_volume_RV", simulation_time, flow_rate_RV, delta_volume_RV);

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
