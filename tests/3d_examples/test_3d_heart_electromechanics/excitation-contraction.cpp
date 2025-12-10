/**
 * @file 	excitation-contraction.cpp
 * @brief 	This is the case studying the electromechanics on a bi-ventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
/** Geometry parameter. */
/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/heart-new.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
/** Parameters and physical properties. */
Vec3d domain_lower_bound(-55.0 * length_scale, -75.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0; /**< Initial particle spacing. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 100 * stress_scale;
std::array<Real, 4> a0 = {Real(496.0) * stress_scale, Real(15196.0) * stress_scale, Real(3283.0) * stress_scale, Real(662.0) * stress_scale};
std::array<Real, 4> b0 = {Real(7.209), Real(20.417), Real(11.176), Real(9.466)};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 0.8;
Real bias_coeff = 0.0;
/** Electrophysiology parameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);

namespace SPH
{
/**
 * Define heart geometry
 */
class Heart : public ComplexShape
{
  public:
    explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
        add<TriangleMeshShapeSTL>(full_path_to_stl_file, translation, length_scale);
    }
};
using FiberDirectionDiffusionRelaxation =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;
/** Imposing diffusion boundary condition */
class DiffusionBCs : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit DiffusionBCs(BodyPartByParticle &body_part, const std::string &species_name)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name)) {};
    virtual ~DiffusionBCs() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);

        Real angle = face_norm.dot(center_norm);
        if (angle >= 0.0)
        {
            phi_[index_i] = 1.0;
        }
        else
        {
            if (pos_[index_i][1] < -getSPHAdaptation().ReferenceSpacing())
                phi_[index_i] = 0.0;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberAndSheetDirections : public LocalDynamics
{
  protected:
    LocallyOrthotropicMuscle &muscle_material_;
    Vecd *pos_;
    Real *phi_;
    Real beta_epi_, beta_endo_;
    Vecd center_line_vector_; // parallel to the ventricular centerline and pointing  apex-to-base

  public:
    explicit ComputeFiberAndSheetDirections(SPHBody &sph_body, const std::string &species_name)
        : LocalDynamics(sph_body),
          muscle_material_(DynamicCast<LocallyOrthotropicMuscle>(this, sph_body_->getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name))
    {
        center_line_vector_ = Vecd(0.0, 1.0, 0.0);
        beta_epi_ = -(70.0 / 180.0) * M_PI;
        beta_endo_ = (80.0 / 180.0) * M_PI;
    };
    virtual ~ComputeFiberAndSheetDirections() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        /**
         * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
         * 		Present  doi.org/10.1016/j.cma.2016.05.031
         */
        /** Probe the face norm from level set field. */
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);
        if (face_norm.dot(center_norm) <= 0.0)
        {
            face_norm = -face_norm;
        }
        /** Compute the centerline's projection on the plane orthogonal to face norm. */
        Vecd circumferential_direction = getCrossProduct(center_line_vector_, face_norm);
        Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
        /** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
        Real beta = (beta_epi_ - beta_endo_) * phi_[index_i] + beta_endo_;
        /** Compute the rotation matrix through Rodrigues rotation formulation. */
        Vecd f_0 = cos(beta) * cd_norm + sin(beta) * getCrossProduct(face_norm, cd_norm) +
                   face_norm.dot(cd_norm) * (1.0 - cos(beta)) * face_norm;

        if (pos_[index_i][1] < -getSPHAdaptation().ReferenceSpacing())
        {
            muscle_material_.local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
            muscle_material_.local_s0_[index_i] = face_norm;
        }
        else
        {
            muscle_material_.local_f0_[index_i] = Vecd::Zero();
            muscle_material_.local_s0_[index_i] = Vecd::Zero();
        }
    };
};
//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
  public:
    MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
    {
        Real l = domain_upper_bound[0] - domain_lower_bound[0];
        Real w = domain_upper_bound[2] - domain_lower_bound[2];
        halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
        resolution_ = 20;
        translation_ = Vec3d(-10.0 * length_scale, -1.0 * dp_0, 0.0);
    }
};
//	application dependent initial condition
class ApplyStimulusCurrentSI : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentSI(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (-30.0 * length_scale <= pos_[index_i][0] && pos_[index_i][0] <= -15.0 * length_scale)
        {
            if (-2.0 * length_scale <= pos_[index_i][1] && pos_[index_i][1] <= 0.0)
            {
                if (-3.0 * length_scale <= pos_[index_i][2] && pos_[index_i][2] <= 3.0 * length_scale)
                {
                    voltage_[index_i] = 0.92;
                }
            }
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentSII : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentSII(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (0.0 <= pos_[index_i][0] && pos_[index_i][0] <= 6.0 * length_scale)
        {
            if (-6.0 * length_scale <= pos_[index_i][1])
            {
                if (12.0 * length_scale <= pos_[index_i][2])
                {
                    voltage_[index_i] = 0.95;
                }
            }
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    observation_points.push_back(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale));
    observation_points.push_back(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale));
    observation_points.push_back(Vecd(0.0, -70.0 * length_scale, 0.0));
    return observation_points;
};
} // namespace SPH

/**
 * The main program.
 */
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	SPHSystem section
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);         // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    //----------------------------------------------------------------------
    //	SPH Particle relaxation for body fitted particle distribution.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        SolidBody herat_model(sph_system, makeShared<Heart>("HeartModel"));
        herat_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
        herat_model.defineClosure<LocallyOrthotropicMuscle, IsotropicDiffusion>(
            ConstructArgs(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0),
            ConstructArgs(diffusion_species_name, diffusion_coeff));
        herat_model.generateParticles<BaseParticles, Lattice>();
        /** topology */
        InnerRelation herat_model_inner(herat_model);
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_particles(herat_model);
        RelaxationStepInner relaxation_step_inner(herat_model_inner);
        //----------------------------------------------------------------------
        //	Relaxation output
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_herat_model_state_to_vtp({herat_model});
        ReloadParticleIO write_particle_reload_files(herat_model);
        write_particle_reload_files.addToReload<Vecd>(herat_model, "Fiber");
        write_particle_reload_files.addToReload<Vecd>(herat_model, "Sheet");
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
        //----------------------------------------------------------------------
        //	Diffusion process to initialize fiber direction
        //----------------------------------------------------------------------
        GetDiffusionTimeStepSize get_time_step_size(herat_model);
        FiberDirectionDiffusionRelaxation diffusion_relaxation(herat_model_inner);
        SimpleDynamics<ComputeFiberAndSheetDirections> compute_fiber_sheet(herat_model, diffusion_species_name);
        BodySurface surface_part(herat_model);
        SimpleDynamics<DiffusionBCs> impose_diffusion_bc(surface_part, diffusion_species_name);
        impose_diffusion_bc.exec();
        write_herat_model_state_to_vtp.addToWrite<Real>(herat_model, diffusion_species_name);
        write_herat_model_state_to_vtp.writeToFile(ite);

        int diffusion_step = 100;
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
    SolidBody mechanics_heart(sph_system, makeShared<Heart>("MechanicalHeart"));
    mechanics_heart.defineMaterial<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? mechanics_heart.generateParticles<BaseParticles, Reload>("HeartModel")
        : mechanics_heart.generateParticles<BaseParticles, Lattice>();

    SolidBody physiology_heart(sph_system, makeShared<Heart>("PhysiologyHeart"));
    AlievPanfilowModel aliev_panfilow_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
    physiology_heart.defineClosure<Solid, MonoFieldElectroPhysiology<LocalDirectionalDiffusion>>(
        Solid(), ConstructArgs(&aliev_panfilow_model, ConstructArgs(diffusion_coeff, bias_coeff, fiber_direction)));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? physiology_heart.generateParticles<BaseParticles, Reload>("HeartModel")
        : physiology_heart.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	SPH Observation section
    //----------------------------------------------------------------------
    ObserverBody voltage_observer(sph_system, "VoltageObserver");
    voltage_observer.generateParticles<ObserverParticles>(createObservationPoints());

    ObserverBody myocardium_observer(sph_system, "MyocardiumObserver");
    myocardium_observer.generateParticles<ObserverParticles>(createObservationPoints());
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
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_excitation(physiology_heart_inner);
    // Time step size calculation.
    GetDiffusionTimeStepSize get_physiology_time_step(physiology_heart);
    // Diffusion process for diffusion body.
    electro_physiology::ElectroPhysiologyDiffusionInnerRK2<LocalDirectionalDiffusion>
        diffusion_relaxation(physiology_heart_inner);
    // Solvers for ODE system.
    electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_heart, aliev_panfilow_model);
    electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_heart, aliev_panfilow_model);
    //	Apply the Iron stimulus.
    SimpleDynamics<ApplyStimulusCurrentSI> apply_stimulus_s1(physiology_heart);
    SimpleDynamics<ApplyStimulusCurrentSII> apply_stimulus_s2(physiology_heart);
    //  Active mechanics.
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration_contraction(mechanics_body_inner);
    InteractionDynamics<CorrectInterpolationKernelWeights> correct_kernel_weights_for_interpolation(mechanics_body_contact);
    /** Interpolate the active contract stress from electrophysiology body. */
    InteractionDynamics<InterpolatingAQuantity<Real>>
        active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
    /** Interpolate the particle position in physiology_heart  from mechanics_heart. */
    // TODO: this is a bug, we should interpolate displacement other than position.
    InteractionDynamics<InterpolatingAQuantity<Vecd>>
        interpolation_particle_position(physiology_heart_contact, "Position", "Position");
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(mechanics_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(mechanics_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> get_mechanics_time_step(mechanics_heart);
    /** Constrain region of the inserted body. */
    MuscleBaseShapeParameters muscle_base_parameters;
    BodyRegionByParticle muscle_base(mechanics_heart, makeShared<TriangleMeshShapeBrick>(muscle_base_parameters, "Holder"));
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(muscle_base);
    //----------------------------------------------------------------------
    //	SPH Output section
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Real>(physiology_heart, "Voltage");
    write_states.addToWrite<Real>(physiology_heart, "GateVariable");
    write_states.addToWrite<Real>(physiology_heart, "ActiveContractionStress");
    write_states.addToWrite<Real>(mechanics_heart, "ActiveContractionStress");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_voltage("Voltage", voltage_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", myocardium_observer_contact);
    //----------------------------------------------------------------------
    //	 Pre-simulation.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correct_configuration_excitation.exec();
    correct_configuration_contraction.exec();
    correct_kernel_weights_for_interpolation.exec();
    /** Output initial states and observations */
    write_states.writeToFile(0);
    write_voltage.writeToFile(0);
    write_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	 Physical parameters for main loop.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
                              << "	dt = " << dt
                              << "	dt_s = " << dt_s << "\n";
                }
                /** Apply stimulus excitation. */
                if (0 <= physical_time && physical_time <= 0.5)
                {
                    apply_stimulus_s1.exec(dt);
                }
                /** Single spiral wave. */
                // if( 60 <= physical_time
                // 	&&  physical_time <= 65)
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
                    dt_s_sum += dt_s;
                }

                ite++;
                dt = get_physiology_time_step.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
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

    if (sph_system.GenerateRegressionData())
    {
        std::cout << "Generate regression data ..." << std::endl;
        write_voltage.generateDataBase(1.0e-3);
        write_displacement.generateDataBase(1.0e-3);
    }
    else
    {
        std::cout << "Regression test ..." << std::endl;
        write_voltage.testResult();
        write_displacement.testResult();
    }

    return 0;
}
