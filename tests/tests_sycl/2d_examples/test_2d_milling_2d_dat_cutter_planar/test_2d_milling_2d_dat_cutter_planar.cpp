/**1
 * @file    milling_2d_dat_cutter_planar.cpp
 * @brief   2D milling demo with a .dat-profile cutter:
 *          - rigid cutter imported from a single 2D .dat contour
 *          - plate + cutter particle relaxation + reload workflow
 *          - deformable J2 plate clamped at both ends
 *          - contact via RepulsionFactor / RepulsionForceCK (UL style)
 *          - rigid cutter driven by Simbody Planar mobilizer (rotation + feed)
 *
 * Workflow:
 *   1) First run:
 *        ./case_name --relax=true
 *      This relaxes BOTH plate and cutter particles and writes reload files.
 *
 *   2) Second run:
 *        ./case_name --relax=false --reload=true
 *      This reloads the relaxed plate/cutter and runs the formal milling simulation.
 */

#include "sphinxsys.h"
#include "general_constraint_ck.h"

using namespace SPH;

//----------------------------------------------------------------------
// 1) Geometry / resolution (2D)
//----------------------------------------------------------------------

// Plate geometry
Real PL = 0.02;                  // Plate length scale in x
Real PH = 0.004;                 // Plate height scale in y
Real dp_0 = PL / 25.0;           // Particle spacing

Real BW = 10.0 * dp_0;           // Margin
Real DL = 10.0 * PL;
Real DH = 30.0 * PL;

// 2D computational domain bounds
Vec2d domain_lower_bound(-DL, -DH);
Vec2d domain_upper_bound(DL, DH);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

// Cutter imported from a single .dat file
std::string cutter_dat_file = "./input/2D.dat";

// Initial cutter placement in simulation space
Vec2d translation_holder_cutter(0.0, -1.4 * PL);

// Plate represented as a simple rectangular box
Vec2d halfsize_holder_plate(3.0 * PL, 1.0 * PL);
Vec2d translation_holder_plate(0.0, -0.051);

// Clamp length at both plate ends
Real clamp_len = 6.0 * dp_0;

//----------------------------------------------------------------------
// 2) Material
//----------------------------------------------------------------------

Real rho0_s = 2700.0;
Real poisson = 0.30;
Real Young_plate = 78.2e9;
Real yield_plate = 0.29e9;

// Numerical sound speed
Real c0 = sqrt(Young_plate / (3.0 * (1.0 - 2.0 * poisson) * rho0_s));

// Time step estimate
Real U_max = 373.0;

//----------------------------------------------------------------------
// 3) Cutter shape imported from .dat
//----------------------------------------------------------------------

class CutterImportModel : public MultiPolygonShape
{
  public:
    explicit CutterImportModel(const std::string &name)
        : MultiPolygonShape(name)
    {
        multi_polygon_.addAPolygonFromFile(cutter_dat_file, GeometricOps::add);
    }
};

//----------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------

int main(int ac, char *av[])
{
    SPHSystem system(system_domain_bounds, dp_0);

    // The two-stage workflow can be controlled by command-line options
    system.setRunParticleRelaxation(true);
    system.handleCommandlineOptions(ac, av);

    // ------------------------------------------------------------
    // Shapes
    // ------------------------------------------------------------
    auto &plate_shape =
        system.addShape<GeometricShapeBox>(
            Transform(translation_holder_plate), halfsize_holder_plate, "Plate");

    auto &cutter_shape =
        system.addShape<CutterImportModel>("Cutter");

    //==================================================================
    // Plate + cutter relaxation branch
    //==================================================================
    if (system.RunParticleRelaxation())
    {
        std::cout << "[Plate + Cutter Relaxation] start...\n";

        // Independent relaxation system
        RelaxationSystem relaxation_system(system_domain_bounds, dp_0);

        // Add both plate and cutter into the relaxation system
        auto &plate_relax  = relaxation_system.addBody<RealBody>(plate_shape);
        auto &cutter_relax = relaxation_system.addBody<SolidBody>(cutter_shape);

        // Level-set definitions
       LevelSetShape &plate_level_set =
          plate_relax.defineBodyLevelSetShape(par_ck)
            .correctLevelSetSign()
            .writeLevelSet();

       LevelSetShape &cutter_level_set =
          cutter_relax.defineBodyLevelSetShape(par_ck)
            .correctLevelSetSign()
            .writeLevelSet();

        // Generate lattice particles
        plate_relax.generateParticles<BaseParticles, Lattice>();
        cutter_relax.generateParticles<BaseParticles, Lattice>();

        // Near-surface body parts
        auto &near_plate_surface  = plate_relax.addBodyPart<NearShapeSurface>();
        auto &near_cutter_surface = cutter_relax.addBodyPart<NearShapeSurface>();

        // Inner relations
        auto &plate_inner_rel  = relaxation_system.addInnerRelation(plate_relax);
        auto &cutter_inner_rel = relaxation_system.addInnerRelation(cutter_relax);

        // Relaxation solver
        SPHSolver relax_solver(relaxation_system);
        auto &relax_main = relax_solver.addParticleMethodContainer(par_ck);
        auto &relax_host = relax_solver.addParticleMethodContainer(par_host);

        // Randomize initial particle positions
        auto &rand_plate =
            relax_host.addStateDynamics<RandomizeParticlePositionCK>(plate_relax);
        auto &rand_cutter =
            relax_host.addStateDynamics<RandomizeParticlePositionCK>(cutter_relax);

        // Configuration updates
        auto &update_plate_cell_linked_list =
            relax_main.addCellLinkedListDynamics(plate_relax);
        auto &update_plate_inner_relation =
            relax_main.addRelationDynamics(plate_inner_rel);

        auto &update_cutter_cell_linked_list =
            relax_main.addCellLinkedListDynamics(cutter_relax);
        auto &update_cutter_inner_relation =
            relax_main.addRelationDynamics(cutter_inner_rel);

        // Relaxation residuals
        auto &plate_relaxation_residual =
            relax_main.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(plate_inner_rel)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(plate_relax, plate_level_set);

        auto &cutter_relaxation_residual =
            relax_main.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(cutter_inner_rel)
                .addPostStateDynamics<LevelsetKernelGradientIntegral>(cutter_relax, cutter_level_set);

        // Particle position updates
        auto &update_plate_particle_position =
            relax_main.addStateDynamics<PositionRelaxationCK>(plate_relax);
        auto &update_cutter_particle_position =
            relax_main.addStateDynamics<PositionRelaxationCK>(cutter_relax);

        // Level-set bounding
        auto &plate_level_set_bounding =
            relax_main.addStateDynamics<LevelsetBounding>(near_plate_surface);
        auto &cutter_level_set_bounding =
            relax_main.addStateDynamics<LevelsetBounding>(near_cutter_surface);

        // Relaxation scaling
        auto &plate_relaxation_scaling =
            relax_main.addReduceDynamics<RelaxationScalingCK>(plate_relax);
        auto &cutter_relaxation_scaling =
            relax_main.addReduceDynamics<RelaxationScalingCK>(cutter_relax);

        // Output
        auto &body_state_recorder =
            relax_main.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(relaxation_system);

        auto &write_reload =
            relax_main.addIODynamics<ReloadParticleIOCK>(
                StdVec<SPHBody *>{&plate_relax, &cutter_relax});

        // Cutter normal direction for later wall-contact usage
        auto &cutter_normal =
            relax_host.addStateDynamics<NormalFromBodyShapeCK>(cutter_relax);
        write_reload.addToReload<Vecd>(cutter_relax, "NormalDirection");

        // Randomize before relaxation
        rand_plate.exec(0.25);
        rand_cutter.exec(0.25);

        body_state_recorder.writeToFile(0);

        int ite_p = 0;
        int ite_max = 2000;
        int output_interval = 100;

        while (ite_p < ite_max)
        {
            // Plate relaxation
            update_plate_cell_linked_list.exec();
            update_plate_inner_relation.exec();
            plate_relaxation_residual.exec();
            Real plate_relaxation_step = plate_relaxation_scaling.exec();
            update_plate_particle_position.exec(plate_relaxation_step);
            plate_level_set_bounding.exec();

            // Cutter relaxation
            update_cutter_cell_linked_list.exec();
            update_cutter_inner_relation.exec();
            cutter_relaxation_residual.exec();
            Real cutter_relaxation_step = cutter_relaxation_scaling.exec();
            update_cutter_particle_position.exec(cutter_relaxation_step);
            cutter_level_set_bounding.exec();

            ite_p += 1;

            if (ite_p % output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "Relaxation steps N = " << ite_p << "\n";
                body_state_recorder.writeToFile(ite_p);
            }
        }

        cutter_normal.exec();
        write_reload.writeToFile();

        std::cout << "[Plate + Cutter Relaxation] done. Reload files written.\n";
        return 0;
    }

    //==================================================================
    // Formal simulation branch
    //==================================================================

    // ------------------------------------------------------------
    // Bodies
    // ------------------------------------------------------------
    auto &plate  = system.addBody<RealBody>(plate_shape);
    auto &cutter = system.addBody<SolidBody>(cutter_shape);

    // ------------------------------------------------------------
    // Materials + particles
    // ------------------------------------------------------------

    // Cutter: reload relaxed particles + normals
    cutter.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Young_plate, poisson);
    cutter.generateParticles<BaseParticles, Reload>(cutter.getName())
     .reloadExtraVariable<Vecd>("NormalDirection");

    // Plate: reload relaxed particles as well
    plate.defineMaterial<J2Plasticity>(rho0_s, c0, Young_plate, poisson, yield_plate);
    plate.generateParticles<BaseParticles, Reload>(plate.getName());

    //------------------------------------------------------------------
    // Observer
    //------------------------------------------------------------------
    StdVec<Vecd> observation_location = {
        Vecd(translation_holder_cutter[0], translation_holder_cutter[1])};
    auto &observer = system.addBody<ObserverBody>("Observer");
    observer.generateParticles<ObserverParticles>(observation_location);

    //------------------------------------------------------------------
    // Relations
    //------------------------------------------------------------------
    auto &cutter_inner = system.addInnerRelation(cutter);
    auto &plate_inner  = system.addInnerRelation(plate);
    auto &plate_cutter_contact = system.addContactRelation(plate, cutter);
    auto &observer_plate_contact = system.addContactRelation(observer, plate);

    //------------------------------------------------------------------
    // Solver + method containers
    //------------------------------------------------------------------
    SPHSolver sph_solver(system);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);

    // Initialize required state variables
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(
        plate, "Velocity", Vec2d::Zero()).exec();
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(
        plate, "ForcePrior", Vec2d::Zero()).exec();

    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(
        cutter, "Velocity", Vec2d::Zero()).exec();
    host_methods.addStateDynamics<VariableAssignment, ConstantValue<Vecd>>(
        cutter, "ForcePrior", Vec2d::Zero()).exec();

    auto &cutter_normal_direction =
        host_methods.addStateDynamics<NormalFromBodyShapeCK>(cutter);

    //------------------------------------------------------------------
    // Fixed regions on both plate ends (2D AlignedBoxByParticle)
    //------------------------------------------------------------------
    Vec2d plate_half = halfsize_holder_plate;
    Vec2d plate_center = translation_holder_plate;

    // Plate outer bounds in world coordinates
    Real x_min = plate_center[0] - plate_half[0];
    Real x_max = plate_center[0] + plate_half[0];
    Real y_min = plate_center[1] - plate_half[1];
    Real y_max = plate_center[1] + plate_half[1];

    // Left clamp region
    Vec2d left_lower(x_min, y_min);
    Vec2d left_upper(x_min + clamp_len, y_max);
    Vec2d left_halfsize = 0.5 * (left_upper - left_lower);
    Vec2d left_center   = 0.5 * (left_upper + left_lower);
    AlignedBox left_box(xAxis, Transform(left_center), left_halfsize);
    AlignedBoxByParticle left_part(plate, left_box);

    // Right clamp region
    Vec2d right_lower(x_max - clamp_len, y_min);
    Vec2d right_upper(x_max, y_max);
    Vec2d right_halfsize = 0.5 * (right_upper - right_lower);
    Vec2d right_center   = 0.5 * (right_upper + right_lower);
    AlignedBox right_box(xAxis, Transform(right_center), right_halfsize);
    AlignedBoxByParticle right_part(plate, right_box);

    auto &fix_left_velocity =
        main_methods.addStateDynamics<ConstantConstraintCK<BodyPartByParticle, Vecd>>(
            left_part, "Velocity", Vec2d::Zero());

    auto &fix_right_velocity =
        main_methods.addStateDynamics<ConstantConstraintCK<BodyPartByParticle, Vecd>>(
            right_part, "Velocity", Vec2d::Zero());

    //------------------------------------------------------------------
    // Configuration update group
    //------------------------------------------------------------------
    ParticleDynamicsGroup update_configuration;
    update_configuration.add(&main_methods.addCellLinkedListDynamics(plate));
    update_configuration.add(&main_methods.addCellLinkedListDynamics(cutter));

    auto &update_plate_inner_and_contact =
        main_methods.addRelationDynamics(plate_inner, plate_cutter_contact);
    update_configuration.add(&update_plate_inner_and_contact);

    auto &update_observer_contact =
        main_methods.addRelationDynamics(observer_plate_contact);

    //------------------------------------------------------------------
    // Plate dynamics
    //------------------------------------------------------------------
    auto &plate_advection_setup =
        main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(plate);
    auto &plate_update_pos =
        main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(plate);

    auto &plate_correction =
        main_methods.addInteractionDynamicsWithUpdate<LinearCorrectionMatrix>(plate_inner);

    ParticleDynamicsGroup plate_shear_force;
    plate_shear_force.add(
        &main_methods.addInteractionDynamics<LinearGradient, Vecd>(plate_inner, "Velocity"));
    plate_shear_force.add(
        &main_methods.addInteractionDynamicsOneLevel<
            continuum_dynamics::ShearIntegration, J2Plasticity>(plate_inner));

    auto &plate_acoustic_1st =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep1stHalf,
            DissipativeRiemannSolverCK, NoKernelCorrectionCK>(plate_inner);

    auto &plate_acoustic_2nd =
        main_methods.addInteractionDynamicsOneLevel<
            fluid_dynamics::AcousticStep2ndHalf,
            DissipativeRiemannSolverCK, NoKernelCorrectionCK>(plate_inner);

    auto &contact_factor =
        main_methods.addInteractionDynamics<solid_dynamics::RepulsionFactor>(plate_cutter_contact);

    auto &contact_force =
        main_methods.addInteractionDynamicsWithUpdate<
            solid_dynamics::RepulsionForceCK, Wall>(plate_cutter_contact);

    //------------------------------------------------------------------
    // Simbody planar cutter driving
    //------------------------------------------------------------------
    Real omega_z = 2.0 * Pi * 200.0;   // Angular velocity [rad/s]
    Real feed_vx = 0.0;
    Real feed_vy = -0.030;             // Feed velocity in y

    Real end_time = 0.009;

    std::cout << "[Drive] omega_z=" << omega_z
              << " feed=(" << feed_vx << "," << feed_vy << ")"
              << " U_max~" << U_max << " end_time=" << end_time << "\n";

    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);

    SolidBodyPartForSimbody cutter_multibody(cutter, cutter_shape);
    SimTK::Body::Rigid cutter_info(*cutter_multibody.body_part_mass_properties_);

    SimTK::MobilizedBody::Planar cutter_mob(
        matter.Ground(),
        SimTK::Transform(SimTKVec3(0.0, 0.0, 0.0)),
        // translation_holder_cutter[0], translation_holder_cutter[1], 0.0)),
        cutter_info,
        SimTK::Transform(SimTKVec3(0.0, 0.0, 0.0)));

    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    auto &constraint_cutter =
        main_methods.addStateDynamics<solid_dynamics::ConstraintBodyPartBySimBodyCK>(
            cutter_multibody, MBsystem, cutter_mob, integ);

    //------------------------------------------------------------------
    // Output / observer
    //------------------------------------------------------------------
    auto &body_state_recorder =
        main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(system);
    body_state_recorder.addToWrite<Real>(plate, "Pressure");
    body_state_recorder.addToWrite<Real>(plate, "Density");

    auto &observer_position =
        main_methods.addObserveRecorder<Vecd>("Position", observer_plate_contact);

    //------------------------------------------------------------------
    // TimeStepper
    //------------------------------------------------------------------
    Real total_physical_time = end_time;
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(total_physical_time);

    auto &plate_advection_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(plate, U_max, 0.2);
    auto &plate_acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(plate, 0.4);

    auto &advection_step =
        time_stepper.addTriggerByInterval(plate_advection_time_step.exec());

    size_t advection_steps = system.RestartStep() + 1;

    int screening_interval = 100;
    int observation_interval = screening_interval * 2;

    auto &state_recording =
        time_stepper.addTriggerByInterval(total_physical_time / 300.0);

    //------------------------------------------------------------------
    // Preparation before the main loop
    //------------------------------------------------------------------
    constraint_cutter.exec();
    cutter_normal_direction.exec();

    update_configuration.exec();
    plate_advection_setup.exec();
    contact_factor.exec();
    plate_correction.exec();

    fix_left_velocity.exec();
    fix_right_velocity.exec();

    body_state_recorder.writeToFile();
    update_observer_contact.exec();
    observer_position.writeToFile(advection_steps);

    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;

    TickCount t0 = TickCount::now();

    //------------------------------------------------------------------
    // Main time-stepping loop
    //------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        TickCount time_instance = TickCount::now();

        Real dt = time_stepper.incrementPhysicalTime(plate_acoustic_time_step);
        Real t_target = time_stepper.getPhysicalTime();

        // (A) Drive the cutter in-plane
        SimTK::State &s_adv = integ.updAdvancedState();
        SimTK::Vec3 u_cmd(omega_z, feed_vx, feed_vy);
        cutter_mob.setU(s_adv, u_cmd);
        MBsystem.realize(s_adv, SimTK::Stage::Velocity);

        if (t_target > integ.getState().getTime())
        {
            integ.stepTo(t_target);
        }

        // (B) Write back cutter pose and update normals
        constraint_cutter.exec();
        cutter_normal_direction.exec();

        // (C) Plate SPH dynamics + contact
        plate_shear_force.exec(dt);
        contact_force.exec();
        plate_acoustic_1st.exec(dt);
        plate_acoustic_2nd.exec(dt);

        fix_left_velocity.exec();
        fix_right_velocity.exec();

        interval_acoustic_step += TickCount::now() - time_instance;

        if (advection_step(plate_advection_time_step))
        {
            advection_steps++;

            plate_update_pos.exec();

            fix_left_velocity.exec();
            fix_right_velocity.exec();

            time_instance = TickCount::now();

            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << advection_steps
                          << "  Time=" << time_stepper.getPhysicalTime()
                          << "  advection_dt=" << advection_step.getInterval()
                          << "  acoustic_dt=" << time_stepper.getGlobalTimeStepSize() << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                update_observer_contact.exec();
                observer_position.writeToFile(advection_steps);
            }

            if (state_recording())
            {
                body_state_recorder.writeToFile();
            }

            interval_output += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            update_configuration.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            plate_advection_setup.exec();
            contact_factor.exec();
            plate_correction.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }

    //------------------------------------------------------------------
    // Summary
    //------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds.\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_advection_step = " << interval_advection_step.seconds() << "\n"
              << "interval_acoustic_step = " << interval_acoustic_step.seconds() << "\n"
              << "interval_updating_configuration = " << interval_updating_configuration.seconds() << "\n";

    if (system.GenerateRegressionData())
    {
        // Add generateDataBase here if needed
    }
    else
    {
        std::cout << "[Regression] baseline missing likely. Skip testResult in dev.\n";
    }

    return 0;
}
