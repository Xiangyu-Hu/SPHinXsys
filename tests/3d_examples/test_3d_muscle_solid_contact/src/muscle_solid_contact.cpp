/**
 * @file 	muscle_solid_contact.cpp
 * @brief 	This is the test for muscle compression with our new contact model.
 * @details A soft body is in contact with rigid moving plate coupled with simbody.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real L = 0.04;
Real PL = 0.1;
Real resolution_ref = L / 12.0;
Real BW = resolution_ref * 4;
Vecd halfsize_myocardium(0.5 * L, 0.5 * L, 0.5 * L);
Vecd translation_myocardium(0.5 * L, 0.0, 0.0);
Vecd halfsize_stationary_plate(0.5 * BW, 0.5 * L + BW, 0.5 * L + BW);
Vecd translation_stationary_plate(-0.5 * BW, 0.0, 0.0);
Vecd halfsize_moving_plate(0.5 * BW, 0.5 * PL, 0.5 * PL);
Vecd translation_moving_plate(L + BW, 0.0, 0.0);

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -0.5 * PL, -0.5 * PL),
                                 Vecd(2.0 * L + BW, 0.5 * PL, 0.5 * PL));

/** For material properties of the solid. */
Real rho0_s = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;

/** Define the myocardium body shape. */
class Myocardium : public ComplexShape
{
  public:
    explicit Myocardium(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(translation_myocardium), halfsize_myocardium);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_stationary_plate), halfsize_stationary_plate);
    }
};
/**
 * @brief define the moving plate shape
 */
class MovingPlate : public ComplexShape
{
  public:
    explicit MovingPlate(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(translation_moving_plate), halfsize_moving_plate);
    }
};
/**
 *  The main program
 */
int main()
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, resolution_ref);
    /** Creat a Myocardium body, corresponding material, particles and reaction model. */
    SolidBody myocardium_body(system, makeShared<Myocardium>("MyocardiumBody"));
    myocardium_body.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    myocardium_body.generateParticles<ParticleGeneratorLattice>();
    /** Plate. */
    SolidBody moving_plate(system, makeShared<MovingPlate>("MovingPlate"));
    moving_plate.defineParticlesAndMaterial<SolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    moving_plate.generateParticles<ParticleGeneratorLattice>();
    /** topology */
    InnerRelation myocardium_body_inner(myocardium_body);
    SurfaceContactRelation myocardium_plate_contact(myocardium_body, {&moving_plate});
    SurfaceContactRelation plate_myocardium_contact(moving_plate, {&myocardium_body});
    /**
     * This section define all numerical methods will be used in this case.
     */
    /** initialize a time step */
    SimpleDynamics<TimeStepInitialization> myocardium_initialize_time_step(myocardium_body);
    SimpleDynamics<TimeStepInitialization> moving_plate_initialize_time_step(moving_plate);
    /** Corrected configuration. */
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(myocardium_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(myocardium_body);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(myocardium_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(myocardium_body_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> myocardium_update_contact_density(myocardium_plate_contact);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> plate_update_contact_density(plate_myocardium_contact);
    InteractionDynamics<solid_dynamics::ContactForce> myocardium_compute_solid_contact_forces(myocardium_plate_contact);
    InteractionDynamics<solid_dynamics::ContactForce> plate_compute_solid_contact_forces(plate_myocardium_contact);
    /** Constrain the holder. */
    BodyRegionByParticle holder(myocardium_body,
                                makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_stationary_plate), halfsize_stationary_plate, "Holder"));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        muscle_damping(0.1, myocardium_body_inner, "Velocity", physical_viscosity);
    /** Output */
    IOEnvironment io_environment(system);
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    /** Simbody interface. */
    /**
     * The multi body system from simbody.
     */
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SimTK::CableTrackerSubsystem cables(MBsystem);
    /** mass properties of the fixed spot. */
    SolidBodyPartForSimbody plate_multibody(moving_plate,
                                            makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_moving_plate), halfsize_moving_plate, "Plate"));
    /** Mass properties of the constrained spot.
     * SimTK::MassProperties(mass, center of mass, inertia)
     */
    SimTK::Body::Rigid rigid_info(*plate_multibody.body_part_mass_properties_);
    SimTK::MobilizedBody::Slider
        plateMBody(matter.Ground(), SimTK::Transform(SimTKVec3(0)), rigid_info, SimTK::Transform(SimTKVec3(0)));
    /** Gravity. */
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(Real(-100.0), 0.0, 0.0));
    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    /** Damper. */
    SimTK::Force::MobilityLinearDamper linear_damper(forces, plateMBody, SimTK::MobilizerUIndex(0), 20.0);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    /** Coupling between SimBody and SPH.*/
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_plate(plate_multibody, MBsystem, plateMBody, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_plate(plate_multibody, MBsystem, plateMBody, integ);
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    /** apply initial condition */
    corrected_configuration.exec();
    write_states.writeToFile(0);
    /** Setup physical parameters. */
    int ite = 0;
    Real end_time = 0.1;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            /** Gravity. */
            myocardium_initialize_time_step.exec();
            moving_plate_initialize_time_step.exec();
            /** Contact model for myocardium. */
            myocardium_update_contact_density.exec();
            myocardium_compute_solid_contact_forces.exec();
            /** Contact model for plate. */
            plate_update_contact_density.exec();
            plate_compute_solid_contact_forces.exec();
            {
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, plateMBody, force_on_plate.exec());
                integ.stepBy(dt);
                constraint_plate.exec();
            }
            /** Stress relaxation and damping. */
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            muscle_damping.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;

            myocardium_body.updateCellLinkedList();
            moving_plate.updateCellLinkedList();

            myocardium_plate_contact.updateConfiguration();
            plate_myocardium_contact.updateConfiguration();
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

    return 0;
}
