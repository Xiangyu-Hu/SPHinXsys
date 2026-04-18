/**
 * @file 	muscle_solid_contact.cpp
 * @brief 	This is the test for muscle compression with our new contact model.
 * @details A soft body is in contact with rigid moving plate coupled with simbody.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.04;
Real PL = 0.1;
Real global_resolution = L / 12.0;
Real BW = global_resolution * 4;
Vecd halfsize_myocardium(0.5 * L, 0.5 * L, 0.5 * L);
Vecd translation_myocardium(0.5 * L, 0.0, 0.0);
Vecd halfsize_stationary_plate(0.5 * BW, 0.5 * L + BW, 0.5 * L + BW);
Vecd translation_stationary_plate(-0.5 * BW, 0.0, 0.0);
Vecd halfsize_moving_plate(0.5 * BW, 0.5 * PL, 0.5 * PL);
Vecd translation_moving_plate(L + BW, 0.0, 0.0);
BoundingBoxd system_domain_bounds(Vecd(-BW, -0.5 * PL, -0.5 * PL),
                                 Vecd(2.0 * L + BW, 0.5 * PL, 0.5 * PL));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_s = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
class Myocardium : public ComplexShape
{
  public:
    explicit Myocardium(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_myocardium), halfsize_myocardium);
        add<GeometricShapeBox>(Transform(translation_stationary_plate), halfsize_stationary_plate);
    }
};

class MovingPlate : public ComplexShape
{
  public:
    explicit MovingPlate(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_moving_plate), halfsize_moving_plate);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    SolidBody myocardium_body(sph_system, makeShared<Myocardium>("MyocardiumBody"));
    myocardium_body.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    myocardium_body.generateParticles<BaseParticles, Lattice>();

    SolidBody moving_plate(sph_system, makeShared<MovingPlate>("MovingPlate"));
    moving_plate.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    moving_plate.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation myocardium_body_inner(myocardium_body);
    SurfaceContactRelation myocardium_plate_contact(myocardium_body, {&moving_plate});
    SurfaceContactRelation plate_myocardium_contact(moving_plate, {&myocardium_body});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    /** Corrected configuration. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(myocardium_body_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(myocardium_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(myocardium_body_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactFactorSummation> myocardium_update_contact_density(myocardium_plate_contact);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> plate_update_contact_density(plate_myocardium_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> myocardium_compute_solid_contact_forces(myocardium_plate_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> plate_compute_solid_contact_forces(plate_myocardium_contact);
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(myocardium_body);
    /** Constrain the holder. */
    GeometricShapeBox holder_shape(Transform(translation_stationary_plate), halfsize_stationary_plate, "Holder");
    BodyRegionByParticle holder(myocardium_body, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        muscle_damping(0.1, myocardium_body_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_myocardium_body_kinetic_energy(myocardium_body);
    //----------------------------------------------------------------------
    //	The multi body system from simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SimTK::CableTrackerSubsystem cables(MBsystem);
    /** mass properties of the fixed spot. */
    GeometricShapeBox moving_plate_shape(Transform(translation_moving_plate), halfsize_moving_plate, "Plate");
    SimpleDynamics<NormalDirectionFromBodyShape> moving_plate_normal_direction(moving_plate);
    SolidBodyPartForSimbody plate_multibody(moving_plate, moving_plate_shape);
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
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    moving_plate_normal_direction.exec();
    corrected_configuration.exec();

    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 0.1;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
                write_myocardium_body_kinetic_energy.writeToFile(ite);
            }
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
            physical_time += dt;

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

    if (sph_system.GenerateRegressionData())
    {
        write_myocardium_body_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_myocardium_body_kinetic_energy.testResult();
    }
    return 0;
}
