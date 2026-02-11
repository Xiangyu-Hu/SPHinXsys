/**
 * @file 	muscle_soft_body_contact.cpp
 * @brief 	This is the test for muscle compression with our new contact model.
 * Different particle resolutions are used for the two soft bodies that are in contact.
 * @author 	Chi Zhang, Bence Rochlitz and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;
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
    moving_plate.defineAdaptationRatios(1.15, 1.5);
    moving_plate.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    moving_plate.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation myocardium_body_inner(myocardium_body);
    InnerRelation moving_plate_inner(moving_plate);
    SurfaceContactRelation myocardium_plate_contact(myocardium_body, {&moving_plate});
    SurfaceContactRelation plate_myocardium_contact(moving_plate, {&myocardium_body});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(-100.0, 0.0, 0.0));
    SimpleDynamics<GravityForce<Gravity>> plate_initialize_constant_gravity(moving_plate, gravity);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(myocardium_body_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_2(moving_plate_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half(myocardium_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(myocardium_body_inner);
    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half_2(moving_plate_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half_2(moving_plate_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactFactorSummation> myocardium_update_contact_density(myocardium_plate_contact);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> plate_update_contact_density(plate_myocardium_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> myocardium_compute_solid_contact_forces(myocardium_plate_contact);
    InteractionWithUpdate<solid_dynamics::ContactForce> plate_compute_solid_contact_forces(plate_myocardium_contact);
    /** Constrain the holder. */
    GeometricShapeBox holder_shape(Transform(translation_stationary_plate), halfsize_stationary_plate, "Holder");
    BodyRegionByParticle holder(myocardium_body, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    /** Add spring constraint on the plate. */
    SimpleDynamics<solid_dynamics::SpringDamperConstraintParticleWise> spring_constraint(moving_plate, Vecd(0.2, 0, 0), 0.01);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> muscle_damping(0.2, myocardium_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>> plate_damping(0.2, moving_plate_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_plate_kinetic_energy(moving_plate);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    corrected_configuration_2.exec();
    plate_initialize_constant_gravity.exec();
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
            if (ite % 50 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
                write_plate_kinetic_energy.writeToFile(ite);
            }
            spring_constraint.exec();
            /** Contact model for myocardium. */
            myocardium_update_contact_density.exec();
            myocardium_compute_solid_contact_forces.exec();
            /** Contact model for plate. */
            plate_update_contact_density.exec();
            plate_compute_solid_contact_forces.exec();

            /** Stress relaxation and damping. */
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            muscle_damping.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            stress_relaxation_first_half_2.exec(dt);
            plate_damping.exec(dt);
            stress_relaxation_second_half_2.exec(dt);

            ite++;
            dt = sph_system.getSmallestTimeStepAmongSolidBodies();
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
        write_plate_kinetic_energy.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_plate_kinetic_energy.testResult();
    }

    return 0;
}
