/**
 * @file passive_cantilever_neohookean.cpp
 * @brief This is the example of cantilever with simple neohookean tissue model
 * @author Bence Rochlitz, Chi Zhang  and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 0.1;
Real PH = 0.04;
Real PW = 0.04;
Real SL = 0.02;
Real resolution_ref = PH / 6.0; /**< Initial particle spacing. */
Real BW = resolution_ref * 4;   /**< Boundary width. */
Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_cantilever(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vecd(-SL, 0, 0), Vecd(PL, PH, PH));
// Observer location
StdVec<Vecd> observation_location = {Vecd(PL, PH, PW)};
/** For material properties of the solid. */
Real rho0_s = 1265.0;           // Gheorghe 2019
Real poisson = 0.45;            // nearly incompressible
Real Youngs_modulus = 5e4;      // Sommer 2015
Real physical_viscosity = 50.0; // physical damping, here we choose the same value as numerical viscosity
Real gravity_g = 9.8;           /**< Value of gravity. */
Real time_to_full_gravity = 0.0;

/** Define the cantilever body. */
class Cantilever : public ComplexShape
{
  public:
    explicit Cantilever(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(translation_cantilever), halfsize_cantilever);
        add<GeometricShapeBox>(Transform(translation_holder), halfsize_holder);
    }
};

/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif /** output environment. */

    /** create a Cantilever body, corresponding material, particles and reaction model. */
    SolidBody cantilever_body(sph_system, makeShared<Cantilever>("CantileverBody"));
    cantilever_body.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    cantilever_body.generateParticles<BaseParticles, Lattice>();
    /** Define Observer. */
    ObserverBody cantilever_observer(sph_system, "CantileverObserver");
    cantilever_observer.generateParticles<ObserverParticles>(observation_location);

    /** topology */
    InnerRelation cantilever_body_inner(cantilever_body);
    ContactRelation cantilever_observer_contact(cantilever_observer, {&cantilever_body});

    //-------- common particle dynamics ----------------------------------------
    IncreaseToFullGravity gravity(Vec3d(0.0, -gravity_g, 0.0), time_to_full_gravity);
    SimpleDynamics<GravityForce<IncreaseToFullGravity>> apply_time_dependent_gravity(cantilever_body, gravity);

    /**
     * This section define all numerical methods will be used in this case.
     */
    /** Corrected configuration. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(cantilever_body_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(cantilever_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(cantilever_body_inner);
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(cantilever_body);
    /** Constrain the holder. */
    GeometricShapeBox holder_shape(Transform(translation_holder), halfsize_holder, "Holder");
    BodyRegionByParticle holder(cantilever_body, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    DampingWithRandomChoice<InteractionSplit<DampingProjectionInner<Vec3d, FixedDampingRate>>>
        muscle_damping(0.1, cantilever_body_inner, "Velocity", physical_viscosity);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", cantilever_observer_contact);
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    write_states.writeToFile(0);
    write_displacement.writeToFile(0);
    /** Setup physical parameters. */
    int ite = 0;
    Real end_time = 1.0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
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
            }

            apply_time_dependent_gravity.exec();
            stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            muscle_damping.exec(dt);
            constraint_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;
        }
        write_displacement.writeToFile(ite);
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
        write_displacement.generateDataBase(1.0e-2);
    }
    else
    {
        write_displacement.testResult();
    }

    return 0;
}
