/**
 * @file passive_cantilever.cpp
 * @brief This is the first example of cantilever
 * @author Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 6.0;
Real PH = 1.0;
Real PW = 1.0;
Real SL = 0.5;
Real global_resolution = PH / 12.0; /**< Initial particle spacing. */
Real BW = global_resolution * 4;    /**< Boundary width. */
Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_cantilever(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vecd(-SL - BW, -BW, -BW),
                                 Vecd(PL + BW, PH + BW, PH + BW));
// Observer location
StdVec<Vecd> observation_location = {Vecd(PL, PH, PW)};
/** For material properties of the solid. */
Real rho0_s = 1100.0;
Real poisson = 0.45;
Real Youngs_modulus = 1.7e7;
Real a = Youngs_modulus / (2.0 * (1.0 + poisson));
Real a_f = 0.0 * a;
std::array<Real, 4> a0 = {a, a_f, 0.0, 0.0};
std::array<Real, 4> b0 = {1.0, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real bulk_modulus = Youngs_modulus / 3.0 / (1.0 - 2.0 * poisson);

/** Define the cantilever shape. */
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
 * application dependent initial condition
 */
class CantileverInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit CantileverInitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][0] > 0.0)
        {
            vel_[index_i][1] = 5.0 * sqrt(3.0);
            vel_[index_i][2] = 5.0;
        }
    };
};
/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    /** create a Cantilever body, corresponding material, particles and reaction model. */
    SolidBody cantilever_body(sph_system, makeShared<Cantilever>("CantileverBody"));
    cantilever_body.defineMaterial<Muscle>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
    cantilever_body.generateParticles<BaseParticles, Lattice>();
    /** Define Observer. */
    ObserverBody cantilever_observer(sph_system, "CantileverObserver");
    cantilever_observer.generateParticles<ObserverParticles>(observation_location);

    /** topology */
    InnerRelation cantilever_body_inner(cantilever_body);
    ContactRelation cantilever_observer_contact(cantilever_observer, {&cantilever_body});

    /**
     * This section define all numerical methods will be used in this case.
     */
    /** Corrected configuration. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(cantilever_body_inner);
    /** active and passive stress relaxation. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(cantilever_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(cantilever_body_inner);
    SimpleDynamics<CantileverInitialCondition> initialization(cantilever_body);
    /** Time step size calculation. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(cantilever_body);
    /** Constrain the holder. */
    GeometricShapeBox holder_shape(Transform(translation_holder), halfsize_holder, "Holder");
    BodyRegionByParticle holder(cantilever_body, holder_shape);
    SimpleDynamics<FixBodyPartConstraint> constraint_holder(holder);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", cantilever_observer_contact);
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    /** apply initial condition */
    initialization.exec();
    corrected_configuration.exec();
    write_states.writeToFile(0);
    write_displacement.writeToFile(0);
    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 3.0;
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
            stress_relaxation_first_half.exec(dt);
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

    write_displacement.testResult();

    return 0;
}
