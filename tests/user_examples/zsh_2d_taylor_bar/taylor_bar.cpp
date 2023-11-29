/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscillation Ring example-update Lagrange           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for understanding SPH method for    *
* solid simulation based on update Lagrange method                            *
* In this case, the constraint of the Ring is implemented with                *
* internal constrained subregion.                                             *
* ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "all_continuum.h"
 
using namespace SPH;
//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.0032; /**< X-direction domain. */
Real PW = 0.0324; /**< Z-direction domain. */
Real particle_spacing_ref = PL / 10.0;
Real SL = particle_spacing_ref * 4.0; /**< Length of the holder is one layer particle. */
Real inner_circle_radius = PL;

Vecd domain_lower_bound(-4.0 * PL, -SL);
Vecd domain_upper_bound(4.0 * PL, 2.0 * PW);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);

//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; /**< Reference density. */
Real poisson = 0.35;  /**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;

Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));

class Wall : public ComplexShape
{
public:
    explicit Wall(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_holder(3.0 * PL, 0.5 * SL);
        Vecd translation_holder(0.0, -0.5 * SL);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_holder), halfsize_holder);
    }
};

/** Define the body. */
class Column : public ComplexShape
{
public:
    explicit Column(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column(0.0, 0.6 * PW);
        Vecd halfsize_column(0.003, 0.5 * PW);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_column), halfsize_column);
    }
};

class InitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit InitialCondition(SPHBody& sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][1] = -227.0;
    }
};

// define an observer body
class ColumnObserverParticleGenerator : public ObserverParticleGenerator
{
public:
    explicit ColumnObserverParticleGenerator(SPHBody& sph_body) : ObserverParticleGenerator(sph_body)
    {
        positions_.push_back(Vecd(0.0, PW));
        positions_.push_back(Vecd(PL, 0.0));
    }
};

int main(int ac, char* av[])
{
    /** Setup the system. Please the make sure the global domain bounds are correctly defined. */
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
#ifdef BOOST_AVAILABLE
    system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(system);

    RealBody column(system, makeShared<Column>("Column"));
    column.defineParticlesAndMaterial<J2PlasticicityParticles, J2Plasticity>(rho0_s, c0, Youngs_modulus, poisson, yield_stress);
    column.generateParticles<ParticleGeneratorLattice>();
    column.addBodyStateForRecording<Real>("VonMisesStress");
    column.addBodyStateForRecording<int>("PlasticIndicator");
    column.addBodyStateForRecording<Real>("Pressure");
    column.addBodyStateForRecording<Real>("Density");
    column.addBodyStateForRecording<Mat3d>("ShearStrain3D");
    column.addBodyStateForRecording<Mat3d>("ShearStress3D");

    SolidBody wall(system, makeShared<Wall>("Wall"));
    wall.defineParticlesAndMaterial<SolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall.generateParticles<ParticleGeneratorLattice>();

    /** Define Observer. */
    ObserverBody my_observer(system, "MyObserver");
    my_observer.generateParticles<ColumnObserverParticleGenerator>();

    /**body relation topology */
    InnerRelation column_inner(column);
    ContactRelation my_observer_contact(my_observer, { &column });
    SurfaceContactRelation column_wall_contact(column, { &wall });
    /**define simple data file input and outputs functions. */
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    //----------------------------------------------------------------------
    //	All numerical methods will be used in this case.
    //----------------------------------------------------------------------
    SimpleDynamics<InitialCondition> initial_condition(column);
    InteractionWithUpdate<KernelCorrectionMatrixInner> corrected_configuration(column_inner);
    Dynamics1Level<continuum_dynamics::Integration1stHalf> column_pressure_relaxation(column_inner);
    Dynamics1Level<continuum_dynamics::ShearStressRelaxationHourglassControlJ2Plasticity> column_shear_stress_relaxation(column_inner, 1);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerDissipativeRiemann> column_density_relaxation(column_inner);
    ReduceDynamics<continuum_dynamics::ContinuumAcousticTimeStepSize> computing_time_step_size(column, 0.3);
    InteractionDynamics<continuum_dynamics::DynamicContactForceWithWall> column_wall_contact_force(column_wall_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
    //----------------------------------------------------------------------
    //	Output
    //----------------------------------------------------------------------
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_velocity("Velocity", io_environment, my_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_displacement("Position", io_environment, my_observer_contact);

    //----------------------------------------------------------------------
    // From here the time stepping begins.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_normal_direction.exec();
    corrected_configuration.exec();
    initial_condition.exec();
    //----------------------------------------------------------------------
    // Setup time-stepping related simulation parameters.
    //----------------------------------------------------------------------
    int ite = 0;
    Real end_time = 1.0e-4;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real output_period = end_time / 100;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_displacement.writeToFile(0);
    write_velocity.writeToFile(0);
    //----------------------------------------------------------------------
    // Main time-stepping loop.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            Real dt = computing_time_step_size.exec();

            if (ite % screen_output_interval == 0)
            {
                std::cout << "N=" << ite << " Time: "
                    << GlobalStaticVariables::physical_time_ << "	dt: "
                    << dt << "\n";

                if (ite != 0 && ite % observation_sample_interval == 0)
                {
                    write_displacement.writeToFile(ite);
                    write_velocity.writeToFile(ite);
                }
            }
            column_wall_contact_force.exec(dt);

            column_pressure_relaxation.exec(dt);
            column_shear_stress_relaxation.exec(dt);
            column_density_relaxation.exec(dt);

            column.updateCellLinkedList();
            column_inner.updateConfiguration();
            column_wall_contact.updateConfiguration();
            corrected_configuration.exec();

            ite++;
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;

            //write_states.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}