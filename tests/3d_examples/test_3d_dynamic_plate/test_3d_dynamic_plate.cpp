/**
 * @file 	test_3d_dynamic_plate.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider deformation of a square plate under stepping load.
 * @author 	Dong Wu and Xiangyu Hu
 * @ref		A New Continuum-Based Thick Shell Finite Element for Soft Biological Tissues in Dynamics:
 *			Part 1 - Preliminary Benchmarking Using Classic Verification Experiments
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 0.254;                                        /** Length of the square plate. */
Real PH = 0.254;                                        /** Width of the square plate. */
Real PT = 0.0127;                                       /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                       /** Pseudo-normal. */
int particle_number = 40;                               /** Particle number in the direction of the length */
Real particle_spacing_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-BW, -BW, -0.5 * particle_spacing_ref),
                                 Vec3d(PL + BW, PH + BW, 0.5 * particle_spacing_ref));
// Observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.5 * PH, 0.0), Vecd(-BW, -BW, 0.0)};
/** For material properties of the solid. */
Real rho0_s = 1.0;                /** Normalized density. */
Real Youngs_modulus = 68.94e9;    /** Normalized Youngs Modulus. */
Real poisson = 0.3;               /** Poisson ratio. */
Real physical_viscosity = 2000.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 2068427.0; /** Total distributed load. */
Real time_to_full_external_force = 0.0;

Real gravitational_acceleration = 0.0;

/** Define application dependent particle generator for thin structure. */
namespace SPH
{
class Plate;
template <>
class ParticleGenerator<SurfaceParticles, Plate> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    virtual void prepareGeometricData() override
    {
        // the plate and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {
            for (int j = 0; j < (particle_number + 2 * BWD); j++)
            {
                Real x = particle_spacing_ref * i - BW + particle_spacing_ref * 0.5;
                Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
                addPositionAndVolumetricMeasure(Vecd(x, y, 0), particle_spacing_ref * particle_spacing_ref);
                addSurfaceProperties(n_0, PT);
            }
        }
    }
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        return base_particles_.ParticlePositions()[index_i][0] < 0.0 ||
               base_particles_.ParticlePositions()[index_i][1] < 0.0 ||
               base_particles_.ParticlePositions()[index_i][0] > PL ||
               base_particles_.ParticlePositions()[index_i][1] > PH;
    };
};
} // namespace SPH

/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    /** create a plate body. */
    SolidBody plate_body(sph_system, makeShared<DefaultShape>("PlateBody"));
    plate_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    plate_body.generateParticles<SurfaceParticles, Plate>();

    /** Define Observer. */
    ObserverBody plate_observer(sph_system, "PlateObserver");
    plate_observer.generateParticles<ObserverParticles>(observation_location);

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation plate_body_inner(plate_body);
    ContactRelation plate_observer_contact(plate_observer, {&plate_body});

    IncreaseToFullGravity time_dependent_external_force(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration), time_to_full_external_force);
    SimpleDynamics<GravityForce<IncreaseToFullGravity>> apply_time_dependent_external_force(plate_body, time_dependent_external_force);

    /**
     * This section define all numerical methods will be used in this case.
     */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(plate_body_inner);

    /** active-passive stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(plate_body_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(plate_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(plate_body);
    /** Constrain the Boundary. */
    BoundaryGeometry boundary_geometry(plate_body);
    SimpleDynamics<FixBodyPartConstraint> constrain_holder(boundary_geometry);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vec3d>(plate_body, "ForcePrior");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_plate_max_displacement("Position", plate_observer_contact);

    /** Apply initial condition. */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    write_states.writeToFile(0);
    write_plate_max_displacement.writeToFile(0);

    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 4.0e-5;
    // Real end_time = 0.8;
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
        Real integral_time = 0.0;
        while (integral_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
            }
            apply_time_dependent_external_force.exec();
            stress_relaxation_first_half.exec(dt);
            constrain_holder.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            physical_time += dt;
        }
        write_plate_max_displacement.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_plate_max_displacement.generateDataBase(0.005);
    }
    else
    {
        write_plate_max_displacement.testResult();
    }

    return 0;
}
