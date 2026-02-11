/**
 * @file 	test_3d_thin_plate.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider the body force applied on a quasi-static square plate.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.ijnonlinmec.2014.04.009, doi.org/10.1201/9780849384165
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real PL = 10.0;                                   /** Length of the square plate. */
Real PH = 10.0;                                   /** Width of the square plate. */
Real PT = 1.0;                                    /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                 /** Pseudo-normal. */
int particle_number = 40;                         /** Particle number in the direction of the length */
Real global_resolution = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;                                      /** Width of the boundary layer measured by number of particles. */
Real BW = global_resolution * (Real)BWD;             /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-BW, -BW, -0.5 * global_resolution),
                                 Vec3d(PL + BW, PH + BW, 0.5 * global_resolution));
// Observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.5 * PH, 0.0)};

/** For material properties of the solid. */
Real rho0_s = 1.0;                 /** Normalized density. */
Real Youngs_modulus = 1.3024653e6; /** Normalized Youngs Modulus. */
Real poisson = 0.3;                /** Poisson ratio. */
Real physical_viscosity = 200.0;   /** physical damping, here we choose the same value as numerical viscosity. */

Real q = 100.0 * Youngs_modulus * 1.0e-4; /** Total distributed load. */
Real time_to_full_external_force = 0.1;

Real gravitational_acceleration = 0.009646;

Real observed_quantity_0(0.0);
Real observed_quantity_n(0.0);
Real displ_max_reference = 1.8687;
TEST(Plate, MaxDisplacement)
{
    Real displ_max = observed_quantity_n - observed_quantity_0;
    EXPECT_NEAR(displ_max, displ_max_reference, displ_max_reference * 0.1);
    std::cout << "displ_max: " << displ_max << std::endl;
}

namespace SPH
{
/** Define application dependent particle generator for thin structure. */
class Plate;
template <>
class ParticleGenerator<SurfaceParticles, Plate> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles) : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    virtual void prepareGeometricData() override
    {
        // the plate and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {
            for (int j = 0; j < (particle_number + 2 * BWD); j++)
            {
                Real x = global_resolution * i - BW + global_resolution * 0.5;
                Real y = global_resolution * j - BW + global_resolution * 0.5;
                addPositionAndVolumetricMeasure(Vecd(x, y, 0.0), global_resolution * global_resolution);
                addSurfaceProperties(n_0, PT);
            }
        }
    }
};
/** Define the boundary geometry. */
class BoundaryGeometryParallelToXAxis : public BodyPartByParticle
{
  public:
    BoundaryGeometryParallelToXAxis(SPHBody &body)
        : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometryParallelToXAxis::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometryParallelToXAxis() {};

  private:
    bool tagManually(size_t index_i)
    {
        return base_particles_.ParticlePositions()[index_i][1] < 0.0 ||
               base_particles_.ParticlePositions()[index_i][1] > PH;
    };
};
class BoundaryGeometryParallelToYAxis : public BodyPartByParticle
{
  public:
    BoundaryGeometryParallelToYAxis(SPHBody &body)
        : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometryParallelToYAxis::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometryParallelToYAxis() {};

  private:
    bool tagManually(size_t index_i)
    {
        return base_particles_.ParticlePositions()[index_i][0] < 0.0 ||
               base_particles_.ParticlePositions()[index_i][0] > PL;
    };
};
} // namespace SPH

/**
 *  The main program
 */
int main(int ac, char *av[])
{
    /** Setup the system. */
    SPHSystem sph_system(system_domain_bounds, global_resolution);

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
    BoundaryGeometryParallelToXAxis boundary_geometry_x(plate_body);
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegionAlongAxis>
        constrain_holder_x(boundary_geometry_x, 0);
    BoundaryGeometryParallelToYAxis boundary_geometry_y(plate_body);
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegionAlongAxis>
        constrain_holder_y(boundary_geometry_y, 1);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        plate_position_damping(0.5, plate_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        plate_rotation_damping(0.5, plate_body_inner, "AngularVelocity", physical_viscosity);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    ObservedQuantityRecording<Vecd> write_plate_max_displacement("Position", plate_observer_contact);

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
    observed_quantity_0 = write_plate_max_displacement.getObservedQuantity()[0][2];

    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 0.8;
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
            constrain_holder_x.exec(dt);
            constrain_holder_y.exec(dt);
            plate_position_damping.exec(dt);
            plate_rotation_damping.exec(dt);
            constrain_holder_x.exec(dt);
            constrain_holder_y.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            physical_time += dt;
        }
        write_plate_max_displacement.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    observed_quantity_n = write_plate_max_displacement.getObservedQuantity()[0][2];

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
