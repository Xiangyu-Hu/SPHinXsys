/**
 * @file 	test_3d_slender_beam.cpp
 * @brief 	This is the benchmark test of the bar.
 * @details  We consider the body force applied on a quasi-static straight beam.
 * @author 	Xipeng Lyu
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
Real PW = 1.0;                                    /** Thickness of the square plate. */
Vec3d n_0 = Vec3d(0.0, 0.0, 1.0);                 /** Pseudo-normal. */
Vec3d b_n_0 = Vec3d(0.0, 1.0, 0.0);               /** Pseudo-binormal. */
int particle_number = 40;                         /** Particle number in the direction of the length */
Real resolution_ref = PL / (Real)particle_number; /** Initial reference particle spacing. */
int BWD = 1;                                      /** Width of the boundary layer measured by number of particles. */
Real BW = resolution_ref * (Real)BWD;             /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-BW, -BW, -0.5 * resolution_ref),
                                 Vec3d(PL + BW, PH + BW, 0.5 * resolution_ref));
// Observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * PL, 0.0, 0.0)};

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
TEST(Beam, MaxDisplacement)
{
    Real displ_max = observed_quantity_n - observed_quantity_0;
    EXPECT_NEAR(displ_max, displ_max_reference, displ_max_reference * 0.1);
    std::cout << "displ_max: " << displ_max << std::endl;
}

namespace SPH
{
/** Define application dependent particle generator for thin structure. */
class Bar;
template <>
class ParticleGenerator<LinearParticles, Bar> : public ParticleGenerator<LinearParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, LinearParticles &linear_particles)
        : ParticleGenerator<LinearParticles>(sph_body, linear_particles)
    {
        sph_body.getSPHAdaptation().getKernel()->reduceOnce();
    };
    virtual void prepareGeometricData() override
    {
        // the beam and boundary
        for (int i = 0; i < (particle_number + 2 * BWD); i++)
        {

            Real x = resolution_ref * i - BW + resolution_ref * 0.5;
            Real y = 0.0;
            Real z = 0.0;
            addPositionAndVolumetricMeasure(Vecd(x, y, z), resolution_ref);
            addLineProperties(n_0, b_n_0, PT, PW);
        }
    }
};
/** Define the boundary geometry. */
class BoundaryGeometryParallelToXAxis : public BodyPartByParticle
{
  public:
    BoundaryGeometryParallelToXAxis(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method =
            std::bind(&BoundaryGeometryParallelToXAxis::tagManually, this, _1);
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);

    /** create a bar body. */
    SolidBody bar_body(sph_system, makeShared<DefaultShape>("BarBody"));
    bar_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    bar_body.generateParticles<LinearParticles, Bar>();

    /** Define Observer. */
    ObserverBody bar_observer(sph_system, "BarObserver");
    bar_observer.generateParticles<ObserverParticles>(observation_location);

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation bar_body_inner(bar_body);
    ContactRelation bar_observer_contact(bar_observer, {&bar_body});

    /** Common particle dynamics. */
    IncreaseToFullGravity time_dependent_external_force(Vec3d(0.0, 0.0, q / (PT * rho0_s) - gravitational_acceleration), time_to_full_external_force);
    SimpleDynamics<GravityForce<IncreaseToFullGravity>> apply_time_dependent_external_force(bar_body, time_dependent_external_force);
    InteractionDynamics<slender_structure_dynamics::BarCorrectConfiguration> corrected_configuration(bar_body_inner);
    /**
     * This section define all numerical methods will be used in this case.
     */
    /** active-passive stress relaxation. */
    Dynamics1Level<slender_structure_dynamics::BarStressRelaxationFirstHalf> stress_relaxation_first_half(bar_body_inner);
    Dynamics1Level<slender_structure_dynamics::BarStressRelaxationSecondHalf> stress_relaxation_second_half(bar_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<slender_structure_dynamics::BarAcousticTimeStepSize> computing_time_step_size(bar_body);
    /** Constrain the Boundary. */
    BoundaryGeometryParallelToXAxis boundary_geometry_x(bar_body);
    SimpleDynamics<slender_structure_dynamics::ConstrainBarBodyRegionAlongAxis>
        constrain_holder_x(boundary_geometry_x, 0);
    BoundaryGeometryParallelToYAxis boundary_geometry_y(bar_body);
    SimpleDynamics<slender_structure_dynamics::ConstrainBarBodyRegionAlongAxis>
        constrain_holder_y(boundary_geometry_y, 1);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        bar_position_damping(0.5, bar_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        bar_rotation_damping(0.5, bar_body_inner, "AngularVelocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        bar_rotation_b_damping(0.5, bar_body_inner, "BinormalAngularVelocity", physical_viscosity);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    ObservedQuantityRecording<Vecd> write_beam_max_displacement("Position", bar_observer_contact);

    /** Apply initial condition. */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();

    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    write_states.writeToFile(0);
    write_beam_max_displacement.writeToFile(0);
    observed_quantity_0 = write_beam_max_displacement.getObservedQuantity()[0][2];

    /** Setup physical parameters. */
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
            bar_position_damping.exec(dt);
            bar_rotation_damping.exec(dt);
            bar_rotation_b_damping.exec(dt);
            constrain_holder_x.exec(dt);
            constrain_holder_y.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            physical_time += dt;
        }
        write_beam_max_displacement.writeToFile(ite);
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    observed_quantity_n = write_beam_max_displacement.getObservedQuantity()[0][2];

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
