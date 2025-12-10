/**
 * @file 	3d_roof.cpp
 * @brief 	This is the benchmark test of the shell.
 * @details  We consider the deformation of a cylindrical surface.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1007/s00466-017-1498-9, doi.org/10.1016/0045-7825(89)90098-4
 */
#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 24.875;                               /** Radius of the inner boundary of the cylinder. */
Real height = 50.0;                                 /** Height of the cylinder. */
Real thickness = 0.25;                              /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0; /** Radius of the mid surface. */
int particle_number = 16;                           /** Particle number in the peripheral direction. */
/** Initial reference particle spacing. */
Real particle_spacing_ref = 2.0 * radius_mid_surface * Pi * 80.0 / 360.0 / (Real)particle_number;
int BWD = 1;                                /** Width of the boundary layer measured by number of particles. */
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec3d(-radius - thickness, 0.0, -radius - thickness),
                                 Vec3d(radius + thickness + BW, height, radius + thickness));
// Observer location
StdVec<Vecd> observation_location = {Vecd(radius_mid_surface * cos((50.0 - 2.0 * 80.0 / particle_number) / 180.0 * Pi),
                                          0.5 * height,
                                          radius_mid_surface *sin((50.0 - 2.0 * 80.0 / particle_number) / 180.0 * Pi))};
/** For material properties of the solid. */
Real rho0_s = 36.0;                          /** Normalized density. */
Real Youngs_modulus = 4.32e8;                /** Normalized Youngs Modulus. */
Real poisson = 0.0;                          /** Poisson ratio. */
Real physical_viscosity = 7.0e3 * thickness; /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.1;
Real gravitational_acceleration = -10.0;

Real observed_quantity_0 = 0.0;
Real observed_quantity_n = 0.0;
Real displ_max_reference = 0.3024;
TEST(Plate, MaxDisplacement)
{
    Real displ_max = observed_quantity_0 - observed_quantity_n;
    EXPECT_NEAR(displ_max, displ_max_reference, displ_max_reference * 0.1);
    std::cout << "displ_max: " << displ_max << std::endl;
}

namespace SPH
{
/** Define application dependent particle generator for thin structure. */
class Cylinder;
template <>
class ParticleGenerator<SurfaceParticles, Cylinder> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    virtual void prepareGeometricData() override
    {
        // the cylinder and boundary
        for (int i = 0; i < particle_number + 1; i++)
        {
            for (int j = 0; j < (height / particle_spacing_ref + 2 * BWD); j++)
            {
                Real x = radius_mid_surface * cos(50.0 / 180.0 * Pi + i * 80.0 / 360.0 * 2 * Pi / (Real)particle_number);
                Real y = particle_spacing_ref * j - BW + particle_spacing_ref * 0.5;
                Real z = radius_mid_surface * sin(50.0 / 180.0 * Pi + i * 80.0 / 360.0 * 2 * Pi / (Real)particle_number);
                addPositionAndVolumetricMeasure(Vecd(x, y, z), particle_spacing_ref * particle_spacing_ref);
                Vecd n_0 = Vec3d(x / radius_mid_surface, 0.0, z / radius_mid_surface);
                addSurfaceProperties(n_0, thickness);
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
        return base_particles_.ParticlePositions()[index_i][1] < 0.0 ||
               base_particles_.ParticlePositions()[index_i][1] > height + 0.5 * particle_spacing_ref;
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
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    /** Create a Cylinder body. */
    SolidBody cylinder_body(sph_system, makeShared<DefaultShape>("CylinderBody"));
    cylinder_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    cylinder_body.generateParticles<SurfaceParticles, Cylinder>();
    /** Define Observer. */
    ObserverBody cylinder_observer(sph_system, "CylinderObserver");
    cylinder_observer.generateParticles<ObserverParticles>(observation_location);

    /** Set body contact map
     *  The contact map gives the data connections between the bodies
     *  basically the the range of bodies to build neighbor particle lists
     */
    InnerRelation cylinder_body_inner(cylinder_body);
    ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder_body});

    /** Common particle dynamics. */
    IncreaseToFullGravity time_dependent_external_force(Vec3d(0.0, 0.0, gravitational_acceleration), time_to_full_external_force);
    SimpleDynamics<GravityForce<IncreaseToFullGravity>> apply_time_dependent_external_force(cylinder_body, time_dependent_external_force);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(cylinder_body_inner);
    /**
     * This section define all numerical methods will be used in this case.
     */
    /** stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(cylinder_body_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(cylinder_body_inner);

    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(cylinder_body);
    BoundaryGeometry boundary_geometry(cylinder_body);
    SimpleDynamics<FixedInAxisDirection> constrain_holder(boundary_geometry, Vecd(0.0, 1.0, 0.0));
    DampingWithRandomChoice<InteractionSplit<DampingProjectionInner<Vec3d, FixedDampingRate>>>
        cylinder_position_damping(0.3, cylinder_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingProjectionInner<Vec3d, FixedDampingRate>>>
        cylinder_rotation_damping(0.3, cylinder_body_inner, "AngularVelocity", physical_viscosity);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_cylinder_max_displacement("Position", cylinder_observer_contact);

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
    write_cylinder_max_displacement.writeToFile(0);
    observed_quantity_0 = write_cylinder_max_displacement.getObservedQuantity()[0][2];

    /** Setup physical parameters. */
    int ite = 0;
    Real end_time = 2.0;
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
            dt = computing_time_step_size.exec();
            apply_time_dependent_external_force.exec();
            stress_relaxation_first_half.exec(dt);

            constrain_holder.exec();
            cylinder_position_damping.exec(dt);
            cylinder_rotation_damping.exec(dt);
            constrain_holder.exec();

            stress_relaxation_second_half.exec(dt);

            ite++;
            integral_time += dt;
            physical_time += dt;
        }
        write_cylinder_max_displacement.writeToFile(ite);
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
        write_cylinder_max_displacement.generateDataBase(0.005);
    }
    else
    {
        write_cylinder_max_displacement.testResult();
    }

    observed_quantity_n = write_cylinder_max_displacement.getObservedQuantity()[0][2];

    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
