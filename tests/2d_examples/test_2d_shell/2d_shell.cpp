/**
 * @file 	2d_shell.cpp
 * @brief 	This is the benchmark test of the 2d shell.
 * @details  We consider large deformation of a 2D shell.
 * @author 	Dong Wu, Chi Zhang and Xiangyu Hu
 * @version  0.1
 */
#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real radius = 24.5;                                            /** Radius of the inner boundary of the cylinder. */
Real thickness = 1.0;                                          /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;            /** Radius of the mid surface. */
int particle_number = 2;                                       /** Particle number in the thickness. */
Real particle_spacing_ref = thickness / (Real)particle_number; /** Initial reference particle spacing. */
int particle_number_mid_surface = 2 * radius_mid_surface * Pi * 80.0 / 360.0 / particle_spacing_ref;
int BWD = 1;                                /** number of boundary particles layers . */
Real BW = particle_spacing_ref * (Real)BWD; /** Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-radius - thickness, 0.0),
                                 Vec2d(radius + thickness, radius + thickness));
// Shell observer location
StdVec<Vecd> observation_location = {Vecd(0.0, radius_mid_surface)};
/** For material properties of the solid. */
Real rho0_s = 3.67346939;         /** Normalized density. */
Real Youngs_modulus = 4.32e7;     /** Normalized Youngs Modulus. */
Real poisson = 0.3;               /** Poisson ratio. */
Real physical_viscosity = 2000.0; /** physical damping, here we choose the same value as numerical viscosity. */

Real time_to_full_external_force = 0.1;
Real gravitational_acceleration = -10000.0;

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
        for (int i = 0; i < particle_number_mid_surface + 2 * BWD; i++)
        {
            Real x = radius_mid_surface *
                     cos(50.0 / 180.0 * Pi + (i + 0.5 - BWD) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
            Real y = radius_mid_surface *
                     sin(50.0 / 180.0 * Pi + (i + 0.5 - BWD) * 80.0 / 360.0 * 2 * Pi / (Real)particle_number_mid_surface);
            addPositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
            Vec2d normal_direction = Vec2d(x / radius_mid_surface, y / radius_mid_surface);
            addSurfaceProperties(normal_direction, thickness);
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
        return base_particles_.ParticlePositions()[index_i][0] < -radius_mid_surface * cos(50.0 / 180.0 * Pi) ||
               base_particles_.ParticlePositions()[index_i][0] > radius_mid_surface * cos(50.0 / 180.0 * Pi);
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
    /**
     * This section define all numerical methods will be used in this case.
     */
    IncreaseToFullGravity time_dependent_external_force(Vec2d(0.0, gravitational_acceleration), time_to_full_external_force);
    SimpleDynamics<GravityForce<IncreaseToFullGravity>> apply_time_dependent_external_force(cylinder_body, time_dependent_external_force);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> corrected_configuration(cylinder_body_inner);
    /** The main shell dynamics model: stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(cylinder_body_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(cylinder_body_inner);

    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(cylinder_body);
    BoundaryGeometry boundary_geometry(cylinder_body);
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> fixed_free_rotate_shell_boundary(boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
        cylinder_position_damping(0.2, cylinder_body_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
        cylinder_rotation_damping(0.2, cylinder_body_inner, "AngularVelocity", physical_viscosity);
    /** Output */
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(cylinder_body, "PseudoNormal");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_cylinder_max_displacement("Position", cylinder_observer_contact);

    /** Apply initial condition. */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();

    write_states.writeToFile(0);
    write_cylinder_max_displacement.writeToFile(0);

    /** Setup physical parameters. */
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
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
            fixed_free_rotate_shell_boundary.exec(dt);
            cylinder_position_damping.exec(dt);
            cylinder_rotation_damping.exec(dt);
            fixed_free_rotate_shell_boundary.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
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
        write_cylinder_max_displacement.generateDataBase(0.05);
    }
    else
    {
        write_cylinder_max_displacement.testResult();
    }

    return 0;
}
