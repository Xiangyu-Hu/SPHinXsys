/**
 * @file 	sliding.cpp
 * @brief 	a 2D elastic cube slides on a rigid slope.
 * @details This is the a case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
#include <gtest/gtest.h>
using namespace SPH; //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real L = 1.0;                     /**< box length. */
const Real DL = 15.0;                   /**< slope length. */
const Real angle = 10.0 * M_PI / 180.0; /**< angle of the slope. */
const Real Lx = DL * cos(angle);
const Real Ly = DL * sin(angle);
const Real resolution_ref = L / 20.0; /**< reference particle spacing. */
const Real resolution_shell = resolution_ref;
const Real thickness_shell = resolution_ref;
const Real BW = resolution_ref * 4; /**< wall width for BCs. */
/** Domain bounds of the system. */
const BoundingBoxd system_domain_bounds(Vec2d(-BW *cos(angle), -Ly - BW * sin(angle)), Vec2d(Lx + BW * cos(angle), L));
// Observer location at box center
const StdVec<Vecd> observation_location = {0.5 * L * Vecd(1.0, 1.0)};
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
const Real rho0_s = 1.0e3;
const Real Youngs_modulus = 1.0e5;
const Real poisson = 0.45;
const Real gravity_g = 9.8;
const Real physical_viscosity = 0.25 * sqrt(rho0_s * Youngs_modulus) * L;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class Cube : public MultiPolygonShape
{
  public:
    explicit Cube(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d halfsize_cube(0.5 * L, 0.5 * L);
        Transform translation_cube(halfsize_cube + 0.65 * (resolution_ref + resolution_shell) * Vec2d::UnitY());
        multi_polygon_.addABox(translation_cube, halfsize_cube, ShapeBooleanOps::add);
    }
};
namespace SPH
{
class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    void prepareGeometricData() override
    {
        Real s0 = -BW + 0.5 * resolution_shell;
        Real x = s0 * cos(angle);
        Real y = -s0 * sin(angle);

        while (x < DL + BW)
        {
            addPositionAndVolumetricMeasure(Vec2d(x, y), resolution_shell);
            addSurfaceProperties(Vec2d(-sin(angle), -cos(angle)), thickness_shell);
            x += resolution_shell * cos(angle);
            y -= resolution_shell * sin(angle);
        }
    }
};
} // namespace SPH
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
void run_simulation()
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles
    //----------------------------------------------------------------------
    SolidBody free_cube(sph_system, makeShared<Cube>("FreeCube"));
    free_cube.defineBodyLevelSetShape();
    free_cube.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    free_cube.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("Wall"));
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<SurfaceParticles, WallBoundary>();

    ObserverBody cube_observer(sph_system, "CubeObserver");
    cube_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation free_cube_inner(free_cube);
    SurfaceContactRelation free_cube_contact(free_cube, {&wall_boundary}, {false});
    ShellInnerRelationWithContactKernel wall_curvature_inner(wall_boundary, free_cube);
    ContactRelation cube_observer_contact(cube_observer, {&free_cube});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /**shell curvature*/
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> wall_average_curvature(wall_curvature_inner);
    Transform transform2d((Rotation2d(-angle)));
    SimpleDynamics<TranslationAndRotation> free_cube_rotation(free_cube, transform2d);
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(free_cube, gravity);
    /** Kernel correction. */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> free_cube_corrected_configuration(free_cube_inner);
    /** stress relaxation for the solid body. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> free_cube_stress_relaxation_first_half(free_cube_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> free_cube_stress_relaxation_second_half(free_cube_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactFactorSummation> free_cube_update_contact_density(free_cube_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> free_cube_compute_solid_contact_forces(free_cube_contact);
    /** Damping*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
        damping(0.5, free_cube_inner, "Velocity", physical_viscosity);
    /** Time step size. */
    ReduceDynamics<solid_dynamics::AcousticTimeStep> free_cube_get_time_step_size(free_cube);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    // analytical solution
    Real sin_theta = sin(angle);
    Real cos_theta = cos(angle);
    // analytical displacement under gravity
    auto get_analytical_displacement = [&](Real time)
    {
        Real a = 0.5 * gravity_g * sin_theta * time * time;
        Real u_x = a * cos_theta;
        Real u_y = a * sin_theta;
        return Vec2d(u_x, u_y);
    };
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    /** Observer and output. */
    ObservedQuantityRecording<Vecd>
        write_free_cube_displacement("Position", cube_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    free_cube_rotation.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    free_cube_corrected_configuration.exec();
    constant_gravity.exec();
    wall_average_curvature.exec();
    free_cube_contact.updateConfiguration();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_free_cube_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 4.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
    Real Dt = 0.1 * output_interval;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: " << dt
                              << "\n";
                }
                free_cube_update_contact_density.exec();
                free_cube_compute_solid_contact_forces.exec();

                free_cube_stress_relaxation_first_half.exec(dt);
                free_cube_stress_relaxation_second_half.exec(dt);

                free_cube.updateCellLinkedList();
                free_cube_contact.updateConfiguration();

                ite++;
                dt = free_cube_get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            write_free_cube_displacement.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    // gtest
    const Vec2d analytical_disp = get_analytical_displacement(physical_time);
    const Vec2d disp = write_free_cube_displacement.getObservedQuantity()[0] - observation_location[0];
    for (int n = 0; n < 2; n++)
        ASSERT_NEAR(abs(disp[n]), abs(analytical_disp[n]), 0.05 * analytical_disp.norm());
}

TEST(sliding_2d, displacement)
{
    run_simulation();
}

int main(int ac, char *av[])
{
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}