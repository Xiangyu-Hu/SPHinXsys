/**
 * @file 	sliding.cpp
 * @brief 	a 2D elastic cube slides on a rigid slope.
 * @details This is the a case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 20.0; /**< box length. */
Real DH = 13.0; /**< box height. */
Real L = 1.0;
Real slop_h = 11.55;
Real global_resolution = L / 10.0; /**< reference particle spacing. */
Real BW = global_resolution * 4;   /**< wall width for BCs. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(25, 15));
// Observer location
StdVec<Vecd> observation_location = {Vecd(7.2, 9.8)};
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e5;
Real poisson = 0.45;
Real gravity_g = 9.8;
Real physical_viscosity = 1000000.0;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> wall_shape{Vecd(0, 0), Vecd(0, slop_h), Vecd(DL, slop_h), Vecd(0, 0)};
        multi_polygon_.addAPolygon(wall_shape, ShapeBooleanOps::add);
    }
};

class Cube : public MultiPolygonShape
{
  public:
    explicit Cube(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> cubic_shape;
        cubic_shape.push_back(Vecd(BW, slop_h + global_resolution));
        cubic_shape.push_back(Vecd(BW, slop_h + L + global_resolution));
        cubic_shape.push_back(Vecd(BW + L, slop_h + L + global_resolution));
        cubic_shape.push_back(Vecd(BW + L, slop_h + global_resolution));
        cubic_shape.push_back(Vecd(BW, slop_h + global_resolution));
        multi_polygon_.addAPolygon(cubic_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles
    //----------------------------------------------------------------------
    SolidBody free_cube(sph_system, makeShared<Cube>("FreeCube"));
    free_cube.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    free_cube.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<BaseParticles, Lattice>();

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
    SurfaceContactRelation free_cube_contact(free_cube, {&wall_boundary});
    ContactRelation cube_observer_contact(cube_observer, {&free_cube});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    Transform transform2d(Rotation2d(-0.5235));
    SimpleDynamics<TranslationAndRotation> wall_boundary_rotation(wall_boundary, transform2d);
    SimpleDynamics<TranslationAndRotation> free_cube_rotation(free_cube, transform2d);
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(free_cube, gravity);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> free_cube_corrected_configuration(free_cube_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> free_cube_stress_relaxation_first_half(free_cube_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> free_cube_stress_relaxation_second_half(free_cube_inner);
    InteractionDynamics<solid_dynamics::ContactFactorSummation> free_cube_update_contact_density(free_cube_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> free_cube_compute_solid_contact_forces(free_cube_contact);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>> damping(0.5, free_cube_inner, "Velocity", physical_viscosity);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> free_cube_get_time_step_size(free_cube);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    /** Observer and output. */
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_free_cube_displacement("Position", cube_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_rotation.exec();
    free_cube_rotation.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    free_cube_corrected_configuration.exec();
    constant_gravity.exec();
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
    Real T0 = 2.5;
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

    if (sph_system.GenerateRegressionData())
    {
        // The lift force at the cylinder is very small and not important in this case.
        write_free_cube_displacement.generateDataBase({1.0e-2, 1.0e-2}, {1.0e-2, 1.0e-2});
    }
    else
    {
        write_free_cube_displacement.testResult();
    }

    return 0;
}
