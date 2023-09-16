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
Real resolution_ref = L / 10.0; /**< reference particle spacing. */
Real BW = resolution_ref * 4;   /**< wall width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(25, 15));
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
        cubic_shape.push_back(Vecd(BW, slop_h + resolution_ref));
        cubic_shape.push_back(Vecd(BW, slop_h + L + resolution_ref));
        cubic_shape.push_back(Vecd(BW + L, slop_h + L + resolution_ref));
        cubic_shape.push_back(Vecd(BW + L, slop_h + resolution_ref));
        cubic_shape.push_back(Vecd(BW, slop_h + resolution_ref));
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles
    //----------------------------------------------------------------------
    SolidBody free_cube(sph_system, makeShared<Cube>("FreeCube"));
    free_cube.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    free_cube.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody cube_observer(sph_system, "CubeObserver");
    cube_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
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
    SimpleDynamics<TimeStepInitialization> free_cube_initialize_timestep(free_cube, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    /** Kernel correction. */
    InteractionWithUpdate<CorrectedConfigurationInner> free_cube_corrected_configuration(free_cube_inner);
    /** Time step size. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> free_cube_get_time_step_size(free_cube);
    /** stress relaxation for the solid body. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> free_cube_stress_relaxation_first_half(free_cube_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> free_cube_stress_relaxation_second_half(free_cube_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> free_cube_update_contact_density(free_cube_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> free_cube_compute_solid_contact_forces(free_cube_contact);
    /** Damping*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        damping(0.5, free_cube_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    /** Observer and output. */
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_free_cube_displacement("Position", io_environment, cube_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    GlobalStaticVariables::physical_time_ = 0.0;
    wall_boundary_rotation.exec();
    free_cube_rotation.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    free_cube_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_free_cube_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                free_cube_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt
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
                GlobalStaticVariables::physical_time_ += dt;
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

    if (sph_system.generate_regression_data_)
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
