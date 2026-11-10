/**
 * @file 	diffusion_RobinBC.cpp
 * @brief 	2D test of diffusion problem with Robin boundary condition.
 * @details This is the first case to implement single-value variables registered on particles.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real global_resolution = H / 100.0;
Real BW = global_resolution * 2.0;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for diffusion properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
std::string species_name = "Phi";
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 100.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real convection = 100.0;
Real T_infinity = 400.0;
//----------------------------------------------------------------------
//	Generate 2D geometrics used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> thermal_domain_edge_points{
    Vecd(0.0, 0.0), Vecd(0.0, H),
    Vecd(L, H), Vecd(L, 0.0), Vecd(0.0, 0.0)};

std::vector<Vecd> left_region_edge_points{
    Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW), Vecd(0.4 * L, H + BW),
    Vecd(0.4 * L, H), Vecd(0.3 * L, H)};

std::vector<Vecd> right_region_edge_points{
    Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW), Vecd(0.7 * L, H + BW),
    Vecd(0.7 * L, H), Vecd(0.6 * L, H)};

std::vector<Vecd> convection_region_edge_points{
    Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0), Vecd(0.55 * L, 0),
    Vecd(0.55 * L, -BW), Vecd(0.45 * L, -BW)};
//----------------------------------------------------------------------
//	Define Shapes used for SPH bodies and body-parts.
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(thermal_domain_edge_points, GeometricOps::add);
    }
};

class DirichletWallBoundary : public MultiPolygonShape
{
  public:
    explicit DirichletWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(left_region_edge_points, GeometricOps::add);
        multi_polygon_.addPolygon(right_region_edge_points, GeometricOps::add);
    }
};

class RobinWallBoundary : public MultiPolygonShape
{
  public:
    explicit RobinWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addPolygon(convection_region_edge_points, GeometricOps::add);
    }
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the middle line. */
    size_t number_of_observation_points = 5;
    Real range_of_measure = L;
    Real start_of_measure = 0;

    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(
            0.5 * L, range_of_measure * Real(i) / Real(number_of_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    diffusion_body.defineMatterMaterial<Solid>();
    diffusion_body.addMaterialProperty<IsotropicDiffusion>(species_name, diffusion_coeff);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary_Dirichlet(sph_system, makeShared<DirichletWallBoundary>("DirichletWallBoundary"));
    wall_boundary_Dirichlet.defineMatterMaterial<Solid>();
    wall_boundary_Dirichlet.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    MultiPolygonShape left_region_shape(MultiPolygon(left_region_edge_points), "LeftRegion");
    BodyRegionByParticle wall_Dirichlet_left_region(wall_boundary_Dirichlet, left_region_shape);

    MultiPolygonShape right_region_shape(MultiPolygon(right_region_edge_points), "RightRegion");
    BodyRegionByParticle wall_Dirichlet_right_region(wall_boundary_Dirichlet, right_region_shape);

    SolidBody wall_boundary_Robin(sph_system, makeShared<RobinWallBoundary>("RobinWallBoundary"));
    wall_boundary_Robin.defineMatterMaterial<Solid>();
    wall_boundary_Robin.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of temperature observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    auto &diffusion_inner = sph_system.addInnerRelation(diffusion_body);
    auto &contact_Dirichlet = sph_system.addContactRelation(diffusion_body, wall_boundary_Dirichlet);
    auto &contact_Robin = sph_system.addContactRelation(diffusion_body, wall_boundary_Robin);
    auto &temperature_observer_contact = sph_system.addContactRelation(temperature_observer, diffusion_body);
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    // Generally, the host methods should be able to run immediately.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &host_methods = sph_solver.getHostMethodContainer();
    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary_Robin).exec(); // host method able to run immediately
    auto &main_methods = sph_solver.getMainMethodContainer();
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defined first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    ParticleDynamicsGroup static_configuration;
    static_configuration.add(&main_methods.addCellLinkedListDynamics(diffusion_body));
    static_configuration.add(&main_methods.addCellLinkedListDynamics(wall_boundary_Dirichlet));
    static_configuration.add(&main_methods.addCellLinkedListDynamics(wall_boundary_Robin));
    static_configuration.add(&main_methods.addRelationDynamics(diffusion_inner, contact_Dirichlet, contact_Robin));
    static_configuration.add(&main_methods.addRelationDynamics(temperature_observer_contact));

    IsotropicDiffusion isotropic_diffusion(species_name, diffusion_coeff);
    auto &thermal_time_step = main_methods.addReturnDynamics<GetDiffusionTimeStepSize>(diffusion_body, &isotropic_diffusion);

    auto &runge_kutta_1st_stage =
        main_methods.addInteractionDynamicsOneLevel<
                        DiffusionRelaxationCK, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>(
                        diffusion_inner, &isotropic_diffusion)
            .addPostContactInteraction<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>(
                contact_Dirichlet, &isotropic_diffusion);//
//            .addPostContactInteraction<InteractionOnly, Robin<IsotropicDiffusion>, NoKernelCorrectionCK>(
//                contact_Robin, &isotropic_diffusion);

    auto &runge_kutta_2nd_stage =
        main_methods.addInteractionDynamicsOneLevel<
                        DiffusionRelaxationCK, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>(
                        diffusion_inner, &isotropic_diffusion)
            .addPostContactInteraction<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>(
                contact_Dirichlet, &isotropic_diffusion);
//            .addPostContactInteraction<InteractionOnly, Robin<IsotropicDiffusion>, NoKernelCorrectionCK>(
//                contact_Robin, &isotropic_diffusion);

    auto &runge_kutta_method = main_methods.addParticleDynamicsGroup();
    runge_kutta_method.add(&runge_kutta_1st_stage).add(&runge_kutta_2nd_stage);
    //----------------------------------------------------------------------
    //	Specify initial condition if necessary.
    //----------------------------------------------------------------------
    auto &diffusion_initial_condition = main_methods.addStateDynamics<
        VariableAssignment, ConstantValue<Real>>(diffusion_body, species_name, initial_temperature);
    auto &left_initial_condition = main_methods.addStateDynamics<
        VariableAssignment, ConstantValue<Real>>(wall_Dirichlet_left_region, species_name, left_temperature);
    auto &right_initial_condition = main_methods.addStateDynamics<
        VariableAssignment, ConstantValue<Real>>(wall_Dirichlet_right_region, species_name, right_temperature);
//    auto &robin_initial_condition = main_methods.addStateDynamics<
//        VariableAssignment, ConstantValue<Real>>(wall_boundary_Robin, species_name, T_infinity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    auto &write_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &observe_temperature = main_methods.addObserveRegression<
        RegressionTestEnsembleAverage, Real>(temperature_observer_contact, species_name);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    static_configuration.exec();

    diffusion_initial_condition.exec();
    left_initial_condition.exec();
    right_initial_condition.exec();
//    robin_initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingleVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 1;
    Real end_time = T0;
    Real Observe_time = 0.01 * end_time;
    Real Output_Time = 0.1 * end_time;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    observe_temperature.writeToFile(ite);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
                if (ite % 500 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << sv_physical_time->getValue() << "	dt: "
                              << dt << "\n";
                }

                runge_kutta_method.exec(dt);

                ite++;
                dt = thermal_time_step.exec();
                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
            }
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        observe_temperature.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TickCount::interval_t tt;
    tt = t4 - t1 - interval;

    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << sv_physical_time->getValue() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        observe_temperature.generateDataBase(1.0e-3, 1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        observe_temperature.testResult();
    }

    return 0;
}
