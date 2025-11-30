/**
 * @file 	diffusion_NeumannBC.cpp
 * @brief 	2D test of diffusion problem with Neumann boundary condition.
 * @details This is the first case to validate multiple boundary conditions.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 100.0;
Real BW = resolution_ref * 2.0;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for diffusion properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
std::string diffusion_species_name = "Phi";
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 100.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real heat_flux = 900.0; // from the Nemann boundary
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

std::vector<Vecd> heat_flux_region_edge_points{
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
        multi_polygon_.addAPolygon(thermal_domain_edge_points, ShapeBooleanOps::add);
    }
};

class DirichletWallBoundary : public MultiPolygonShape
{
  public:
    explicit DirichletWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(left_region_edge_points, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(right_region_edge_points, ShapeBooleanOps::add);
    }
};

class NeumannWallBoundary : public MultiPolygonShape
{
  public:
    explicit NeumannWallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(heat_flux_region_edge_points, ShapeBooleanOps::add);
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
        Vec2d point_coordinate(0.5 * L, range_of_measure * Real(i) /
                                                Real(number_of_observation_points - 1) +
                                            start_of_measure);
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    diffusion_body.defineMaterial<Solid>();
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_Dirichlet(sph_system, makeShared<DirichletWallBoundary>("DirichletWallBoundary"));
    wall_Dirichlet.defineMaterial<Solid>();
    wall_Dirichlet.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_Neumann(sph_system, makeShared<NeumannWallBoundary>("NeumannWallBoundary"));
    wall_Neumann.defineMaterial<Solid>();
    wall_Neumann.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    MultiPolygonShape left_region_shape(MultiPolygon(left_region_edge_points), "LeftRegion");
    BodyRegionByParticle wall_Dirichlet_left_region(wall_Dirichlet, left_region_shape);

    MultiPolygonShape right_region_shape(MultiPolygon(right_region_edge_points), "RightRegion");
    BodyRegionByParticle wall_Dirichlet_right_region(wall_Dirichlet, right_region_shape);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> diffusion_body_contact_Dirichlet(diffusion_body, {&wall_Dirichlet});
    Contact<> diffusion_body_contact_Neumann(diffusion_body, {&wall_Neumann});
    Contact<> temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the configuration dynamics, such as update cell linked list,
    // update body relations, are defiend first.
    // Then the geometric models or simple objects without data dependencies,
    // such as gravity, initialized normal direction.
    // After that, the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> diffusion_body_cell_linked_list(diffusion_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_Dirichlet_cell_linked_list(wall_Dirichlet);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_Neumann_cell_linked_list(wall_Neumann);

    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>, Contact<>>
        water_block_update_complex_relation(
            diffusion_body_inner, diffusion_body_contact_Dirichlet, diffusion_body_contact_Neumann);
    UpdateRelation<MainExecutionPolicy, Contact<>> observer_contact_relation(temperature_observer_contact);

    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_Neumann);
    StateDynamics<MainExecutionPolicy, VariableAssignment<ConstantValue<Real>, SPHBody>>
        diffusion_initial_condition(diffusion_body, diffusion_species_name, initial_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<ConstantValue<Real>, BodyRegionByParticle>>
        left_initial_condition(wall_Dirichlet_left_region, diffusion_species_name, left_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<ConstantValue<Real>, BodyRegionByParticle>>
        right_initial_condition(wall_Dirichlet_right_region, diffusion_species_name, right_temperature);
    StateDynamics<MainExecutionPolicy, VariableAssignment<ConstantValue<Real>, SPHBody>>
        wall_Neumann_initial_condition(wall_Neumann, diffusion_species_name + "Flux", heat_flux);

    IsotropicDiffusion isotropic_diffusion(diffusion_species_name, diffusion_coeff);
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, &isotropic_diffusion);
    RungeKuttaSequence<InteractionDynamicsCK<
        MainExecutionPolicy,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta1stStage, IsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Neumann<IsotropicDiffusion>, NoKernelCorrectionCK>>,
        DiffusionRelaxationCK<
            Inner<OneLevel, RungeKutta2ndStage, IsotropicDiffusion, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Dirichlet<IsotropicDiffusion>, NoKernelCorrectionCK>,
            Contact<InteractionOnly, Neumann<IsotropicDiffusion>, NoKernelCorrectionCK>>>>
        diffusion_relaxation_rk2(DynamicsArgs(diffusion_body_inner, &isotropic_diffusion),
                                 DynamicsArgs(diffusion_body_contact_Dirichlet, &isotropic_diffusion),
                                 DynamicsArgs(diffusion_body_contact_Neumann, &isotropic_diffusion));
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> write_states(sph_system);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<MainExecutionPolicy, Real, RestoringCorrection>>
        write_solid_temperature(diffusion_species_name, temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();

    diffusion_body_cell_linked_list.exec();
    wall_Dirichlet_cell_linked_list.exec();
    wall_Neumann_cell_linked_list.exec();

    water_block_update_complex_relation.exec();
    observer_contact_relation.exec();

    diffusion_initial_condition.exec();
    left_initial_condition.exec();
    right_initial_condition.exec();
    wall_Neumann_initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
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
    write_solid_temperature.writeToFile(ite);
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

                diffusion_relaxation_rk2.exec(dt);

                ite++;
                dt = get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                sv_physical_time->incrementValue(dt);
            }
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        write_solid_temperature.writeToFile(ite);
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
        write_solid_temperature.generateDataBase(1.0e-3, 1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_solid_temperature.testResult();
    }

    return 0;
}