/**
 * @file regression_test.cpp
 * @brief This is a test case based on diffusion and using computing kernel, which can be used to
 * validate the generation of the converged database in a regression test.
 * It can be run successfully (using CMake's CTest) in Linux system installed with Python 3.
 * @author 	Bo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
using namespace SPH;   // namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and simulation setup.
//----------------------------------------------------------------------
Real L = 0.2;
Real H = 0.2;
Real resolution_ref = H / 40.0;
Real BW = resolution_ref * 4;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 1.0e-3;
Real bias_coeff = 0.0;
Real alpha = Pi / 4.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
Real initial_temperature = 0.0;
Real high_temperature = 1.0;
Real low_temperature = 0.0;
//----------------------------------------------------------------------
//	Cases-dependent 2D geometries
//----------------------------------------------------------------------
MultiPolygon createDiffusionDomain()
{
    // thermal solid domain geometry.
    std::vector<Vecd> diffusion_domain;
    diffusion_domain.push_back(Vecd(-BW, -BW));
    diffusion_domain.push_back(Vecd(-BW, H + BW));
    diffusion_domain.push_back(Vecd(L + BW, H + BW));
    diffusion_domain.push_back(Vecd(L + BW, -BW));
    diffusion_domain.push_back(Vecd(-BW, -BW));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(diffusion_domain, ShapeBooleanOps::add);
    return multi_polygon;
}

MultiPolygon createInnerDomain()
{
    // thermal solid inner domain geometry.
    std::vector<Vecd> inner_domain;
    inner_domain.push_back(Vecd(0.0, 0.0));
    inner_domain.push_back(Vecd(0.0, H));
    inner_domain.push_back(Vecd(L, H));
    inner_domain.push_back(Vecd(L, 0.0));
    inner_domain.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(inner_domain, ShapeBooleanOps::add);

    return multi_polygon;
}

MultiPolygon createLeftSideBoundary()
{
    // left isothermal boundary geometry.
    std::vector<Vecd> left_boundary;
    left_boundary.push_back(Vecd(-BW, -BW));
    left_boundary.push_back(Vecd(-BW, H + BW));
    left_boundary.push_back(Vecd(0.0, H));
    left_boundary.push_back(Vecd(0.0, 0.0));
    left_boundary.push_back(Vecd(-BW, -BW));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(left_boundary, ShapeBooleanOps::add);

    return multi_polygon;
}

MultiPolygon createOtherSideBoundary()
{
    // other side isothermal boundary geometry.
    std::vector<Vecd> other_boundaries;
    other_boundaries.push_back(Vecd(-BW, -BW));
    other_boundaries.push_back(Vecd(0.0, 0.0));
    other_boundaries.push_back(Vecd(L, 0.0));
    other_boundaries.push_back(Vecd(L, H));
    other_boundaries.push_back(Vecd(0.0, L));
    other_boundaries.push_back(Vecd(-BW, H + BW));
    other_boundaries.push_back(Vecd(L + BW, H + BW));
    other_boundaries.push_back(Vecd(L + BW, -BW));
    other_boundaries.push_back(Vecd(-BW, -BW));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(other_boundaries, ShapeBooleanOps::add);

    return multi_polygon;
}

StdVec<Vecd> createObservationPoints()
{
    /** A line of measuring points at the middle line. */
    size_t number_of_observation_points = 11;
    Real range_of_measure = L - BW;
    Real start_of_measure = BW;

    StdVec<Vecd> observation_points;
    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(0.5 * L, start_of_measure + range_of_measure * (Real)i / (Real)(number_of_observation_points - 1));
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
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Create body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<MultiPolygonShape>(createDiffusionDomain(), "DiffusionBody"));
    diffusion_body.defineClosure<Solid, DirectionalDiffusion>(
        Solid(), ConstructArgs(diffusion_species_name, diffusion_coeff, bias_coeff, bias_direction));
    diffusion_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Observer body
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    BodyRegionByParticle left_boundary(diffusion_body, makeShared<MultiPolygonShape>(createLeftSideBoundary()));
    BodyRegionByParticle other_boundary(diffusion_body, makeShared<MultiPolygonShape>(createOtherSideBoundary()));
    BodyRegionByParticle inner_domain(diffusion_body, makeShared<MultiPolygonShape>(createInnerDomain(), "InnerDomain"));
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    Inner<> diffusion_body_inner(diffusion_body);
    Contact<> observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    // Define SPH solver with particle methods and execution policies.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
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
    auto &diffusion_body_cell_linked_list = main_methods.addCellLinkedListDynamics(diffusion_body);
    auto &water_block_update_complex_relation = main_methods.addRelationDynamics(diffusion_body_inner);
    auto &observer_update_contact_relation = main_methods.addRelationDynamics(observer_contact);

    auto &correct_configuration = main_methods.addInteractionDynamics<LinearCorrectionMatrixInner>(diffusion_body_inner);
    auto &diffusion_initial_condition = main_methods.addStateDynamics<VariableAssignment, ConstantValue<Real>>(
        diffusion_body, diffusion_species_name, initial_temperature);

    GetDiffusionTimeStepSize get_time_step_size(diffusion_body);
    auto &diffusion_relaxation_rk2 =
        main_methods.addRK2Sequence<DiffusionRelaxationCK, DirectionalDiffusion, LinearCorrectionCK>(diffusion_body_inner);

    auto &left_boundary_condition =
        main_methods.addStateDynamics<ConstantConstraintCK, Real>(left_boundary, diffusion_species_name, high_temperature);
    auto &other_boundary_condition =
        main_methods.addStateDynamics<ConstantConstraintCK, Real>(other_boundary, diffusion_species_name, low_temperature);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations of the simulation.
    //	Regression tests are also defined here.
    //----------------------------------------------------------------------
    auto &write_states = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_solid_temperature =
        main_methods.addObserveRegression<RegressionTestEnsembleAverage, Real>(
            diffusion_species_name, observer_contact);
    auto &write_solid_average_temperature_part =
        main_methods.addReduceRegression<RegressionTestDynamicTimeWarping, QuantityAverage, Real>(
            inner_domain, diffusion_species_name);
    //----------------------------------------------------------------------
    //	Define time stepper with end and start time.
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(20.0);
    size_t iteration_steps = 0;
    auto &state_recording = time_stepper.addTriggerByInterval(0.2);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    diffusion_body_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    observer_update_contact_relation.exec();

    correct_configuration.exec();
    diffusion_initial_condition.exec();
    left_boundary_condition.exec();
    other_boundary_condition.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile();
    write_solid_temperature.writeToFile(iteration_steps);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (!time_stepper.isEndTime())
    {
        Real dt = time_stepper.incrementPhysicalTime(get_time_step_size);
        diffusion_relaxation_rk2.exec(dt);
        left_boundary_condition.exec();
        other_boundary_condition.exec();
        iteration_steps++;

        if (iteration_steps % 10 == 0)
        {
            std::cout << "N=" << iteration_steps << " Time: " << time_stepper.getPhysicalTime()
                      << "	dt: " << time_stepper.getGlobalTimeStepSize() << "\n";
        }

        if (iteration_steps % 100 == 0)
        {
            write_solid_temperature.writeToFile(iteration_steps);
            write_solid_average_temperature_part.writeToFile(iteration_steps);
        }

        if (state_recording())
        {
            write_states.writeToFile();
        }
    }
    //----------------------------------------------------------------------
    //	@ensemble_average_method.
    //	The first argument is the threshold of mean value convergence.
    //	The second argument is the threshold of variance convergence.
    //----------------------------------------------------------------------
    write_solid_temperature.generateDataBase(0.001, 0.001);
    //----------------------------------------------------------------------
    //	@dynamic_time_warping_method.
    //	The value is the threshold of dynamic time warping (dtw) distance.
    //----------------------------------------------------------------------
    write_solid_average_temperature_part.generateDataBase(0.001);

    return 0;
}