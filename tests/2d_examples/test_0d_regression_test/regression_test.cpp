/**
 * @file egression_test.cpp
 * @brief This is a test case based on diffusion, which can be used to
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
Real global_resolution = H / 40.0;
Real BW = global_resolution * 4;
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
//----------------------------------------------------------------------
// Define extra classes which are used in the main program.
// These classes are defined under the namespace of SPH.
//----------------------------------------------------------------------
namespace SPH
{
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][0] >= 0 && pos_[index_i][0] <= L && pos_[index_i][1] >= 0 && pos_[index_i][1] <= H)
        {
            phi_[index_i] = initial_temperature;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
using DiffusionBodyRelaxation =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<CorrectedKernelGradientInner>, BaseDiffusion>>;

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

} // namespace SPH
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
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
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation diffusion_body_inner_relation(diffusion_body);
    ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> correct_configuration(diffusion_body_inner_relation);

    DiffusionBodyRelaxation diffusion_relaxation(diffusion_body_inner_relation);

    GetDiffusionTimeStepSize get_time_step_size(diffusion_body);
    BodyRegionByParticle left_boundary(diffusion_body, makeShared<MultiPolygonShape>(createLeftSideBoundary()));
    SimpleDynamics<ConstantConstraint<BodyRegionByParticle, Real>> left_boundary_condition(left_boundary, diffusion_species_name, high_temperature);
    BodyRegionByParticle other_boundary(diffusion_body, makeShared<MultiPolygonShape>(createOtherSideBoundary()));
    SimpleDynamics<ConstantConstraint<BodyRegionByParticle, Real>> other_boundary_condition(other_boundary, diffusion_species_name, low_temperature);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations of the simulation.
    //	Regression tests are also defined here.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
        write_solid_temperature(diffusion_species_name, temperature_observer_contact);
    BodyRegionByParticle inner_domain(diffusion_body, makeShared<MultiPolygonShape>(createInnerDomain(), "InnerDomain"));
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<Average<QuantitySummation<Real, BodyPartByParticle>>>>
        write_solid_average_temperature_part(inner_domain, diffusion_species_name);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correct_configuration.exec();
    setup_diffusion_initial_condition.exec();
    left_boundary_condition.exec();
    other_boundary_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real T0 = 20.0;
    Real end_time = T0;
    Real Output_Time = 0.1 * end_time;
    Real Observe_time = 0.1 * Output_Time;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states.writeToFile(0);
    write_solid_temperature.writeToFile(ite);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < Output_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
                if (ite % 1 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                diffusion_relaxation.exec(dt);
                left_boundary_condition.exec();
                other_boundary_condition.exec();
                ite++;
                dt = get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                if (ite % 100 == 0)
                {
                    write_solid_temperature.writeToFile(ite);
                    write_solid_average_temperature_part.writeToFile(ite);
                    write_states.writeToFile(ite);
                }
            }
        }

        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
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