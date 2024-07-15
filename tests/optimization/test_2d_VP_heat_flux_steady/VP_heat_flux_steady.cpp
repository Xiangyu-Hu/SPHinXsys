/**
 * @file 	VP_heat_flux_steady.cpp
 * @brief 	This is the steady test for the heat flux problem.
 * @author 	Bo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
#include <gtest/gtest.h>
using namespace SPH; // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 50.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 0.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real heat_flux = 2000.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
    std::vector<Vecd> thermalDomainShape;
    thermalDomainShape.push_back(Vecd(0.0, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, H));
    thermalDomainShape.push_back(Vecd(L, H));
    thermalDomainShape.push_back(Vecd(L, 0.0));
    thermalDomainShape.push_back(Vecd(0.0, 0.0));

    return thermalDomainShape;
}

std::vector<Vecd> createBoundaryDomain()
{
    std::vector<Vecd> boundaryDomain;
    boundaryDomain.push_back(Vecd(-BW, -BW));
    boundaryDomain.push_back(Vecd(-BW, H + BW));
    boundaryDomain.push_back(Vecd(L + BW, H + BW));
    boundaryDomain.push_back(Vecd(L + BW, -BW));
    boundaryDomain.push_back(Vecd(-BW, -BW));

    return boundaryDomain;
}

std::vector<Vecd> heat_flux_boundary{
    Vecd(0.45 * L, H - resolution_ref), Vecd(0.45 * L, H), Vecd(0.55 * L, H),
    Vecd(0.55 * L, H - resolution_ref), Vecd(0.45 * L, H - resolution_ref)};
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
namespace SPH
{
class DiffusionBody : public MultiPolygonShape
{
  public:
    explicit DiffusionBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createBoundaryDomain(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::sub);
    }
};

MultiPolygon createHeatFluxBoundary()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(heat_flux_boundary, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DiffusionBodyInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit DiffusionBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          pos_(*particles_->getVariableDataByName<Vecd>("Position")),
          phi_(*particles_->registerSharedVariable<Real>("Phi")){};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = 550.0 + 50.0 * rand_uniform(0.0, 1.0);
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &phi_;
};

class WallBoundaryInitialCondition : public LocalDynamics, public DataDelegateSimple
{
  public:
    explicit WallBoundaryInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
          pos_(*particles_->getVariableDataByName<Vecd>("Position")),
          phi_(*particles_->registerSharedVariable<Real>("Phi")),
          heat_flux_(*(particles_->getVariableDataByName<Real>("HeatFlux"))) {}

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;
        if (pos_[index_i][1] < 0.0 && pos_[index_i][0] > 0.3 * L && pos_[index_i][0] < 0.4 * L)
        {
            phi_[index_i] = left_temperature;
        }
        if (pos_[index_i][1] < 0.0 && pos_[index_i][0] > 0.6 * L && pos_[index_i][0] < 0.7 * L)
        {
            phi_[index_i] = right_temperature;
        }
        if (pos_[index_i][1] > H && pos_[index_i][0] > 0.45 * L && pos_[index_i][0] < 0.55 * L)
        {
            heat_flux_[index_i] = heat_flux;
        }
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &phi_, &heat_flux_;
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the middle line. */
    size_t number_of_observation_points = 11;
    Real range_of_measure = L;
    Real start_of_measure = 0;

    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(0.5 * L, range_of_measure * Real(i) / Real(number_of_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
} // namespace SPH

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
TEST(test_optimization, test_problem4_non_optimization)
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
    LocalIsotropicDiffusion *local_isotropic_diffusion =
        diffusion_body.defineMaterial<LocalIsotropicDiffusion>("Phi", "Phi", diffusion_coeff);
    diffusion_body.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------  ------------------------------------------
    //	Particle and body creation of temperature observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation diffusion_body_inner(diffusion_body);
    ContactRelation diffusion_body_contact(diffusion_body, {&wall_boundary});
    ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation diffusion_body_complex(diffusion_body_inner, diffusion_body_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> diffusion_body_normal_direction(diffusion_body);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    InteractionSplit<TemperatureSplittingByPDEWithBoundary<Real>>
        temperature_splitting(diffusion_body_inner, diffusion_body_contact, "Phi");
    GetDiffusionTimeStepSize get_time_step_size(diffusion_body, *local_isotropic_diffusion);
    SimpleDynamics<DiffusionBodyInitialCondition> setup_diffusion_initial_condition(diffusion_body);
    SimpleDynamics<WallBoundaryInitialCondition> setup_boundary_condition(wall_boundary);

    ReduceDynamics<Average<QuantitySummation<Real, SPHBody>>>
        calculate_averaged_temperature(diffusion_body, "Phi");
    BodyRegionByParticle heat_flux_region(diffusion_body, makeShared<MultiPolygonShape>(createHeatFluxBoundary(), "HeatFluxRegion"));
    ReduceDynamics<Average<QuantitySummation<Real, BodyPartByParticle>>>
        calculate_boundary_averaged_temperature(heat_flux_region, "Phi");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    RestartIO restart_io(sph_system);
    ObservedQuantityRecording<Real> write_solid_temperature("Phi", temperature_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    setup_diffusion_initial_condition.exec();
    setup_boundary_condition.exec();
    diffusion_body_normal_direction.exec();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        diffusion_body.updateCellLinkedList();
        diffusion_body_complex.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = sph_system.RestartStep();
    Real T0 = 10;
    Real End_Time = T0;
    int restart_output_interval = 1000;
    Real dt = 0.0;
    Real current_averaged_temperature = 0.0;
    Real current_averaged_boundary_temperature = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TickCount::interval_t interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    std::string filefullpath_nonopt_temperature =
        sph_system.getIOEnvironment().output_folder_ + "/" + "nonopt_temperature.dat";
    std::ofstream out_file_nonopt_temperature(filefullpath_nonopt_temperature.c_str(), std::ios::app);
    std::string filefullpath_nonopt_boundary_temperature =
        sph_system.getIOEnvironment().output_folder_ + "/" + "nonopt_boundary_temperature.dat";
    std::ofstream out_file_nonopt_boundary_temperature(filefullpath_nonopt_boundary_temperature.c_str(), std::ios::app);

    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        dt = get_time_step_size.exec();
        if (ite % 500 == 0)
        {
            write_states.writeToFile(ite);
            write_solid_temperature.writeToFile(ite);

            current_averaged_temperature = calculate_averaged_temperature.exec();
            out_file_nonopt_temperature << std::fixed << std::setprecision(12) << ite << "   " << current_averaged_temperature << "\n";

            current_averaged_boundary_temperature = calculate_boundary_averaged_temperature.exec();
            out_file_nonopt_boundary_temperature << std::fixed << std::setprecision(12) << ite << "   " << current_averaged_boundary_temperature << "\n";

            std::cout << "N= " << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
            std::cout << "The averaged temperature is " << current_averaged_temperature << std::endl;
            std::cout << "The averaged boundary temperature is " << current_averaged_boundary_temperature << std::endl;
        }

        temperature_splitting.exec(dt);
        ite++;
        GlobalStaticVariables::physical_time_ += dt;

        if (ite % restart_output_interval == 0)
        {
            restart_io.writeToFile(ite);
        }
    }
    TickCount t4 = TickCount::now();
    TickCount::interval_t tt;
    tt = t4 - t1;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    std::cout << "Total physical time for computation: " << GlobalStaticVariables::physical_time_ << " seconds." << std::endl;

    EXPECT_NEAR(442.74, calculate_averaged_temperature.exec(), 0.01);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}