#include "sph_system.hpp"

#include "base_body_relation.h"
#include "geometric_shape.h"
#include "io_environment.h"
#include "predefined_bodies.h"

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>
#define TBB_PARALLEL true
#ifdef BOOST_AVAILABLE
#include "boost/program_options.hpp"
namespace po = boost::program_options;
#endif

namespace SPH
{
//=================================================================================================//
namespace
{
SharedPtr<tbb::global_control> &getTbbGlobalControlHolder()
{
    // Intentionally keep this alive until process termination to avoid
    // static destruction-order issues inside oneTBB global control teardown.
    static SharedPtr<tbb::global_control> *holder = new SharedPtr<tbb::global_control>();
    return *holder;
}
} // namespace
//=================================================================================================//
SPHSystem::SPHSystem(BoundingBoxd system_domain_bounds, Real global_resolution, size_t number_of_threads)
    : SPHSystem(true, system_domain_bounds, global_resolution, number_of_threads)
{
    writeSystemDomainShapeToVtp();
}
//=================================================================================================//
SPHSystem::SPHSystem(bool is_physical, BoundingBoxd system_domain_bounds,
                     Real global_resolution, size_t number_of_threads)
    : system_name_("SPHSystem"),
      system_bounds_(system_domain_bounds.expand(global_resolution * 4)),
      global_resolution_(global_resolution),
      is_physical_(is_physical), run_particle_relaxation_(false), reload_particles_(false),
      restart_step_(0), generate_regression_data_(false), state_recording_(true)
{
    IO::initEnvironment();
    IO::initLogger();
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level_));
    sv_physical_time_ = registerSystemVariable<Real>("PhysicalTime", 0.0);
    IO::getLogger()->info("The reference resolution of the SPHSystem is {}.", global_resolution_);
    getTbbGlobalControlHolder() = std::make_shared<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism, number_of_threads);
}
//=================================================================================================//
SPHSystem::~SPHSystem() = default;
//=================================================================================================//
void SPHSystem::writeSystemDomainShapeToVtp()
{
    GeometricShapeBox domain_shape(system_bounds_, system_name_ + "Domain");
    domain_shape.writeGeometricShapeBoxToVtp();
}
//=================================================================================================//
void SPHSystem::setLogLevel(size_t log_level)
{
    if (log_level < 0 || log_level > 6)
    {
        std::cerr << "Log level must be between 0 and 6.\n";
        exit(1);
    }

    log_level_ = log_level;
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level_));
}
//=================================================================================================//
void SPHSystem::addRealBody(RealBody *real_body)
{
    real_bodies_.push_back(real_body);
}
//=================================================================================================//

void SPHSystem::initializeSystemCellLinkedLists()
{
    for (auto &body : real_bodies_)
    {
        DynamicCast<RealBody>(this, body)->updateCellLinkedList();
    }
}
//=================================================================================================//
void SPHSystem::initializeSystemConfigurations()
{
    for (auto &body : sph_bodies_)
    {
        StdVec<SPHRelation *> &body_relations = body->getBodyRelations();
        for (size_t i = 0; i < body_relations.size(); i++)
        {
            body_relations[i]->updateConfiguration();
        }
    }
}
//=================================================================================================//
#ifdef BOOST_AVAILABLE
SPHSystem *SPHSystem::handleCommandlineOptions(int ac, char *av[])
{
    try
    {

        po::options_description desc("Allowed options");
        desc.add_options()("help", "produce help message");
        desc.add_options()("relax", po::value<bool>(), "Particle relaxation.");
        desc.add_options()("reload", po::value<bool>(), "Particle reload from input file.");
        desc.add_options()("regression", po::value<bool>(), "Regression test.");
        desc.add_options()("state_recording", po::value<bool>(), "State recording in output folder.");
        desc.add_options()("restart_step", po::value<int>(), "Run form a restart file.");
        desc.add_options()("log_level", po::value<int>(), "Output log level (0-6). "
                                                          "0: trace, 1: debug, 2: info, 3: warning, 4: error, 5: critical, 6: off");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            exit(0);
        }

        if (vm.count("relax"))
        {
            run_particle_relaxation_ = vm["relax"].as<bool>();
            std::cout << "Particle relaxation was set to "
                      << vm["relax"].as<bool>() << ".\n";
        }
        else
        {
            std::cout << "Particle relaxation was set to default ("
                      << run_particle_relaxation_ << ").\n";
        }

        if (run_particle_relaxation_)
        {
            IO::getEnvironment().reinitializeReloadFolder();
        }

        if (vm.count("reload"))
        {
            reload_particles_ = vm["reload"].as<bool>();
            std::cout << "Particle reload from input file was set to "
                      << vm["reload"].as<bool>() << ".\n";
        }
        else
        {
            std::cout << "Particle reload from input file was set to default ("
                      << reload_particles_ << ").\n";
        }

        if (vm.count("regression"))
        {
            generate_regression_data_ = vm["regression"].as<bool>();
            std::cout << "Generate regression test data set was set to "
                      << vm["regression"].as<bool>() << ".\n";
        }
        else
        {
            std::cout << "Generate regression test data was set to default ("
                      << generate_regression_data_ << ").\n";
        }

        if (vm.count("state_recording"))
        {
            state_recording_ = vm["state_recording"].as<bool>();
            std::cout << "State recording was set to "
                      << vm["state_recording"].as<bool>() << ".\n";
        }
        else
        {
            std::cout << "State recording was set to default ("
                      << state_recording_ << ").\n";
        }

        if (vm.count("restart_step"))
        {
            restart_step_ = vm["restart_step"].as<int>();
            std::cout << "Restart step was set to "
                      << vm["restart_step"].as<int>() << ".\n";
        }
        else
        {
            std::cout << "Restart inactivated, i.e. restart_step ("
                      << restart_step_ << ").\n";
        }

        if (vm.count("log_level"))
        {
            log_level_ = vm["log_level"].as<int>();
            if (log_level_ < 0 || log_level_ > 6)
            {
                std::cerr << "Log level must be between 0 and 6.\n";
                exit(1);
            }
            std::cout << "Log level was set to " << log_level_ << ".\n";
            spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level_));
        }
        else
        {
            std::cout << "Log level was set to default (info).\n";
        }
    }
    catch (std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    catch (...)
    {
        std::cerr << "Exception of unknown type!\n";
    }

    return this;
}
#endif
//=================================================================================================//
RelaxationSystem::RelaxationSystem(
    BoundingBoxd system_domain_bounds, Real global_resolution, size_t number_of_threads)
    : SPHSystem(false, system_domain_bounds, global_resolution, number_of_threads)
{
    system_name_ = "RelaxationSystem";
    writeSystemDomainShapeToVtp();
}
//=================================================================================================//
} // namespace SPH
