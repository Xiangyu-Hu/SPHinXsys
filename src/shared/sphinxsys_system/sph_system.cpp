#include "sph_system.h"

#include "all_body_relations.h"
#include "base_body.h"
#include "elastic_dynamics.h"

namespace SPH
{
//=================================================================================================//
SPHSystem::SPHSystem(BoundingBox system_domain_bounds, Real resolution_ref, size_t number_of_threads)
    : system_domain_bounds_(system_domain_bounds),
      resolution_ref_(resolution_ref),
      tbb_global_control_(tbb::global_control::max_allowed_parallelism, number_of_threads),
      io_environment_(nullptr), run_particle_relaxation_(false), reload_particles_(false),
      restart_step_(0), generate_regression_data_(false), state_recording_(true) {}
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
        for (size_t i = 0; i < body->body_relations_.size(); i++)
        {
            body->body_relations_[i]->updateConfiguration();
        }
    }
}
//=================================================================================================//
Real SPHSystem::getSmallestTimeStepAmongSolidBodies(Real CFL)
{
    Real dt = Infinity;
    for (size_t i = 0; i < solid_bodies_.size(); i++)
    {
        ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(*solid_bodies_[i], CFL);
        Real dt_temp = computing_time_step_size.exec();
        if (dt_temp < dt)
            dt = dt_temp;
    }
    return dt;
}
//=================================================================================================//
#ifdef BOOST_AVAILABLE
void SPHSystem::handleCommandlineOptions(int ac, char *av[])
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
}
#endif
//=================================================================================================//
} // namespace SPH
