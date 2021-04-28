/**
 * @file sph_system.cpp
 * @brief 	Definition of all system level functions
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */

#include "sph_system.h"

#include "base_body.h"
#include "body_relation.h"

namespace SPH
{
	//=================================================================================================//
	SPHSystem::SPHSystem(BoundingBox system_domain_bounds,
		Real resolution_ref, size_t number_of_threads) :
		system_domain_bounds_(system_domain_bounds),
		resolution_ref_(resolution_ref), 
		tbb_global_control_(tbb::global_control::max_allowed_parallelism, number_of_threads),
		restart_step_(0), run_particle_relaxation_(false),
		reload_particles_(false)
	{
		output_folder_ = "./output";
		if (!fs::exists(output_folder_))
		{
			fs::create_directory(output_folder_);
		}

		restart_folder_ = "./restart";
		if (!fs::exists(restart_folder_))
		{
			fs::create_directory(restart_folder_);
		}

		reload_folder_ = "./reload";
	}
	//=================================================================================================//
	void SPHSystem::addABody(SPHBody* sph_body)
	{
		bodies_.push_back(sph_body);
	}
	//=================================================================================================//
	void SPHSystem::addARealBody(RealBody* real_body)
	{
		real_bodies_.push_back(real_body);
	}
	//=================================================================================================//
	void SPHSystem::addAFictitiousBody(FictitiousBody* fictitious_body)
	{
		fictitious_bodies_.push_back(fictitious_body);
	}
	//=================================================================================================//
	void SPHSystem::initializeSystemCellLinkedLists()
	{
		for (auto &body : real_bodies_)
		{
			dynamic_cast<RealBody*>(body)->updateCellLinkedList();
		}
	}
	//=================================================================================================//
	void SPHSystem::initializeSystemConfigurations()
	{
		for (auto& body : bodies_)
		{
			for (size_t i = 0; i < body->body_relations_.size(); i++)
			{
				body->body_relations_[i]->updateConfiguration();
			}

		}
	}
	//=================================================================================================//
	void SPHSystem::handleCommandlineOptions(int ac, char* av[])
	{
		try {

			po::options_description desc("Allowed options");
			desc.add_options()
				("help", "produce help message")
				("r", po::value<bool>(), "Particle relaxation.")
				("i", po::value<bool>(), "Particle reload from input file.")
				;

			po::variables_map vm;
			po::store(po::parse_command_line(ac, av, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << desc << "\n";
				exit(0);
			}

			if (vm.count("r")) {
				run_particle_relaxation_ = vm["r"].as<bool>();
				cout << "Particle relaxation was set to "
					 << vm["r"].as<bool>() << ".\n";
			}
			else {
				cout << "Particle relaxation was set to default (" 
				     << run_particle_relaxation_ <<").\n";
			}
			if (vm.count("i")) {
				reload_particles_ = vm["i"].as<bool>();
				cout << "Particle reload from input file was set to "
					 << vm["i"].as<bool>() << ".\n";
			}
			else {
				cout << "Particle reload from input file was set to default (" 
				     << reload_particles_ << ").\n";
			}
		}
		catch (std::exception & e) {
			cerr << "error: " << e.what() << "\n";
			exit(1);
		}
		catch (...) {
			cerr << "Exception of unknown type!\n";
		}

	}
	//=================================================================================================//
}
