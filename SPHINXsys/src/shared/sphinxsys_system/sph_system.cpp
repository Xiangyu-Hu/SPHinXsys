/**
 * @file sph_system.cpp
 * @brief 	Definatioin of all the functions decleared in spy_system.h
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */

#include "sph_system.h"
#include "base_body.h"
#include "particle_generator_lattice.h"

namespace SPH
{
	//===============================================================//
	SPHSystem::SPHSystem(Vecd lower_bound, Vecd upper_bound,
		Real particle_spacing_ref, int number_of_threads)
		: lower_bound_(lower_bound), upper_bound_(upper_bound),
		particle_spacing_ref_(particle_spacing_ref), tbb_init_(number_of_threads),
		restart_step_(0), run_particle_relaxation_(false),
		reload_particles_(false)
	{
	}
	//===============================================================//
	SPHSystem::~SPHSystem()
	{

	}
	//===============================================================//
	void SPHSystem::AddBody(SPHBody* body)
	{
		bodies_.push_back(body);
	}
	//===============================================================//
	void SPHSystem::AddRealBody(SPHBody* body)
	{
		real_bodies_.push_back(body);
	}
	//===============================================================//
	void SPHSystem::AddFictitiousBody(SPHBody* body)
	{
		fictitious_bodies_.push_back(body);
	}
	//===============================================================//
	void SPHSystem::SetBodyTopology(SPHBodyTopology* body_topology)
	{
		body_topology_ = body_topology;
		for (size_t i = 0; i < body_topology_->size(); i++)
		{
			for (auto& body : bodies_) {
				if (body == body_topology_->at(i).first) {
					body->SetContactMap(body_topology_->at(i));
					body->AllocateMeoemryCellLinkedList();
					body->AllocateMemoriesForConfiguration();
				}
			}
		}
	}
	//===============================================================//
	void SPHSystem::InitializeSystemCellLinkedLists()
	{
		for (auto &body : bodies_)
		{
			body->UpdateCellLinkedList();
		}
	}
	//===============================================================//
	void SPHSystem::InitializeSystemConfigurations()
	{
		for (size_t i = 0; i < body_topology_->size(); i++)
		{
			SPHBody *body = body_topology_->at(i).first;
			body->BuildInnerConfiguration();
			body->BuildContactConfiguration();
		}
	}
	//===============================================================//
	void SPHSystem::SetupSPHSimulation()
	{
		InitializeSystemCellLinkedLists();
		InitializeSystemConfigurations();
	}
	//===============================================================//
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
				cout << "Particle relaxation was set to default(false).\n";
			}
			if (vm.count("i")) {
				reload_particles_ = vm["i"].as<bool>();
				cout << "Particle reload from input file was set to "
					<< vm["i"].as<bool>() << ".\n";
			}
			else {
				cout << "Particle reload from input file was set to default(false).\n";
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
	//===============================================================//
}