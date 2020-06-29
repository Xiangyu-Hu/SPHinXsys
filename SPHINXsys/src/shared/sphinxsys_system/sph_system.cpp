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
	//=================================================================================================//
	SPHBodySystem::SPHBodySystem(SPHSystem* sph_system)
		: sph_system_(sph_system)
	{
	}
	//=================================================================================================//
	bool SPHBodySystem::addABody(SPHBody* sph_body) {
		bool new_sph_body = true;
		for (size_t i = 0; i != sph_bodies_.size(); ++i)
			if (sph_body == sph_bodies_[i]) new_sph_body = false;
		if (new_sph_body) {
			sph_bodies_.push_back(sph_body);
		}
		return new_sph_body;
	}
	//=================================================================================================//
	void SPHBodySystem::addBodies(SPHBodyVector sph_bodies) {
		for (size_t i = 0; i != sph_bodies.size(); ++i)
			addABody(sph_bodies[i]);
	}
	//=================================================================================================//
	void SPHBodySystem::addSPHBodyContactRealtion(SPHBodyContactRealtion* sph_body_contact_realtion)
	{
		/** add a new sph bodies to the system */
		addABody(sph_body_contact_realtion->body_);
		addBodies(sph_body_contact_realtion->contact_bodies_);
		sph_body_contact_realtions_.push_back(sph_body_contact_realtion);
	}
	//=================================================================================================//
	bool SPHBodyCollisionSystem::addABody(SPHBody* sph_body)
	{
		bool new_sph_body = SPHBodySystem::addABody(sph_body);
		if (new_sph_body) {
			sph_bodies_.push_back(sph_body);
			sph_body_lower_bounds_.push_back(BodyLowerBound(sph_body));
			sph_body_upper_bounds_.push_back(BodyUpperBound(sph_body));
		}
		return new_sph_body;
	}
	//=================================================================================================//
	void SPHBodyCollisionSystem::updateBodyBound()
	{
		for (size_t i = 0; i != sph_bodies_.size(); ++i) {
			if (!sph_bodies_[i]->is_static_) {
				sph_bodies_[i]->setBodyLowerBound(sph_body_lower_bounds_[i].parallel_exec());
				sph_bodies_[i]->setBodyUpperBound(sph_body_upper_bounds_[i].parallel_exec());
			}
		}
	}
	void SPHBodySystem::updateCellLinkedList()
	{
		for (size_t i = 0; i != sph_bodies_.size(); ++i) {
			if (!sph_bodies_[i]->is_static_) {
				sph_bodies_[i]->UpdateCellLinkedList();
			}
		}
	}
	//=================================================================================================//
	SPHSystem::SPHSystem(Vecd lower_bound, Vecd upper_bound,
		Real particle_spacing_ref, int number_of_threads)
		: sph_body_systems_(NULL), lower_bound_(lower_bound), upper_bound_(upper_bound),
		tbb_init_(number_of_threads), particle_spacing_ref_(particle_spacing_ref),
		restart_step_(0), run_particle_relaxation_(false),
		reload_particles_(false)
	{
	}
	//=================================================================================================//
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
					body->AllocateMemoriesForContactConfiguration();
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
				cout << "Particle relaxation was set to default (" << run_particle_relaxation_ <<").\n";
			}
			if (vm.count("i")) {
				reload_particles_ = vm["i"].as<bool>();
				cout << "Particle reload from input file was set to "
					<< vm["i"].as<bool>() << ".\n";
			}
			else {
				cout << "Particle reload from input file was set to default (" << reload_particles_ << ").\n";
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