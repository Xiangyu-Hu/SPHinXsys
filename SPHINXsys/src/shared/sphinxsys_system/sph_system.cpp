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
}