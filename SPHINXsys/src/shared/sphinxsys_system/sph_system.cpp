/**
 * @file sph_system.cpp
 * @brief 	Definatioin of all the functions decleared in spy_system.h
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */
#include "sph_system.h"
#include "all_kernels.h"
#include "base_body.h"
#include "particle_generator_lattice.h"

namespace SPH
{
	//===============================================================//
	SPHSystem::SPHSystem(Vecd lower_bound, Vecd upper_bound,
		Real particle_spacing_ref, Real smoothinglength_ratio)
		: lower_bound_(lower_bound), upper_bound_(upper_bound),
		particle_spacing_ref_(particle_spacing_ref),
		smoothinglength_ratio_(smoothinglength_ratio)
	{
	}
	//===============================================================//
	SPHSystem::~SPHSystem()
	{

	}
	//===============================================================//
	Kernel* SPHSystem::GenerateAKernel(Real particle_spacing)
	{	
		return new KernelWendlandC2(particle_spacing*smoothinglength_ratio_);
		//return new KernelHyperbolic(particle_spacing*smoothinglength_ratio_);
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
	}
	//===============================================================//
	void SPHSystem::CreateParticelsForAllBodies()
	{
		for (auto &body : bodies_)
		{
			body->CreateParticelsInSpecificManner();
		}
	}
	//===============================================================//
	void SPHSystem::InitializeAllRealBodies()
	{
		for (auto &body : real_bodies_)
		{
			dynamic_cast<RealBody*>
				(body->PointToThisObject())->InitializeLocalMaterialProperties();
			dynamic_cast<RealBody*>
				(body->PointToThisObject())->InitialCondition();
		}
	}
	//===============================================================//
	void SPHSystem::InitializeSystemCellLinkedLists()
	{
		for (auto &body : bodies_)
		{
			body->AllocateMeoemryCellLinkedList();
			body->UpdateCellLinkedList();
		}
	}
	//===============================================================//
	void SPHSystem::InitializeSystemConfigurations()
	{
		for (size_t i = 0; i < body_topology_->size(); i++)
		{
			SPHBody *body = body_topology_->at(i).first;
			body->SetContactMap(body_topology_->at(i));
			body->AllocateMemoriesForConfiguration();
			body->BuildInnerConfiguration();
			body->BuildContactConfiguration();
		}
	}
	//===============================================================//
	void SPHSystem::SetupSPHSimulation()
	{
		CreateParticelsForAllBodies();
		InitializeSystemCellLinkedLists();
		InitializeSystemConfigurations();
		InitializeAllRealBodies();
	}
	//===============================================================//
	void SPHSystem::ReInitializeAllRealBodiesFromRestart()
	{
		for (auto &body : real_bodies_)
		{
			dynamic_cast<RealBody*>
				(body->PointToThisObject())->InitialConditionFromRestartFile();
		}
	}
	//===============================================================//
	void SPHSystem::ResetSPHSimulationFromRestart()
	{
		ReInitializeAllRealBodiesFromRestart();
	}
	//===============================================================//
}