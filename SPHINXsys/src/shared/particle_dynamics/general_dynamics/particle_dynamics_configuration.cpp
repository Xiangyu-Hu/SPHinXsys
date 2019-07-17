#include "particle_dynamics_configuration.h"
#include "sph_system.h"
#include "base_body.h"

namespace SPH {

	ParticleDynamicsConfiguration::ParticleDynamicsConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
	//===============================================================//
	void ParticleDynamicsConfiguration::exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
		body_->UpdateContactConfiguration();
	}
	//===============================================================//
	void ParticleDynamicsConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
		body_->UpdateContactConfiguration();
	}
	//===============================================================//
	ParticleDynamicsContactConfiguration
		::ParticleDynamicsContactConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
	//===============================================================//
	void ParticleDynamicsContactConfiguration::exec(Real dt)
	{
		body_->UpdateContactConfiguration();
	}
	//===============================================================//
	void ParticleDynamicsContactConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateContactConfiguration();
	}
	//===============================================================//
	ParticleDynamicsInteractionConfiguration
		::ParticleDynamicsInteractionConfiguration(SPHBody *body, 
			SPHBodyVector interacting_bodies)
		: ParticleDynamics<void, SPHBody>(body), interacting_bodies_(interacting_bodies)
	{
		//building interacting topology
		interacting_bodies_ = interacting_bodies;
	}
	//===============================================================//
	void ParticleDynamicsInteractionConfiguration::exec(Real dt)
	{
		body_->UpdateInteractionConfiguration(interacting_bodies_);
	}
	//===============================================================//
	void ParticleDynamicsInteractionConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInteractionConfiguration(interacting_bodies_);
	}
	//===============================================================//
	ParticleDynamicsInnerConfiguration
		::ParticleDynamicsInnerConfiguration(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
	//===============================================================//
	void ParticleDynamicsInnerConfiguration::exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
	}
	//===============================================================//
	void ParticleDynamicsInnerConfiguration::parallel_exec(Real dt)
	{
		body_->UpdateInnerConfiguration();
	}
	//===============================================================//
	ParticleDynamicsCellLinkedList
		::ParticleDynamicsCellLinkedList(SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body)
	{

	}
	//===============================================================//
	void ParticleDynamicsCellLinkedList::exec(Real dt)
	{
		body_->UpdateCellLinkedList();
	}
	//===============================================================//
	void ParticleDynamicsCellLinkedList::parallel_exec(Real dt)
	{
		body_->UpdateCellLinkedList();
	}
	//===============================================================//
	ParticleDynamicsByCellParticleLists
		::ParticleDynamicsByCellParticleLists(SPHSystem *system, SPHBody *body)
		: ParticleDynamics<void, SPHBody>(body), system_(system)
	{

	}
	//===============================================================//
	void ParticleDynamicsByCellParticleLists::exec(Real dt)
	{
		cout << "\n ParticleDynamicsByCellParticleLists not done yet. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
	void ParticleDynamicsByCellParticleLists::parallel_exec(Real dt)
	{
		cout << "\n ParticleDynamicsByCellParticleLists not done yet. Exit the program! \n";
		exit(0);
	}
}