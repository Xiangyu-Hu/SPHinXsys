/**
 * @file 	observer_body.cpp
 * @brief 	This is the class for bodies recording the state of the flow or solid
 * 			in given locations. This body has no inner configuration so that no
 * 			cell linked list is requires. However, it has contact configuration to
 * 			the body it is obesering.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "observer_body.h"
#include "observer_particles.h"
#include "mesh_cell_linked_list.h"
#include "sph_system.h"

namespace SPH {
	//===========================================================//
	ObserverBody::ObserverBody(SPHSystem &system, string body_name,
		ObserverParticles &observer_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, observer_particles, refinement_level, op),
		observer_particles_(observer_particles)
	{

	}
	//===========================================================//
	ObserverLagrangianBody::ObserverLagrangianBody(SPHSystem &system, string body_name,
		ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op)
		: ObserverBody(system, body_name, observer_particles, refinement_level, op)
	{

	}
	//===========================================================//
	void ObserverLagrangianBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceContactConfiguration(*this);
	}
	//===========================================================//
	ObserverEulerianBody::ObserverEulerianBody(SPHSystem &system, string body_name,
		ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op)
		: ObserverBody(system, body_name, observer_particles, refinement_level, op)
	{

	}
	//===========================================================//
	void ObserverEulerianBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===========================================================//
}