/**
 * @file 	observer_body.cpp
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
		int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{

	}
	//===========================================================//
	ObserverLagrangianBody::ObserverLagrangianBody(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: ObserverBody(system, body_name, refinement_level, op)
	{

	}
	//===========================================================//
	void ObserverLagrangianBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceContactConfiguration(*this);
	}
	//===========================================================//
	ObserverEulerianBody::ObserverEulerianBody(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: ObserverBody(system, body_name, refinement_level, op)
	{

	}
	//===========================================================//
	void ObserverEulerianBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===========================================================//
}