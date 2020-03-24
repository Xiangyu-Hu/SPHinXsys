/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_body.h"
#include "sph_system.h"
#include "mesh_cell_linked_list.h"

namespace SPH {
	//===============================================================//
	FluidBody::FluidBody(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, refinement_level, 1.3, op)
	{

	}
	//===============================================================//
	void FluidBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this, current_inner_configuration_);
	}
	//===============================================================//
	void FluidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
}