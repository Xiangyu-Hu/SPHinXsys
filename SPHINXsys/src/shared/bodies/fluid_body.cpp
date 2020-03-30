/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_body.h"
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
		base_mesh_cell_linked_list_->UpdateInnerConfiguration(current_configuration_);
	}
	//===============================================================//
	void FluidBody::BuildContactConfiguration()
	{
		base_mesh_cell_linked_list_->UpdateContactConfiguration();
	}
	//===============================================================//
}