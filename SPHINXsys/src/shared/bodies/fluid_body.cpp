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
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this);
	}
	//===============================================================//
	void FluidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
	FluidBodyPart::FluidBodyPart(FluidBody *fluid_body, string fluid_body_part_name)
		: BodyPartByCell(fluid_body, fluid_body_part_name), fluid_body_(fluid_body)
	{

	}
	//===============================================================//
}