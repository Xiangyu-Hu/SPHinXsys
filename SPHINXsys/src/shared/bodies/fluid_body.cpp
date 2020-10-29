/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_body.h"
#include "mesh_cell_linked_list.h"

namespace SPH {
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, string body_name,
		int refinement_level, ParticleGenerator* particle_generator)
		: RealBody(system, body_name, refinement_level, 1.3, particle_generator),
		iteration_count_(0)
	{

	}
	//=================================================================================================//
	void FluidBody::updateCellLinkedList()
	{
		//sorting is carried out once for 100 iterations
		if(iteration_count_% 100 == 0) mesh_cell_linked_list_->sortingParticleData();
		iteration_count_++;
		mesh_cell_linked_list_->UpdateCellLists();
	}
	//=================================================================================================//
}
