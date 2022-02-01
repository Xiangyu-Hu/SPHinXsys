/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "fluid_body.h"
#include "cell_linked_list.h"

namespace SPH {
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, const std::string &body_name,
					 SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: RealBody(system, body_name, sph_adaptation_ptr),
		iteration_count_(0) {}
	//=================================================================================================//
	void FluidBody::updateCellLinkedList()
	{
		//sorting is carried out once for 100 iterations
		if (iteration_count_ % 100 == 0) sortParticleWithCellLinkedList();
		iteration_count_++;
		cell_linked_list_->UpdateCellLists();
	}
	//=================================================================================================//
	EulerianFluidBody::EulerianFluidBody(SPHSystem &system, const std::string &body_name,
					 SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: RealBody(system, body_name, sph_adaptation_ptr) {}
	//=================================================================================================//
}
