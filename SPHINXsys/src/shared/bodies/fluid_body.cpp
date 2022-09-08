/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "fluid_body.h"

#include "base_material.h"
#include "fluid_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: RealBody(system, shape_ptr), iteration_count_(0) {}
	//=================================================================================================//
	void FluidBody::updateCellLinkedList()
	{
		// sorting is carried out once for 100 iterations
		if (iteration_count_ % 100 == 0)
			sortParticleWithCellLinkedList();
		iteration_count_++;
		cell_linked_list_->UpdateCellLists();
	}
	//=================================================================================================//
	EulerianFluidBody::EulerianFluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: RealBody(system, shape_ptr)
	{
		defineAdaptation<SPHAdaptation>(1.3);
	}
	//=================================================================================================//
}
