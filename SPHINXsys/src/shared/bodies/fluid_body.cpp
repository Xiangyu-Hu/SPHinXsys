/**
 * @file 	fluid_body.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "fluid_body.h"
#include "cell_linked_list.h"

namespace SPH {
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, std::string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: RealBody(system, body_name, particle_adaptation, particle_generator),
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
	EulerianFluidBody::EulerianFluidBody(SPHSystem &system, std::string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: RealBody(system, body_name, particle_adaptation, particle_generator) {}
	//=================================================================================================//
}
