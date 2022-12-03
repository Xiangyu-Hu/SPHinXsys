#include "fluid_body.h"

#include "base_material.h"
#include "fluid_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: RealBody(system, shape_ptr) {}
	//=================================================================================================//
	EulerianFluidBody::EulerianFluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: FluidBody(system, shape_ptr)
	{
		defineAdaptation<SPHAdaptation>(1.3);
	}
	//=================================================================================================//
}
