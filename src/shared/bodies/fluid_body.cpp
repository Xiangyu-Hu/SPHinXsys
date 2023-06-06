#include "fluid_body.h"

#include "base_material.h"
 
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	FluidBody::FluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: RealBody(system, shape_ptr) {}
	//=================================================================================================//
}
