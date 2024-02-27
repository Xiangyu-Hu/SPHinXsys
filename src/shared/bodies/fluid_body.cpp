#include "fluid_body.h"

#include "base_material.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
FluidBody::FluidBody(SPHSystem &system, SharedPtr<Shape> initial_shape_ptr)
    : RealBody(system, initial_shape_ptr) {}
//=================================================================================================//
} // namespace SPH
