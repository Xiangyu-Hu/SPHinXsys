#include "granular_body.h"
#include "base_material.h"
#include "fluid_particles.h"
#include "cell_linked_list.h"
namespace SPH
{
	GranularBody::GranularBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: FluidBody(system, shape_ptr) {}
}