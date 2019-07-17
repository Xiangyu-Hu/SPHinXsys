#include "fluid_body.h"
#include "mesh_cell_linked_list.h"
#include "weakly_compressible_fluid.h"
#include "fluid_particles.h"
#include "sph_system.h"
#include "base_kernel.h"

namespace SPH {
	//===============================================================//
	FluidBody::FluidBody(SPHSystem &system, string body_name,
		Fluid &material, FluidParticles &fluid_particles,
			int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, material, fluid_particles,
			refinement_level, 1.3, op), fluid_material_(material),
		fluid_particles_(fluid_particles)
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
	FluidBodyPart
		::FluidBodyPart(FluidBody *fluid_body, string fluid_body_part_name)
		: EulerianBodyPart(fluid_body, fluid_body_part_name),
		fluid_body_(fluid_body), fluid_body_part_region_(fluid_body_part_name)
	{

	}
}