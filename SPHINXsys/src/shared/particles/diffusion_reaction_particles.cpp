/**
 * @file 	diffusion_reaction_particles.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "diffusion_reaction_particles.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	ElectroPhysiologyParticles::
		ElectroPhysiologyParticles(
			SPHBody &sph_body,
			SharedPtr<DiffusionReaction<SolidParticles, Solid>> diffusion_reaction_material_ptr,
			SharedPtr<ParticleGenerator> particle_generator_ptr)
		: DiffusionReactionParticles<SolidParticles, Solid>(sph_body,diffusion_reaction_material_ptr, particle_generator_ptr)
	{
	}
	//=================================================================================================//
	ElectroPhysiologyReducedParticles::
		ElectroPhysiologyReducedParticles(
			SPHBody &sph_body,
			SharedPtr<DiffusionReaction<SolidParticles, Solid>> diffusion_reaction_material_ptr,
			SharedPtr<ParticleGenerator> particle_generator_ptr)
		: DiffusionReactionParticles<SolidParticles, Solid>(sph_body,diffusion_reaction_material_ptr, particle_generator_ptr)
	{
	}
	//=================================================================================================//
}