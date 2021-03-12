/**
 * @file 	diffusion_reaction_particles.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "diffusion_reaction_particles.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	ElectroPhysiologyParticles
		::ElectroPhysiologyParticles(SPHBody* body, 
			DiffusionReactionMaterial<SolidParticles, Solid>* diffusion_reaction_material)
		: DiffusionReactionParticles<SolidParticles, Solid>(body, diffusion_reaction_material)
	{
	}
	//=================================================================================================//
}
