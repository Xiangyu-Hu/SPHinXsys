/**
 * @file 	diffusion_reaction_particles.cpp
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
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