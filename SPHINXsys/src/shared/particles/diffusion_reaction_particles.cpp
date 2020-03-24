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
	DiffusionReactionData
		::DiffusionReactionData(size_t number_of_species) 
	{
		for (size_t i = 0; i < number_of_species; ++i) 
		{
			species_n_.push_back(0.0);
			dspecies_dt_.push_back(0.0);
			species_s_.push_back(0.0);
		}
	};
	//=================================================================================================//
	ElectroPhysiologyParticles
		::ElectroPhysiologyParticles(SPHBody* body, BaseMaterial* base_material) 
		: DiffusionReactionParticles<SolidParticles, Solid>(body, base_material)
	{
	}
	//=================================================================================================//
}