/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	diffusion_reaction_particles.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef DIFFUSION_REACTION_PARTICLES_H
#define DIFFUSION_REACTION_PARTICLES_H

#include "base_particles.h"
#include "base_body.h"
#include "base_material.h"
#include "diffusion_reaction.h"

namespace SPH
{

	/**
	 * @class DiffusionReactionParticles
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template <class BaseParticlesType, class BaseMaterialType = BaseMaterial, int NUM_SPECIES = 1>
	class DiffusionReactionParticles : public BaseParticlesType
	{
	protected:
		size_t number_of_species_;			 /**< Total number of diffusion and reaction species . */
		size_t number_of_diffusion_species_; /**< Total number of diffusion species . */
		std::map<std::string, size_t> species_indexes_map_;

	public:
		StdVec<StdLargeVec<Real>> species_n_;	 /**< array of diffusion/reaction scalars */
		StdVec<StdLargeVec<Real>> diffusion_dt_; /**< array of the time derivative of diffusion species */
		DiffusionReaction<BaseMaterialType, NUM_SPECIES> &diffusion_reaction_material_;

		DiffusionReactionParticles(SPHBody &sph_body,
								   DiffusionReaction<BaseMaterialType, NUM_SPECIES> *diffusion_reaction_material)
			: BaseParticlesType(sph_body, diffusion_reaction_material),
			  number_of_species_(diffusion_reaction_material->NumberOfSpecies()),
			  number_of_diffusion_species_(diffusion_reaction_material->NumberOfSpeciesDiffusion()),
			  species_indexes_map_(diffusion_reaction_material->SpeciesIndexMap()),
			  diffusion_reaction_material_(*diffusion_reaction_material)
		{
			species_n_.resize(number_of_species_);
			diffusion_dt_.resize(number_of_diffusion_species_);
		};
		virtual ~DiffusionReactionParticles(){};

		std::map<std::string, size_t> SpeciesIndexMap() { return species_indexes_map_; };

		virtual void initializeOtherVariables() override
		{
			BaseParticlesType::initializeOtherVariables();

			std::map<std::string, size_t>::iterator itr;
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr)
			{
				// Register a specie.
				this->registerVariable(species_n_[itr->second], itr->first);
				/** the scalars will be sorted if particle sorting is called, Note that we call a template function from a template class. */
				this->template registerSortableVariable<Real>(itr->first);
				/** add species to basic output particle data. */
				this->template addVariableToWrite<Real>(itr->first);
			}

			constexpr int type_index = DataTypeIndex<Real>::value;
			for (size_t m = 0; m < number_of_diffusion_species_; ++m)
			{
				// register reactive change rate terms without giving variable name
				std::get<type_index>(this->all_particle_data_).push_back(&diffusion_dt_[m]);
				diffusion_dt_[m].resize(this->real_particles_bound_, Real(0));
			}
		};

		virtual DiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES> *ThisObjectPtr() override { return this; };
	};
}
#endif // DIFFUSION_REACTION_PARTICLES_H