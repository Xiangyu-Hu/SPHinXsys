/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	diffusion_reaction_particles.h
* @brief 	This is the derived class of diffusion reaction particles.
* @author	Xiangyu Huand Chi Zhang
*/


#ifndef DIFFUSION_REACTION_PARTICLES_H
#define DIFFUSION_REACTION_PARTICLES_H



#include "base_particles.h"
#include "base_body.h"
#include "base_material.h"
#include "diffusion_reaction.h"
#include "xml_engine.h"
#include "solid_particles.h"

using namespace std;
namespace SPH {

	/**
	 * @class DiffusionReactionParticles
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template<class BaseParticlesType = BaseParticles, class BaseMaterialType = BaseMaterial>
	class DiffusionReactionParticles : public BaseParticlesType
	{
	protected:
		size_t number_of_species_;				/**< Total number of diffusion and reaction species . */
		size_t number_of_diffusion_species_;	/**< Total number of diffusion species . */
		map<string, size_t> species_indexes_map_;
	public:
		StdVec<StdLargeVec<Real>> species_n_;	/**< array of diffusion/reaction scalars */
		StdVec<StdLargeVec<Real>> diffusion_dt_;/**< array of the time derivative of diffusion species */

		DiffusionReactionParticles(SPHBody* body, 
			DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>* diffusion_reaction_material)
			: BaseParticlesType(body, diffusion_reaction_material)
		{
			diffusion_reaction_material->assignDiffusionReactionParticles(this);
			species_indexes_map_ = diffusion_reaction_material->SpeciesIndexMap();
			number_of_species_ = diffusion_reaction_material->NumberOfSpecies();
			species_n_.resize(number_of_species_);

			map<string, size_t>::iterator itr;
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) 
			{
				//Register a species. Note that we call a template function from a template class
				this->template registerAVariable<indexScalar, Real>(species_n_[itr->second], itr->first, true);
				//the scalars will be sorted if particle sorting is called
				this->sortable_scalars_.push_back(&species_n_[itr->second]);
			}

			number_of_diffusion_species_ = diffusion_reaction_material->NumberOfSpeciesDiffusion();
			diffusion_dt_.resize(number_of_diffusion_species_);
			for (size_t m = 0; m < number_of_diffusion_species_; ++m)
			{
				//----------------------------------------------------------------------
				//	register reactive change rate terms without giving variable name
				//----------------------------------------------------------------------
				std::get<indexScalar>(this->all_particle_data_).push_back(&diffusion_dt_[m]);
				diffusion_dt_[m].resize(this->real_particles_bound_, Real(0));
			}
		};
		virtual ~DiffusionReactionParticles() {};

		map<string, size_t> SpeciesIndexMap() { return  species_indexes_map_; };

		/** Write particle data in PLT format for Tecplot. */
		virtual void writeParticlesToPltFile(ofstream& output_file) override
		{
			size_t total_real_particles = this->total_real_particles_;
			map<string, size_t>::iterator itr;
			
			if(Vecd(0).size()==2)
			{
				output_file << " VARIABLES = \" x \", \"y\", \"ID \" " << " ";
			}
			if(Vecd(0).size()==3)
			{
				output_file << " VARIABLES = \" x \", \"y\" , \"z\", \" ID \" " << " ";
			}

			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) 
			{
				output_file << ","<< "\"" << itr->first << " \" " << "";
			}
			output_file << "\n";

			for (size_t i = 0; i != total_real_particles; ++i)
			{
				for(int j = 0; j != Vecd(0).size(); ++j)
				{
					output_file << this->pos_n_[i][j] << "  ";
				}

				output_file << i << " ";

				for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr) 
				{
					size_t k = itr->second;
					output_file << species_n_[i][k] << " ";
				}

				output_file << "\n";
			}
		};
		virtual void writeParticlesToXmlForRestart(std::string& filefullpath) override {};
		virtual void readParticleFromXmlForRestart(std::string& filefullpath) override {};
		virtual DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>*  
			pointToThisObject() override { return this; };
	};

	/**
	 * @class ElectroPhysiologyParticles
	 * @brief A group of particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyParticles 
		: public DiffusionReactionParticles<SolidParticles, Solid>
	{
	public:
		ElectroPhysiologyParticles(SPHBody* body, 
			DiffusionReactionMaterial<SolidParticles, Solid>* diffusion_reaction_material);
		virtual ~ElectroPhysiologyParticles() {};
		virtual ElectroPhysiologyParticles* pointToThisObject() override { return this; };

	};
}
#endif //DIFFUSION_REACTION_PARTICLES_H