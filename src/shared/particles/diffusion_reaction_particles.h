/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	diffusion_reaction_particles.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_REACTION_PARTICLES_H
#define DIFFUSION_REACTION_PARTICLES_H

#include "base_body.h"
#include "base_material.h"
#include "base_particles.h"
#include "diffusion_reaction.h"

namespace SPH
{

/**
 * @class DiffusionReactionParticles
 * @brief A group of particles with diffusion or/and reactions particle data.
 */
template <class BaseParticlesType, class DiffusionReactionMaterialType>
class DiffusionReactionParticles : public BaseParticlesType
{
  protected:
    StdVec<std::string> all_species_names_;
    std::map<std::string, size_t> species_indexes_map_;
    StdVec<StdLargeVec<Real> *> diffusion_species_;
    StdVec<StdLargeVec<Real> *> gradient_species_;
    StdVec<StdLargeVec<Real> *> reactive_species_;

  public:
    StdVec<StdLargeVec<Real>> all_species_; /**< array of diffusion/reaction scalars */
    DiffusionReactionMaterialType &diffusion_reaction_material_;
    static constexpr int NumReactiveSpecies = DiffusionReactionMaterialType::NumReactiveSpecies;
    typedef DiffusionReactionMaterialType DiffusionReactionMaterial;

    DiffusionReactionParticles(SPHBody &sph_body, DiffusionReactionMaterialType *diffusion_reaction_material)
        : BaseParticlesType(sph_body, diffusion_reaction_material),
          all_species_names_(diffusion_reaction_material->AllSpeciesNames()),
          species_indexes_map_(diffusion_reaction_material->AllSpeciesIndexMap()),
          diffusion_reaction_material_(*diffusion_reaction_material)
    {
        all_species_.resize(all_species_names_.size());
        const IndexVector &diffusion_species_indexes = diffusion_reaction_material_.DiffusionSpeciesIndexes();
        const IndexVector &gradient_species_indexes = diffusion_reaction_material_.GradientSpeciesIndexes();
        for (size_t i = 0; i != diffusion_species_indexes.size(); ++i)
        {
            diffusion_species_.push_back(&all_species_[diffusion_species_indexes[i]]);
            gradient_species_.push_back(&all_species_[gradient_species_indexes[i]]);
        }

        const IndexVector &reactive_species_indexes = diffusion_reaction_material_.ReactiveSpeciesIndexes();
        for (size_t i = 0; i != reactive_species_indexes.size(); ++i)
        {
            reactive_species_.push_back(&all_species_[reactive_species_indexes[i]]);
        }
    };
    virtual ~DiffusionReactionParticles(){};

    StdVec<std::string> &AllSpeciesNames() { return all_species_names_; };
    std::map<std::string, size_t> AllSpeciesIndexMap() { return species_indexes_map_; };
    StdVec<StdLargeVec<Real> *> &DiffusionSpecies() { return diffusion_species_; };
    StdVec<StdLargeVec<Real> *> &GradientSpecies() { return gradient_species_; };
    StdVec<StdLargeVec<Real> *> &ReactiveSpecies() { return reactive_species_; };

    virtual void initializeOtherVariables() override
    {
        BaseParticlesType::initializeOtherVariables();

        std::map<std::string, size_t>::iterator itr;
        for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr)
        {
            // Register a specie.
            this->registerVariable(all_species_[itr->second], itr->first);
            /** the scalars will be sorted if particle sorting is called, Note that we call a template function from a template class. */
            this->template registerSortableVariable<Real>(itr->first);
            /** add species to basic output particle data. */
            this->template addVariableToWrite<Real>(itr->first);
        }
    };

    virtual DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType> *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // DIFFUSION_REACTION_PARTICLES_H