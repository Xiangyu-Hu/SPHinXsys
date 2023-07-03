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
 * @file 	diffusion_reaction.h
 * @brief 	Describe the diffusive and reaction in which
 *          the dynamics is characterized by diffusion equation and reactive source terms.
 *			Typical physical processes are diffusion, heat conduction
 *			and chemical and biological reactions.
 */

#ifndef DIFFUSION_REACTION_H
#define DIFFUSION_REACTION_H

#include "base_material.h"

#include <functional>
#include <map>
using namespace std::placeholders;

namespace SPH
{
/**
 * @class BaseDiffusion
 * @brief diffusion property abstract base class.
 */
class BaseDiffusion : public BaseMaterial
{
  public:
    BaseDiffusion(size_t diffusion_species_index, size_t gradient_species_index)
        : BaseMaterial(), diffusion_species_index_(diffusion_species_index),
          gradient_species_index_(gradient_species_index)
    {
        material_type_name_ = "BaseDiffusion";
    };
    virtual ~BaseDiffusion(){};

    size_t diffusion_species_index_;
    size_t gradient_species_index_;

    virtual Real getReferenceDiffusivity() = 0;
    virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd &direction_from_j_to_i) = 0;
};

/**
 * @class IsotropicDiffusion
 * @brief isotropic diffusion property.
 */
class IsotropicDiffusion : public BaseDiffusion
{
  protected:
    Real diff_cf_; /**< diffusion coefficient. */

  public:
    IsotropicDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
                       Real diff_cf = 1.0)
        : BaseDiffusion(diffusion_species_index, gradient_species_index),
          diff_cf_(diff_cf)
    {
        material_type_name_ = "IsotropicDiffusion";
    };
    virtual ~IsotropicDiffusion(){};

    virtual Real getReferenceDiffusivity() override { return diff_cf_; };
    virtual Real getInterParticleDiffusionCoff(size_t particle_i, size_t particle_j, Vecd &direction_from_j_to_i) override
    {
        return diff_cf_;
    };
};

/**
 * @class DirectionalDiffusion
 * @brief Diffusion is biased along a specific direction.
 */
class DirectionalDiffusion : public IsotropicDiffusion
{
  protected:
    Vecd bias_direction_;          /**< Reference bias direction. */
    Real bias_diff_cf_;            /**< The bias diffusion coefficient along the fiber direction. */
    Matd transformed_diffusivity_; /**< The transformed diffusivity with inverse Cholesky decomposition. */

    void initializeDirectionalDiffusivity(Real diff_cf, Real bias_diff_cf, Vecd bias_direction);

  public:
    DirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
                         Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
        : IsotropicDiffusion(diffusion_species_index, gradient_species_index, diff_cf),
          bias_direction_(bias_direction), bias_diff_cf_(bias_diff_cf),
          transformed_diffusivity_(Matd::Identity())
    {
        material_type_name_ = "DirectionalDiffusion";
        initializeDirectionalDiffusivity(diff_cf, bias_diff_cf, bias_direction);
    };
    virtual ~DirectionalDiffusion(){};

    virtual Real getReferenceDiffusivity() override
    {
        return SMAX(diff_cf_, diff_cf_ + bias_diff_cf_);
    };

    virtual Real getInterParticleDiffusionCoff(size_t particle_index_i,
                                               size_t particle_index_j, Vecd &inter_particle_direction) override
    {
        Vecd grad_ij = transformed_diffusivity_ * inter_particle_direction;
        return 1.0 / grad_ij.squaredNorm();
    };
};

/**
 * @class LocalDirectionalDiffusion
 * @brief Diffusion is biased along a specific direction.
 */
class LocalDirectionalDiffusion : public DirectionalDiffusion
{
  protected:
    StdLargeVec<Vecd> local_bias_direction_;
    StdLargeVec<Matd> local_transformed_diffusivity_;

  public:
    LocalDirectionalDiffusion(size_t diffusion_species_index, size_t gradient_species_index,
                              Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
        : DirectionalDiffusion(diffusion_species_index, gradient_species_index, diff_cf, bias_diff_cf, bias_direction)
    {
        material_type_name_ = "LocalDirectionalDiffusion";
    };
    virtual ~LocalDirectionalDiffusion(){};

    virtual void registerReloadLocalParameters(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getInterParticleDiffusionCoff(size_t particle_index_i, size_t particle_index_j, Vecd &inter_particle_direction) override
    {
        Matd trans_diffusivity = getAverageValue(local_transformed_diffusivity_[particle_index_i], local_transformed_diffusivity_[particle_index_j]);
        Vecd grad_ij = trans_diffusivity * inter_particle_direction;
        return 1.0 / grad_ij.squaredNorm();
    };
};

/**
 * @class BaseReactionModel
 * @brief Base class for all reaction models.
 */
template <int NUM_SPECIES>
class BaseReactionModel
{
  public:
    static constexpr int NumSpecies = NUM_SPECIES;
    typedef std::array<Real, NUM_SPECIES> LocalSpecies;
    typedef std::array<std::string, NUM_SPECIES> SpeciesNames;
    typedef std::function<Real(LocalSpecies &)> ReactionFunctor;
    StdVec<ReactionFunctor> get_production_rates_;
    StdVec<ReactionFunctor> get_loss_rates_;

    BaseReactionModel()
    {
        std::cout << "\n Error: default constructor for non-empty reaction model!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    };

    explicit BaseReactionModel(const SpeciesNames &species_names)
        : species_names_(species_names)
    {
        if constexpr (NUM_SPECIES == 0)
        {
            reaction_model_ = "EmptyReactionModel";
        }
        else
        {
            reaction_model_ = "BaseReactionModel";
            for (size_t i = 0; i != species_names.size(); ++i)
            {
                species_indexes_map_.insert(make_pair(species_names[i], i));
            }
        }
    };
    virtual ~BaseReactionModel(){};
    SpeciesNames &getSpeciesNames() { return species_names_; };

  protected:
    std::string reaction_model_;
    SpeciesNames species_names_;
    std::map<std::string, size_t> species_indexes_map_;
};
/** explicit specialization for empty reaction model */
template <>
inline BaseReactionModel<0>::BaseReactionModel() : reaction_model_("EmptyReactionModel"){};
using NoReaction = BaseReactionModel<0>;

/**
 * @class DiffusionReaction
 * @brief Complex material for diffusion or/and reactions.
 */
template <class BaseMaterialType = BaseMaterial, int NUM_REACTIVE_SPECIES = 0>
class DiffusionReaction : public BaseMaterialType
{
  public:
    static constexpr int NumReactiveSpecies = NUM_REACTIVE_SPECIES;

  private:
    UniquePtrsKeeper<BaseDiffusion> diffusion_ptr_keeper_;
    SharedPtrKeeper<BaseReactionModel<NUM_REACTIVE_SPECIES>> reaction_ptr_keeper_;

  protected:
    typedef std::array<std::string, NUM_REACTIVE_SPECIES> ReactiveSpeciesNames;
    StdVec<std::string> all_species_names_;
    BaseReactionModel<NUM_REACTIVE_SPECIES> &reaction_model_;
    std::map<std::string, size_t> all_species_indexes_map_;
    StdVec<BaseDiffusion *> all_diffusions_;
    IndexVector reactive_species_indexes_;
    IndexVector diffusion_species_indexes_;
    IndexVector gradient_species_indexes_;

  public:
    /** Constructor for material with diffusion and reaction. */
    template <typename... MaterialArgs>
    DiffusionReaction(const StdVec<std::string> &all_species_names,
                      SharedPtr<BaseReactionModel<NUM_REACTIVE_SPECIES>> reaction_model_ptr,
                      MaterialArgs &&...material_args)
        : BaseMaterialType(std::forward<MaterialArgs>(material_args)...),
          all_species_names_(all_species_names),
          reaction_model_(reaction_ptr_keeper_.assignRef(reaction_model_ptr))
    {
        BaseMaterialType::material_type_name_ =
            NUM_REACTIVE_SPECIES == 0 ? "Diffusion" : "DiffusionReaction";

        for (size_t i = 0; i != all_species_names.size(); ++i)
        {
            all_species_indexes_map_.insert(make_pair(all_species_names[i], i));
        }

        for (size_t i = 0; i != NUM_REACTIVE_SPECIES; ++i)
        {
            const ReactiveSpeciesNames &reactive_species_names = reaction_model_.getSpeciesNames();
            size_t reactive_species_index = all_species_indexes_map_[reactive_species_names[i]];
            if (reactive_species_index != all_species_indexes_map_.size())
            {
                reactive_species_indexes_.push_back(reactive_species_index);
            }
            else
            {
                std::cout << "\n Error: reactive species '" << reactive_species_names[i] << "' not defined!" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
        }
    };
    virtual ~DiffusionReaction(){};
    StdVec<std::string> &AllSpeciesNames() { return all_species_names_; };
    std::map<std::string, size_t> AllSpeciesIndexMap() { return all_species_indexes_map_; };
    IndexVector &ReactiveSpeciesIndexes() { return reactive_species_indexes_; };
    IndexVector &DiffusionSpeciesIndexes() { return diffusion_species_indexes_; };
    IndexVector &GradientSpeciesIndexes() { return gradient_species_indexes_; };
    StdVec<BaseDiffusion *> &AllDiffusions() { return all_diffusions_; };
    BaseReactionModel<NUM_REACTIVE_SPECIES> &ReactionModel() { return reaction_model_; };

    virtual void registerReloadLocalParameters(BaseParticles *base_particles) override
    {
        BaseMaterialType::registerReloadLocalParameters(base_particles);
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            all_diffusions_[k]->registerReloadLocalParameters(base_particles);
    };

    virtual void initializeLocalParameters(BaseParticles *base_particles) override
    {
        BaseMaterialType::initializeLocalParameters(base_particles);
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            all_diffusions_[k]->initializeLocalParameters(base_particles);
    };

    /**
     * @brief Get diffusion time step size. Here, I follow the reference:
     * https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf
     */
    Real getDiffusionTimeStepSize(Real smoothing_length)
    {
        Real diff_coff_max = 0.0;
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            diff_coff_max = SMAX(diff_coff_max, all_diffusions_[k]->getReferenceDiffusivity());
        return 0.5 * smoothing_length * smoothing_length / diff_coff_max / Real(Dimensions);
    };

    /** Initialize a diffusion material. */
    template <class DiffusionType, typename... ConstructorArgs>
    void initializeAnDiffusion(const std::string &diffusion_species_name,
                               const std::string &gradient_species_name, ConstructorArgs &&...args)
    {
        size_t diffusion_species_index = all_species_indexes_map_[diffusion_species_name];
        size_t gradient_species_index = all_species_indexes_map_[gradient_species_name];
        diffusion_species_indexes_.push_back(diffusion_species_index);
        gradient_species_indexes_.push_back(gradient_species_index);

        all_diffusions_.push_back(
            diffusion_ptr_keeper_.createPtr<DiffusionType>(
                diffusion_species_index, gradient_species_index, std::forward<ConstructorArgs>(args)...));
    };

    virtual DiffusionReaction<BaseMaterialType, NUM_REACTIVE_SPECIES> *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // DIFFUSION_REACTION_H