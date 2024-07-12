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
 * @file diffusion_reaction.h
 * @brief Describe the diffusive and reaction in which
 * the dynamics is characterized by diffusion equation and reactive source terms.
 * Typical physical processes are diffusion, heat conduction
 * and chemical and biological reactions.
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
    BaseDiffusion(const std::string &diffusion_species_name,
                  const std::string &gradient_species_name);
    BaseDiffusion(const std::string &species_name);
    virtual ~BaseDiffusion(){};

    std::string DiffusionSpeciesName() { return diffusion_species_name_; };
    std::string GradientSpeciesName() { return gradient_species_name_; };
    Real getDiffusionTimeStepSize(Real smoothing_length);
    virtual Real getReferenceDiffusivity() = 0;
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) = 0;
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) = 0;

  protected:
    std::string diffusion_species_name_;
    std::string gradient_species_name_;
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
    IsotropicDiffusion(const std::string &diffusion_species_name,
                       const std::string &gradient_species_name,
                       Real diff_cf = 1.0);
    IsotropicDiffusion(const std::string &species_name, Real diff_cf = 1.0);
    virtual ~IsotropicDiffusion(){};

    virtual Real getReferenceDiffusivity() override { return diff_cf_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return diff_cf_; }
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return diff_cf_;
    };
};

/**
 * @class LocalIsotropicDiffusion
 * @brief diffusion coefficient is locally different (k is not uniformly distributed).
 * TODO: The difference between algebraic and geometric average should be identified.
 */
class LocalIsotropicDiffusion : public IsotropicDiffusion
{
  protected:
    Real *local_diffusivity_;

  public:
    LocalIsotropicDiffusion(const std::string &diffusion_species_name,
                            const std::string &gradient_species_name,
                            Real diff_cf = 1.0);
    LocalIsotropicDiffusion(const std::string &species_name, Real diff_cf = 1.0);
    virtual ~LocalIsotropicDiffusion(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getReferenceDiffusivity() override { return diff_cf_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return local_diffusivity_[index_i]; };
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
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
    DirectionalDiffusion(const std::string &diffusion_species_name,
                         const std::string &gradient_species_name,
                         Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
    DirectionalDiffusion(const std::string &species_name,
                         Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
    virtual ~DirectionalDiffusion(){};

    virtual Real getReferenceDiffusivity() override
    {
        return SMAX(diff_cf_, diff_cf_ + bias_diff_cf_);
    };

    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        Vecd grad_ij = transformed_diffusivity_ * e_ij;
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
    Vecd *local_bias_direction_;
    Matd *local_transformed_diffusivity_;

  public:
    LocalDirectionalDiffusion(const std::string &diffusion_species_name,
                              const std::string &gradient_species_name,
                              Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
    LocalDirectionalDiffusion(const std::string &species_name,
                              Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
    virtual ~LocalDirectionalDiffusion(){};

    virtual void registerLocalParameters(BaseParticles *base_particles) override;
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        Matd trans_diffusivity = getAverageValue(local_transformed_diffusivity_[index_i], local_transformed_diffusivity_[index_j]);
        Vecd grad_ij = trans_diffusivity * e_ij;
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

/**
 * @class ReactionDiffusion
 * @brief Complex material for reaction and diffusion.
 */
template <class ReactionType, class DiffusionType>
class ReactionDiffusion : public BaseMaterial
{
  public:
    static constexpr int NumReactiveSpecies = ReactionType::NumSpecies;

  private:
    UniquePtrsKeeper<DiffusionType> diffusion_ptrs_keeper_;

  protected:
    ReactionType &reaction_model_;
    StdVec<DiffusionType *> all_diffusions_;

  public:
    /** Constructor for material with diffusion and reaction. */
    ReactionDiffusion(ReactionType &reaction_model)
        : BaseMaterial(), reaction_model_(reaction_model)
    {
        material_type_name_ = "ReactionDiffusion";
    };
    virtual ~ReactionDiffusion(){};
    StdVec<DiffusionType *> AllDiffusions() { return all_diffusions_; };
    ReactionType &ReactionModel() { return reaction_model_; };

    virtual void registerLocalParameters(BaseParticles *base_particles) override
    {
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            all_diffusions_[k]->registerLocalParameters(base_particles);
    };

    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override
    {
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            all_diffusions_[k]->registerLocalParametersFromReload(base_particles);
    };

    virtual void initializeLocalParameters(BaseParticles *base_particles) override
    {
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            all_diffusions_[k]->initializeLocalParameters(base_particles);
    };

    /**
     * @brief Get diffusion time step size. Here, I follow the reference:
     * https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf
     */
    Real getDiffusionTimeStepSize(Real smoothing_length)
    {
        Real dt = MaxReal;
        for (size_t k = 0; k < all_diffusions_.size(); ++k)
            dt = SMIN(dt, all_diffusions_[k]->getDiffusionTimeStepSize(smoothing_length));
        return dt;
    };

    template <typename... Args>
    void addDiffusion(const std::string &diffusion_species_name,
                      const std::string &gradient_species_name, Args &&...args)
    {
        auto species_names = reaction_model_.getSpeciesNames();

        if (std::find(species_names.begin(), species_names.end(), diffusion_species_name) != std::end(species_names) &&
            std::find(species_names.begin(), species_names.end(), gradient_species_name) != std::end(species_names))
        {
            all_diffusions_.push_back(
                diffusion_ptrs_keeper_.template createPtr<DiffusionType>(
                    diffusion_species_name, gradient_species_name, std::forward<Args>(args)...));
        }
        else
        {
            std::cout << "\n Error: diffusion species '" << diffusion_species_name
                      << "' or gradient species '" << gradient_species_name
                      << "' not defined!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    };
};
} // namespace SPH
#endif // DIFFUSION_REACTION_H