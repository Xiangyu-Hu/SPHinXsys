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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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

#include "base_data_type_package.h"
#include "particle_functors.h"

#include <functional>
#include <map>
using namespace std::placeholders;

namespace SPH
{
class BaseParticles;

class AbstractDiffusion
{
  public:
    typedef Real DataType;
    AbstractDiffusion() {};
    virtual ~AbstractDiffusion() {};
    virtual StdVec<AbstractDiffusion *> AllDiffusions() = 0;
    virtual Real getDiffusionTimeStepSize(Real smoothing_length) = 0;
    virtual void registerLocalParameters(BaseParticles *base_particles) {};
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) {};
    virtual void initializeLocalParameters(BaseParticles *base_particles) {};
};

class BaseDiffusion : public AbstractDiffusion
{
  public:
    BaseDiffusion(const std::string &diffusion_species_name,
                  const std::string &gradient_species_name, Real cv);
    BaseDiffusion(const std::string &species_name, Real cv);
    virtual ~BaseDiffusion() {};

    std::string DiffusionSpeciesName() { return diffusion_species_name_; };
    std::string GradientSpeciesName() { return gradient_species_name_; };
    Real getVolumeCapacity() { return cv_; };
    virtual Real getDiffusionTimeStepSize(Real smoothing_length) override;
    virtual Real getReferenceDiffusivity() = 0;
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) = 0;
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) = 0;
    virtual StdVec<AbstractDiffusion *> AllDiffusions() override { return {this}; };

    class InverseVolumetricCapacity
    {
        Real cv1_;

      public:
        InverseVolumetricCapacity() : cv1_(0) {};
        InverseVolumetricCapacity(BaseDiffusion &encloser) : cv1_(1.0 / encloser.cv_) {};
        template <class ExecutionPolicy>
        InverseVolumetricCapacity(const ExecutionPolicy &ex_policy, BaseDiffusion &encloser)
            : InverseVolumetricCapacity(encloser){};
        Real operator()(size_t index_i) { return cv1_; };
    };

  protected:
    std::string diffusion_species_name_;
    std::string gradient_species_name_;
    Real cv_; // volumetric capacity
};

/**
 * @class IsotropicDiffusion
 * @brief isotropic diffusion property.
 */
class IsotropicDiffusion : public BaseDiffusion
{
  protected:
    Real d_coeff_; /**< diffusion coefficient. */

  public:
    IsotropicDiffusion(const std::string &diffusion_species_name,
                       const std::string &gradient_species_name,
                       Real d_coeff = 1.0, Real cv = 1.0);
    IsotropicDiffusion(const std::string &species_name, Real d_coeff = 1.0, Real cv = 1.0);
    explicit IsotropicDiffusion(ConstructArgs<std::string, Real> args);
    virtual ~IsotropicDiffusion() {};

    virtual Real getReferenceDiffusivity() override { return d_coeff_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return d_coeff_; }
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return d_coeff_;
    };

    class InterParticleDiffusionCoeff
    {
        Real d_coeff_;

      public:
        InterParticleDiffusionCoeff() : d_coeff_(0) {};
        InterParticleDiffusionCoeff(IsotropicDiffusion &encloser)
            : d_coeff_(encloser.d_coeff_) {};
        template <class ExecutionPolicy>
        InterParticleDiffusionCoeff(const ExecutionPolicy &ex_policy, IsotropicDiffusion &encloser)
            : InterParticleDiffusionCoeff(encloser){};
        Real operator()(size_t index_i, size_t index_j, const Vecd &e_ij) { return d_coeff_; };
        Real operator()(size_t index_i, size_t index_j) { return d_coeff_; };
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
    Real diff_max_; /**< maximum diffusion coefficient. */
    Real *local_diffusivity_;

  public:
    LocalIsotropicDiffusion(const std::string &diffusion_species_name,
                            const std::string &gradient_species_name,
                            Real diff_background, Real diff_max, Real cv = 1.0);
    LocalIsotropicDiffusion(const std::string &species_name, Real diff_background, Real diff_max, Real cv = 1.0);
    explicit LocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args);
    virtual ~LocalIsotropicDiffusion() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getReferenceDiffusivity() override { return diff_max_; };
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
    Real bias_d_coeff_;            /**< The bias diffusion coefficient along the fiber direction. */
    Matd transformed_diffusivity_; /**< The transformed diffusivity with inverse Cholesky decomposition. */

    void initializeDirectionalDiffusivity(Real d_coeff, Real bias_d_coeff_, Vecd bias_direction);

  public:
    DirectionalDiffusion(const std::string &diffusion_species_name,
                         const std::string &gradient_species_name,
                         Real d_coeff, Real bias_d_coeff_,
                         Vecd bias_direction, Real cv = 1.0);
    DirectionalDiffusion(const std::string &species_name,
                         Real d_coeff, Real bias_d_coeff_,
                         Vecd bias_direction, Real cv = 1.0);
    explicit DirectionalDiffusion(ConstructArgs<std::string, Real, Real, Vecd> args);
    virtual ~DirectionalDiffusion() {};

    virtual Real getReferenceDiffusivity() override
    {
        return SMAX(d_coeff_, d_coeff_ + bias_d_coeff_);
    };

    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        Vecd grad_ij = transformed_diffusivity_ * e_ij;
        return 1.0 / grad_ij.squaredNorm();
    };

    class InterParticleDiffusionCoeff
    {
        Matd transformed_diffusivity_;

      public:
        InterParticleDiffusionCoeff() : transformed_diffusivity_(ZeroData<Matd>::value) {};
        InterParticleDiffusionCoeff(DirectionalDiffusion &encloser)
            : transformed_diffusivity_(encloser.transformed_diffusivity_) {};
        template <class ExecutionPolicy>
        InterParticleDiffusionCoeff(const ExecutionPolicy &ex_policy, DirectionalDiffusion &encloser)
            : InterParticleDiffusionCoeff(encloser){};
        Real operator()(size_t index_i, size_t index_j, const Vecd &e_ij)
        {
            Vecd grad_ij = transformed_diffusivity_ * e_ij;
            return 1.0 / grad_ij.squaredNorm();
        };
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
                              Real d_coeff, Real bias_d_coeff_,
                              Vecd bias_direction, Real cv = 1.0);
    LocalDirectionalDiffusion(const std::string &species_name,
                              Real d_coeff, Real bias_d_coeff_,
                              Vecd bias_direction, Real cv = 1.0);
    virtual ~LocalDirectionalDiffusion() {};

    virtual void registerLocalParameters(BaseParticles *base_particles) override;
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) override;
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;

    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        Matd trans_diffusivity = getAverageValue(local_transformed_diffusivity_[index_i],
                                                 local_transformed_diffusivity_[index_j]);
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
    virtual ~BaseReactionModel() {};
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
class ReactionDiffusion : public AbstractDiffusion
{
  public:
    static constexpr int NumReactiveSpecies = ReactionType::NumSpecies;

  private:
    UniquePtrsKeeper<DiffusionType> diffusion_ptrs_keeper_;

  protected:
    ReactionType *reaction_model_;
    StdVec<AbstractDiffusion *> all_diffusions_;

  public:
    /** Constructor for material with diffusion and reaction. */
    ReactionDiffusion(ReactionType *reaction_model)
        : AbstractDiffusion(), reaction_model_(reaction_model) {};
    virtual ~ReactionDiffusion() {};
    ReactionType &ReactionModel() { return *reaction_model_; };
    virtual StdVec<AbstractDiffusion *> AllDiffusions() override { return all_diffusions_; };

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
    virtual Real getDiffusionTimeStepSize(Real smoothing_length) override
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
        auto species_names = reaction_model_->getSpeciesNames();

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