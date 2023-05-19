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
 * @file 	particle_dynamics_diffusion_reaction.h
 * @brief 	This is the particle dynamics applicable for all type bodies
 * 			TODO: there is an issue on applying corrected configuration for contact bodies..
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_H
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_H

#include "all_particle_dynamics.h"
#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"

namespace SPH
{
template <class ParticlesType>
using DiffusionReactionSimpleData = DataDelegateSimple<ParticlesType>;

template <class ParticlesType>
using DiffusionReactionInnerData = DataDelegateInner<ParticlesType>;

template <class ParticlesType, class ContactParticlesType>
using DiffusionReactionContactData =
    DataDelegateContact<ParticlesType, ContactParticlesType>;

/**
 * @class DiffusionReactionInitialCondition
 * @brief Pure abstract class for initial conditions
 */
template <class ParticlesType>
class DiffusionReactionInitialCondition
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionReactionInitialCondition(SPHBody &sph_body);
    virtual ~DiffusionReactionInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

/**
 * @class GetDiffusionTimeStepSize
 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
 */
template <class ParticlesType>
class GetDiffusionTimeStepSize
    : public BaseDynamics<Real>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit GetDiffusionTimeStepSize(SPHBody &sph_body);
    virtual ~GetDiffusionTimeStepSize(){};

    virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };

  protected:
    Real diff_time_step_;
};

/**
 * @class DiffusionRelaxationInner
 * @brief Compute the diffusion relaxation process of all species
 */
template <class ParticlesType>
class DiffusionRelaxationInner
    : public LocalDynamics,
      public DiffusionReactionInnerData<ParticlesType>
{
  protected:
    typedef typename ParticlesType::DiffusionReactionMaterial Material;
    Material &material_;
    StdVec<BaseDiffusion *> &all_diffusions_;
    StdVec<StdLargeVec<Real> *> &diffusion_species_;
    StdVec<StdLargeVec<Real> *> &gradient_species_;
    StdVec<StdLargeVec<Real> *> diffusion_dt_;

    void initializeDiffusionChangeRate(size_t particle_i);
    void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij);
    virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

  public:
    typedef ParticlesType InnerParticlesType;
    typedef BaseInnerRelation BodyRelationType;

    explicit DiffusionRelaxationInner(BaseInnerRelation &inner_relation);
    virtual ~DiffusionRelaxationInner(){};
    StdVec<BaseDiffusion *> &AllDiffusions() { return material_.AllDiffusions(); };
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class DiffusionRelaxationContact
 * @brief Base class for diffusion relaxation process between two contact bodies.
 */
template <class ParticlesType, class ContactParticlesType>
class DiffusionRelaxationContact
    : public LocalDynamics,
      public DiffusionReactionContactData<ParticlesType, ContactParticlesType>
{
  protected:
    typedef typename ParticlesType::DiffusionReactionMaterial Material;
    Material &material_;
    StdVec<BaseDiffusion *> &all_diffusions_;
    StdVec<StdLargeVec<Real> *> &diffusion_species_;
    StdVec<StdLargeVec<Real> *> &gradient_species_;
    StdVec<StdLargeVec<Real> *> diffusion_dt_;

  public:
    StdVec<StdVec<StdLargeVec<Real> *>> contact_gradient_species_;

    typedef ParticlesType InnerParticlesType;
    typedef BaseContactRelation BodyRelationType;

    explicit DiffusionRelaxationContact(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxationContact(){};
    StdVec<BaseDiffusion *> &AllDiffusions() { return material_.AllDiffusions(); };
};

/**
 * @class DiffusionRelaxationDirichlet
 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
 */
template <class ParticlesType, class ContactParticlesType>
class DiffusionRelaxationDirichlet
    : public DiffusionRelaxationContact<ParticlesType, ContactParticlesType>
{
  protected:
    void getDiffusionChangeRateDirichletContact(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
                                                const StdVec<StdLargeVec<Real> *> &gradient_species_k);

  public:
    explicit DiffusionRelaxationDirichlet(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxationDirichlet(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class InitializationRK
 * @brief Initialization of a runge-kutta integration scheme.
 */
template <class ParticlesType>
class InitializationRK : public LocalDynamics,
                         public DiffusionReactionSimpleData<ParticlesType>
{
  protected:
    typename ParticlesType::DiffusionReactionMaterial &material_;
    StdVec<BaseDiffusion *> &all_diffusions_;
    StdVec<StdLargeVec<Real> *> &diffusion_species_;
    StdVec<StdLargeVec<Real>> &diffusion_species_s_;

  public:
    InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &diffusion_species_s);
    virtual ~InitializationRK(){};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class SecondStageRK2
 * @brief The second stage of the 2nd-order Runge-Kutta scheme.
 */
template <class FirstStageType>
class SecondStageRK2 : public FirstStageType
{
  protected:
    StdVec<StdLargeVec<Real>> &diffusion_species_s_;
    virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

  public:
    template <typename... ContactArgsType>
    SecondStageRK2(typename FirstStageType::BodyRelationType &body_relation,
                   StdVec<StdLargeVec<Real>> &diffusion_species_s, ContactArgsType &&...contact_args)
        : FirstStageType(body_relation, std::forward<ContactArgsType>(contact_args)...),
          diffusion_species_s_(diffusion_species_s){};
    virtual ~SecondStageRK2(){};
};

template <class FirstStageType>
class DiffusionRelaxationRK2 : public BaseDynamics<void>
{
  protected:
    /** Intermediate Value */
    StdVec<StdLargeVec<Real>> diffusion_species_s_;
    SimpleDynamics<InitializationRK<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
    InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
    InteractionWithUpdate<SecondStageRK2<FirstStageType>> rk2_2nd_stage_;
    StdVec<BaseDiffusion *> all_diffusions_;

  public:
    template <typename... ContactArgsType>
    explicit DiffusionRelaxationRK2(typename FirstStageType::BodyRelationType &body_relation, ContactArgsType &&...contact_args)
        : BaseDynamics<void>(body_relation.getSPHBody()),
          rk2_initialization_(body_relation.getSPHBody(), diffusion_species_s_),
          rk2_1st_stage_(body_relation, std::forward<ContactArgsType>(contact_args)...),
          rk2_2nd_stage_(body_relation, diffusion_species_s_, std::forward<ContactArgsType>(contact_args)...),
          all_diffusions_(rk2_1st_stage_.AllDiffusions())
    {
        diffusion_species_s_.resize(all_diffusions_.size());
        StdVec<std::string> &all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
        for (size_t i = 0; i != all_diffusions_.size(); ++i)
        {
            // Register diffusion species intermediate
            size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
            std::string &diffusion_species_name = all_species_names[diffusion_species_index];
            rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
        }
    }
    virtual ~DiffusionRelaxationRK2(){};

    virtual void exec(Real dt = 0.0) override;
};

struct UpdateAReactionSpecies
{
    Real operator()(Real input, Real production_rate, Real loss_rate, Real dt) const
    {
        return input * exp(-loss_rate * dt) + production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + TinyReal);
    };
};

/**
 * @class BaseReactionRelaxation
 * @brief Base class for computing the reaction process of all species
 */
template <class ParticlesType>
class BaseReactionRelaxation
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
    static constexpr int NumReactiveSpecies = ParticlesType::NumReactiveSpecies;
    typedef std::array<Real, NumReactiveSpecies> LocalSpecies;
    StdVec<StdLargeVec<Real> *> &reactive_species_;
    BaseReactionModel<NumReactiveSpecies> &reaction_model_;
    UpdateAReactionSpecies updateAReactionSpecies;
    void loadLocalSpecies(LocalSpecies &local_species, size_t index_i);
    void applyGlobalSpecies(LocalSpecies &local_species, size_t index_i);

  public:
    explicit BaseReactionRelaxation(SPHBody &sph_body);
    virtual ~BaseReactionRelaxation(){};

  protected:
    void advanceForwardStep(size_t index_i, Real dt);
    void advanceBackwardStep(size_t index_i, Real dt);
};

/**
 * @class ReactionRelaxationForward
 * @brief Compute the reaction process of all species by forward splitting
 */
template <class ParticlesType>
class ReactionRelaxationForward
    : public BaseReactionRelaxation<ParticlesType>
{
  public:
    ReactionRelaxationForward(SPHBody &sph_body)
        : BaseReactionRelaxation<ParticlesType>(sph_body){};
    virtual ~ReactionRelaxationForward(){};
    void update(size_t index_i, Real dt = 0.0) { this->advanceForwardStep(index_i, dt); };
};

/**
 * @class ReactionRelaxationBackward
 * @brief Compute the reaction process of all species by backward splitting
 */
template <class ParticlesType>
class ReactionRelaxationBackward
    : public BaseReactionRelaxation<ParticlesType>
{
  public:
    explicit ReactionRelaxationBackward(SPHBody &sph_body)
        : BaseReactionRelaxation<ParticlesType>(sph_body){};
    virtual ~ReactionRelaxationBackward(){};
    void update(size_t index_i, Real dt = 0.0) { this->advanceBackwardStep(index_i, dt); };
};

/**
 * @class DiffusionReactionSpeciesConstraint
 * @brief set boundary condition for diffusion problem
 */
template <class DynamicsIdentifier, class ParticlesType>
class DiffusionReactionSpeciesConstraint
    : public BaseLocalDynamics<DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    DiffusionReactionSpeciesConstraint(DynamicsIdentifier &identifier, const std::string &species_name)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()),
          phi_(this->particles_->diffusion_reaction_material_.AllSpeciesIndexMap()[species_name]),
          species_(this->particles_->all_species_[phi_]){};
    virtual ~DiffusionReactionSpeciesConstraint(){};

  protected:
    size_t phi_;
    StdLargeVec<Real> &species_;
};

/**
 * @class DiffusionBasedMapping
 * @brief Mapping inside of body according to diffusion.
 * This is a abstract class to be override for case specific implementation
 */
template <class ParticlesType>
class DiffusionBasedMapping
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionBasedMapping(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          DiffusionReactionSimpleData<ParticlesType>(sph_body),
          pos_(this->particles_->pos_), all_species_(this->particles_->all_species_){};
    virtual ~DiffusionBasedMapping(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

/**
 * @class 	SpeciesSummation
 * @brief 	Computing the total averaged parameter on the whole diffusion body.
 * 			TODO: need a test using this method
 */
template <class DynamicsIdentifier, class ParticlesType>
class SpeciesSummation
    : public BaseLocalDynamicsReduce<Real, ReduceSum<Real>, DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  protected:
    StdVec<StdLargeVec<Real>> &all_species_;
    size_t phi_;

  public:
    SpeciesSummation(DynamicsIdentifier &identifier, const std::string &species_name)
        : BaseLocalDynamicsReduce<Real, ReduceSum<Real>, DynamicsIdentifier>(identifier, Real(0)),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()),
          all_species_(this->particles_->all_species_),
          phi_(this->particles_->diffusion_reaction_material_.AllSpeciesIndexMap()[species_name])
    {
        this->quantity_name_ = "DiffusionReactionSpeciesAverage";
    };
    virtual ~SpeciesSummation(){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return all_species_[phi_][index_i];
    };
};
} // namespace SPH
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_H