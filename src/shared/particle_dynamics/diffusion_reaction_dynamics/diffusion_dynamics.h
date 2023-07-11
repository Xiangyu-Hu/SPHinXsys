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
 * @file    diffusion_dynamics.h
 * @brief   These are particle dynamics applicable for all type of particles.
 * @author  Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_H
#define DIFFUSION_DYNAMICS_H

#include "general_diffusion_reaction_dynamics.h"

namespace SPH
{
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
 * @class BaseDiffusionRelaxation
 * @brief Base for compute the diffusion of all species
 */
template <class ParticlesType>
class BaseDiffusionRelaxation
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  protected:
    typedef typename ParticlesType::DiffusionReactionMaterial Material;
    Material &material_;
    StdVec<BaseDiffusion *> &all_diffusions_;
    StdVec<StdLargeVec<Real> *> &diffusion_species_;
    StdVec<StdLargeVec<Real> *> &gradient_species_;
    StdVec<StdLargeVec<Real> *> diffusion_dt_;

  public:
    typedef ParticlesType InnerParticlesType;
    explicit BaseDiffusionRelaxation(SPHBody &sph_body);
    virtual ~BaseDiffusionRelaxation(){};
    StdVec<BaseDiffusion *> &AllDiffusions() { return material_.AllDiffusions(); };
};

class KernelGradientInner
{
  public:
    KernelGradientInner(BaseParticles *inner_particles){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return dW_ijV_j * e_ij;
    };
};

class CorrectedKernelGradientInner
{
    StdLargeVec<Matd> &B_;

  public:
    CorrectedKernelGradientInner(BaseParticles *inner_particles)
        : B_(*inner_particles->getVariableByName<Matd>("CorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + B_[index_j]) * e_ij;
    };
};

/**
 * @class DiffusionRelaxationInner
 * @brief Compute the diffusion relaxation process of all species
 */
template <class ParticlesType, class KernelGradientType = KernelGradientInner>
class DiffusionRelaxationInner
    : public BaseDiffusionRelaxation<ParticlesType>,
      public DataDelegateInner<ParticlesType, DataDelegateEmptyBase>
{
  protected:
    KernelGradientType kernel_gradient_;
    void initializeDiffusionChangeRate(size_t particle_i);
    void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij);
    virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

  public:
    typedef BaseInnerRelation BodyRelationType;
    explicit DiffusionRelaxationInner(BaseInnerRelation &inner_relation);
    virtual ~DiffusionRelaxationInner(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class KernelGradientContact
{
  public:
    KernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return dW_ijV_j * e_ij;
    };
};

class CorrectedKernelGradientContact
{
    StdLargeVec<Matd> &B_;
    StdLargeVec<Matd> &contact_B_;

  public:
    CorrectedKernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles)
        : B_(*inner_particles->getVariableByName<Matd>("CorrectionMatrix")),
          contact_B_(*contact_particles->getVariableByName<Matd>("CorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
    };
};

/**
 * @class DiffusionRelaxationContact
 * @brief Base class for diffusion relaxation process between two contact bodies.
 */
template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact>
class BaseDiffusionRelaxationContact
    : public BaseDiffusionRelaxation<ParticlesType>,
      public DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>
{
  protected:
    StdVec<StdVec<std::string>> contact_gradient_species_names_;
    StdVec<KernelGradientType> contact_kernel_gradients_;

  public:
    typedef BaseContactRelation BodyRelationType;

    explicit BaseDiffusionRelaxationContact(BaseContactRelation &contact_relation);
    virtual ~BaseDiffusionRelaxationContact(){};
};

/**
 * @class DiffusionRelaxationDirichlet
 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
 */
template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact>
class DiffusionRelaxationDirichlet
    : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>
{
  protected:
    StdVec<StdVec<StdLargeVec<Real> *>> contact_gradient_species_;
    void getDiffusionChangeRateDirichletContact(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<StdLargeVec<Real> *> &gradient_species_k);

  public:
    explicit DiffusionRelaxationDirichlet(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxationDirichlet(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class DiffusionRelaxationNeumann
 * @brief Contact diffusion relaxation with Neumann boundary condition.
 */
template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact>
class DiffusionRelaxationNeumann
    : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Real> *> contact_heat_flux_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;

  protected:
    void getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j,
                                              Real surface_area_ij_Neumann, StdLargeVec<Real> &heat_flux_k);

  public:
    explicit DiffusionRelaxationNeumann(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxationNeumann(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class RelaxationOfAllDiffusionSpeciesRobinContact
 * @brief Contact diffusion relaxation with Robin boundary condition.
 */
template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact>
class DiffusionRelaxationRobin
    : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Real> *> contact_convection_;
    StdVec<Real *> contact_T_infinity_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;

  protected:
    void getDiffusionChangeRateRobinContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, StdLargeVec<Real> &convection_k, Real &T_infinity_k);

  public:
    explicit DiffusionRelaxationRobin(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxationRobin(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class InitializationRK
 * @brief Initialization of a runge-kutta integration scheme.
 */
template <class ParticlesType>
class InitializationRK : public BaseDiffusionRelaxation<ParticlesType>
{
  protected:
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

/**
 * @class DiffusionRelaxationRK2
 * @brief The 2nd-order runge-kutta integration scheme.
 * A intermediate state for species is introduced here to achieve multi-step integration.
 */
template <class FirstStageType>
class DiffusionRelaxationRK2 : public BaseDynamics<void>
{
  protected:
    StdVec<StdLargeVec<Real>> diffusion_species_s_; /**< Intermediate state */
    SimpleDynamics<InitializationRK<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
    InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
    InteractionWithUpdate<SecondStageRK2<FirstStageType>> rk2_2nd_stage_;
    StdVec<BaseDiffusion *> all_diffusions_;

  public:
    template <typename... ContactArgsType>
    explicit DiffusionRelaxationRK2(typename FirstStageType::BodyRelationType &body_relation,
                                    ContactArgsType &&...contact_args);
    virtual ~DiffusionRelaxationRK2(){};

    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_H