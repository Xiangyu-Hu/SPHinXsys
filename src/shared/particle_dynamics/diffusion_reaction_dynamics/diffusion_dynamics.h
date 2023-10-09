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

template <typename... InteractionTypes>
class DiffusionRelaxation;

template <template <typename... T> class DataDelegationType, class ParticlesType, class... ContactParticlesType>
class DiffusionRelaxation<Base, DataDelegationType<ParticlesType, ContactParticlesType...>>
    : public LocalDynamics,
      public DataDelegationType<ParticlesType, ContactParticlesType...>
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

    template <class DynamicsIdentifier>
    explicit DiffusionRelaxation(DynamicsIdentifier &identifier);
    virtual ~DiffusionRelaxation(){};
    StdVec<BaseDiffusion *> &AllDiffusions() { return material_.AllDiffusions(); };
    void update(size_t index_i, Real dt = 0.0);
};

class KernelGradientInner
{
  public:
    explicit KernelGradientInner(BaseParticles *inner_particles){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return dW_ijV_j * e_ij;
    };
};

class CorrectedKernelGradientInner
{
    StdLargeVec<Matd> &B_;

  public:
    explicit CorrectedKernelGradientInner(BaseParticles *inner_particles)
        : B_(*inner_particles->getVariableByName<Matd>("KernelCorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + B_[index_j]) * e_ij;
    };
};

/**
 * @class DiffusionRelaxationInner
 * @brief Compute the diffusion relaxation process of all species
 */
template <class ParticlesType, class KernelGradientType>
class DiffusionRelaxation<Inner<ParticlesType, KernelGradientType>>
    : public DiffusionRelaxation<Base, DiffusionReactionInnerData<ParticlesType>>
{
  protected:
    KernelGradientType kernel_gradient_;

  public:
    typedef BaseInnerRelation BodyRelationType;
    explicit DiffusionRelaxation(BaseInnerRelation &inner_relation);
    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
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
        : B_(*inner_particles->getVariableByName<Matd>("KernelCorrectionMatrix")),
          contact_B_(*contact_particles->getVariableByName<Matd>("KernelCorrectionMatrix")){};
    Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
    {
        return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
    };
};

template <class ParticlesType, class ContactParticlesType, class ContactKernelGradientType>
class DiffusionRelaxation<BaseContact, ParticlesType, ContactParticlesType, ContactKernelGradientType>
    : public DiffusionRelaxation<Base, DiffusionReactionContactData<ParticlesType, ContactParticlesType>>
{
  protected:
    StdVec<StdVec<std::string>> contact_gradient_species_names_;
    StdVec<ContactKernelGradientType> contact_kernel_gradients_;

  public:
    typedef BaseContactRelation BodyRelationType;

    explicit DiffusionRelaxation(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxation(){};
};

template <typename... ControlTypes>
class Dirichlet; /**< Contact interaction with Dirichlet boundary condition */

template <typename... ContactParameters>
class DiffusionRelaxation<Dirichlet<ContactParameters...>>
    : public DiffusionRelaxation<BaseContact, ContactParameters...>
{
  protected:
    StdVec<StdVec<StdLargeVec<Real> *>> contact_gradient_species_;
    void getDiffusionChangeRateDirichlet(
        size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij,
        const StdVec<StdLargeVec<Real> *> &gradient_species_k);

  public:
    explicit DiffusionRelaxation(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxation(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
};

template <typename... ControlTypes>
class Neumann; /**< Contact interaction with Neumann boundary condition */

template <typename... ContactParameters>
class DiffusionRelaxation<Neumann<ContactParameters...>>
    : public DiffusionRelaxation<BaseContact, ContactParameters...>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Real> *> contact_heat_flux_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;

  protected:
    void getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j,
                                       Real surface_area_ij_Neumann, StdLargeVec<Real> &heat_flux_k);

  public:
    explicit DiffusionRelaxation(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxation(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <typename... ControlTypes>
class Robin; /**< Contact interaction with Robin boundary condition */

template <typename... ContactParameters>
class DiffusionRelaxation<Robin<ContactParameters...>>
    : public DiffusionRelaxation<BaseContact, ContactParameters...>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Real> *> contact_convection_;
    StdVec<Real *> contact_T_infinity_;
    StdVec<StdLargeVec<Vecd> *> contact_n_;

  protected:
    void getDiffusionChangeRateRobin(
        size_t particle_i, size_t particle_j, Real surface_area_ij_Robin,
        StdLargeVec<Real> &convection_k, Real &T_infinity_k);

  public:
    explicit DiffusionRelaxation(BaseContactRelation &contact_relation);
    virtual ~DiffusionRelaxation(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class InitializationRK
 * @brief Initialization of a runge-kutta integration scheme.
 */
template <class ParticlesType>
class InitializationRK : public DiffusionRelaxation<Base, DiffusionReactionSimpleData<ParticlesType>>
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

  public:
    template <typename... ContactArgsType>
    SecondStageRK2(typename FirstStageType::BodyRelationType &body_relation,
                   StdVec<StdLargeVec<Real>> &diffusion_species_s, ContactArgsType &&...contact_args)
        : FirstStageType(body_relation, std::forward<ContactArgsType>(contact_args)...),
          diffusion_species_s_(diffusion_species_s){};
    virtual ~SecondStageRK2(){};

    void update(size_t index_i, Real dt = 0.0);
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