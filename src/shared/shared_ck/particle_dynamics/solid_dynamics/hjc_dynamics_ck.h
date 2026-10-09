// SPDX-License-Identifier: Apache-2.0
#ifndef HJC_DYNAMICS_CK_H
#define HJC_DYNAMICS_CK_H

#include "hjc_solid.hpp"
#include "structure_dynamics.hpp"

namespace SPH
{
namespace solid_dynamics
{
/** Incremental HJC stress update with the total Lagrangian pair force. */
template <typename...>
class HJCIntegration1stHalfCK;

template <typename... Parameters>
class HJCIntegration1stHalfCK<Inner<OneLevel, Parameters...>>
    : public Interaction<Inner<Parameters...>>, public BaseStructureIntegration1stHalf
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit HJCIntegration1stHalfCK(DynamicsIdentifier &identifier);

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t i, Real dt);

      protected:
        HJCSolid::ConstituteKernel constitute_;
        Real rho0_;
        Real *rho_;
        Vecd *pos_, *vel_;
        Matd *B_, *F_, *dF_dt_, *stress_on_particle_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t i, Real dt);

      protected:
        Real damping_scale_;
        Real *Vol0_;
        Vecd *pos_, *vel_, *force_;
        Matd *F_, *stress_on_particle_;
    };

  protected:
    HJCSolid &material_;
    Real h_ref_;
};

/** Includes the evolving EOS tangent and the same acceleration bound as CPU HJC. */
class HJCAcousticTimeStepCK : public LocalDynamicsReduce<ReduceMin<Real>>
{
  public:
    explicit HJCAcousticTimeStepCK(SPHBody &body, Real acousticCFL = Real(0.2));

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Real reduce(size_t i, Real dt = 0);

      protected:
        Real h_ref_, acousticCFL_;
        Real *rho_, *mass_, *modulus_;
        Vecd *vel_, *force_, *force_prior_;
    };

  protected:
    Real h_ref_, acousticCFL_;
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_modulus_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_force_, *dv_force_prior_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif
