#ifndef JOHNSON_COOK_DYNAMICS_CK_H
#define JOHNSON_COOK_DYNAMICS_CK_H

#include "johnson_cook_solid.hpp"
#include "structure_dynamics.hpp"

namespace SPH
{
namespace solid_dynamics
{
/** History-dependent JC stress with the classic total-Lagrangian pair force. */
template <typename...>
class JohnsonCookIntegration1stHalfCK;

template <typename... Parameters>
class JohnsonCookIntegration1stHalfCK<Inner<OneLevel, Parameters...>>
    : public Interaction<Inner<Parameters...>>, public BaseStructureIntegration1stHalf
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit JohnsonCookIntegration1stHalfCK(DynamicsIdentifier &identifier);

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t i, Real dt);

      protected:
        JohnsonCookSolid::ConstituteKernel constitute_;
        Real rho0_;
        Real *rho_;
        Vecd *pos_, *vel_;
        Matd *B_, *F_, *dF_dt_, *stress_on_particle_;
        int *active_;
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
        int *active_;
    };

  protected:
    JohnsonCookSolid &material_;
    Real h_ref_;
    DiscreteVariable<int> *dv_active_;
};

/** Failed particles retain their mass and transmit no stress or damping. */
template <typename...>
class JohnsonCookIntegration2ndHalfCK;

template <typename... Parameters>
class JohnsonCookIntegration2ndHalfCK<Inner<OneLevel, Parameters...>>
    : public StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>
{
    using Base = StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit JohnsonCookIntegration2ndHalfCK(DynamicsIdentifier &identifier);

    class InteractKernel : public Base::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t i, Real dt);

      protected:
        int *active_;
    };

  protected:
    DiscreteVariable<int> *dv_active_;
};

/** Same per-particle acoustic and acceleration bounds as the classic JC case. */
class JohnsonCookAcousticTimeStepCK : public LocalDynamicsReduce<ReduceMin<Real>>
{
  public:
    explicit JohnsonCookAcousticTimeStepCK(SPHBody &body, Real cfl = Real(0.3));

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Real reduce(size_t i, Real dt = 0);

      protected:
        Real h_min_, c0_, cfl_;
        Real *mass_;
        Vecd *vel_, *force_, *force_prior_;
    };

  protected:
    Real h_min_, c0_, cfl_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_force_, *dv_force_prior_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // JOHNSON_COOK_DYNAMICS_CK_H
