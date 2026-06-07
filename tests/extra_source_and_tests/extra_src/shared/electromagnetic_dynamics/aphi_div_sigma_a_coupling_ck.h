#ifndef APHI_DIV_SIGMA_A_COUPLING_CK_H
#define APHI_DIV_SIGMA_A_COUPLING_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... RelationTypes>
class AphiDivSigmaACouplingCK;

/** OPERATOR default: omega * sigma_ij * g_ij dot (A_i - A_j) on LhsPhi. */
template <typename... Parameters>
class AphiDivSigmaACouplingCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiDivSigmaACouplingCK(Inner<Parameters...> &inner_relation, Real omega,
                                     const AphiVariableNames &variable_names);
    template <typename FirstArg, typename SecondArg>
    explicit AphiDivSigmaACouplingCK(DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiDivSigmaACouplingCK(parameters.identifier_, std::get<0>(parameters.others_),
                                  std::get<1>(parameters.others_)){};
    virtual ~AphiDivSigmaACouplingCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *sigma_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;
        Real omega_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;
    DiscreteVariable<Real> *dv_sigma_;
};

/** Contact pairwise div(sigma A): omega * sigma_ij * g_ij dot (A_i - A_j) on LhsPhi. */
template <typename... Parameters>
class AphiDivSigmaACouplingCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiDivSigmaACouplingCK(Contact<Parameters...> &contact_relation, Real omega,
                                     const AphiVariableNames &variable_names);
    template <typename FirstArg, typename SecondArg>
    explicit AphiDivSigmaACouplingCK(DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg> parameters)
        : AphiDivSigmaACouplingCK(parameters.identifier_, std::get<0>(parameters.others_),
                                  std::get<1>(parameters.others_)){};
    virtual ~AphiDivSigmaACouplingCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        Real *sigma_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *contact_sigma_;
        Vecd *contact_a_real_;
        Vecd *contact_a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;
        Real omega_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;
    DiscreteVariable<Real> *dv_sigma_;
    StdVec<DiscreteVariable<Real> *> dv_contact_sigma_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_a_real_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_a_imag_;
};

/**
 * DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — constant sigma only.
 * LhsPhiReal += omega*sigma*DivAImag, LhsPhiImag -= omega*sigma*DivAReal.
 */
class AphiDivSigmaAConstSigmaDiagnosticCouplingCK : public LocalDynamics
{
  public:
    explicit AphiDivSigmaAConstSigmaDiagnosticCouplingCK(SPHBody &sph_body, Real omega,
                                                         const AphiVariableNames &variable_names,
                                                         const AphiDiagnosticNames &diagnostic_names);
    virtual ~AphiDivSigmaAConstSigmaDiagnosticCouplingCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *sigma_;
        Real *div_a_real_;
        Real *div_a_imag_;
        Real *lhs_phi_real_;
        Real *lhs_phi_imag_;
        Real omega_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_div_a_real_;
    DiscreteVariable<Real> *dv_div_a_imag_;
    DiscreteVariable<Real> *dv_lhs_phi_real_;
    DiscreteVariable<Real> *dv_lhs_phi_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIV_SIGMA_A_COUPLING_CK_H
