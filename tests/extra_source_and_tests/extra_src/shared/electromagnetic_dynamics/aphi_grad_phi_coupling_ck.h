#ifndef APHI_GRAD_PHI_COUPLING_CK_H
#define APHI_GRAD_PHI_COUPLING_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... RelationTypes>
class AphiGradPhiCouplingCK;

/** OPERATOR default: uncorrected pairwise sigma_i * g_ij * (phi_i - phi_j) on LhsA. */
template <typename... Parameters>
class AphiGradPhiCouplingCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiGradPhiCouplingCK(Inner<Parameters...> &inner_relation, const AphiVariableNames &variable_names);
    template <typename FirstArg>
    explicit AphiGradPhiCouplingCK(DynamicsArgs<Inner<Parameters...>, FirstArg> parameters)
        : AphiGradPhiCouplingCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AphiGradPhiCouplingCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *sigma_;
        Real *phi_real_;
        Real *phi_imag_;
        Vecd *lhs_a_real_;
        Vecd *lhs_a_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_lhs_a_real_;
    DiscreteVariable<Vecd> *dv_lhs_a_imag_;
    DiscreteVariable<Real> *dv_sigma_;
};

/**
 * DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — constant sigma.
 * Adds sigma_i * LinearGradient(phi) to LhsA. Requires phi*Gradient fields from LinearGradient.
 */
class AphiGradPhiDiagnosticCouplingCK : public LocalDynamics
{
  public:
    explicit AphiGradPhiDiagnosticCouplingCK(SPHBody &sph_body, const AphiVariableNames &variable_names);
    virtual ~AphiGradPhiDiagnosticCouplingCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *sigma_;
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
        Vecd *lhs_a_real_;
        Vecd *lhs_a_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    DiscreteVariable<Vecd> *dv_lhs_a_real_;
    DiscreteVariable<Vecd> *dv_lhs_a_imag_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GRAD_PHI_COUPLING_CK_H
