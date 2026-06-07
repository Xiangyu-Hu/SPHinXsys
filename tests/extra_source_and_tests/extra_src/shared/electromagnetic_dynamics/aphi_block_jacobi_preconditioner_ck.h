#ifndef APHI_BLOCK_JACOBI_PRECONDITIONER_CK_H
#define APHI_BLOCK_JACOBI_PRECONDITIONER_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "interaction_algorithms_ck.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

enum class AphiBlockJacobiPreconditionerKind
{
    /** Laplace/reaction/phi blocks inverted separately (ignores A--phi grad/div coupling). */
    Decoupled,
    /** Per-particle 8x8 coupled inverse including local grad(phi) and div(sigma A) diagonal coupling. */
    CoupledPointBlock8x8
};

inline const char *AphiBlockJacobiPreconditionerKindName(AphiBlockJacobiPreconditionerKind kind)
{
    switch (kind)
    {
    case AphiBlockJacobiPreconditionerKind::Decoupled:
        return "Decoupled";
    case AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8:
        return "CoupledPointBlock8x8";
    default:
        return "Unknown";
    }
}

struct AphiBlockJacobiDiagonalNames
{
    std::string laplace_a_diag = "JacobiLaplaceADiag";
    std::string laplace_phi_diag = "JacobiLaplacePhiDiag";
    /** sigma_i * sum_j g_ij for local grad(phi) coupling onto A. */
    std::string grad_phi_coupling = "JacobiGradPhiCoupling";
    /** omega * sum_j sigma_ij g_ij for local div(sigma A) coupling onto phi. */
    std::string div_a_coupling = "JacobiDivACoupling";
    /** Uncorrected pairwise local 3x3 block for -grad(div A) (real/imag use the same block). */
    std::string graddiv_a_block = "JacobiGradDivABlock";
};

/** Optional per-particle diagnostics for CoupledPointBlock8x8 apply (Stage 10A). */
struct AphiBlockJacobiDiagnosticNames
{
    std::string fallback_flag = "Jacobi8x8FallbackFlag";
    std::string min_pivot = "Jacobi8x8MinPivot";
};

class RegisterAphiBlockJacobiDiagnosticFieldsCK : public LocalDynamics
{
  public:
    explicit RegisterAphiBlockJacobiDiagnosticFieldsCK(SPHBody &sph_body,
                                                       const AphiBlockJacobiDiagnosticNames &diag_names =
                                                           AphiBlockJacobiDiagnosticNames{});
    virtual ~RegisterAphiBlockJacobiDiagnosticFieldsCK() = default;
};

template <typename... RelationTypes>
class AphiComputeBlockJacobiDiagonalCK;

/** Sum neighbor Laplace weights for block-Jacobi diagonal (ignores grad/div coupling). */
template <typename... Parameters>
class AphiComputeBlockJacobiDiagonalCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiComputeBlockJacobiDiagonalCK(Inner<Parameters...> &inner_relation, const AphiMaterialNames &material_names,
                                              Real omega, const AphiLhsAssemblyOptions &options,
                                              const AphiBlockJacobiDiagonalNames &diag_names = AphiBlockJacobiDiagonalNames{},
                                              Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiComputeBlockJacobiDiagonalCK(
        DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiComputeBlockJacobiDiagonalCK(parameters.identifier_, std::get<0>(parameters.others_),
                                           std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiComputeBlockJacobiDiagonalCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *sigma_;
        Real *nu_;
        Real *laplace_a_diag_;
        Real *laplace_phi_diag_;
        Vecd *grad_phi_coupling_;
        Vecd *div_a_coupling_;
        Matd *graddiv_a_block_;
        bool use_a_divergence_penalty_;
        Real omega_;
        AphiLhsTermFlags terms_;
        bool use_phi_gauge_penalty_;
        Real phi_gauge_penalty_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    AphiLhsAssemblyOptions options_;
    AphiBlockJacobiDiagonalNames diag_names_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    DiscreteVariable<Real> *dv_laplace_a_diag_;
    DiscreteVariable<Real> *dv_laplace_phi_diag_;
    DiscreteVariable<Vecd> *dv_grad_phi_coupling_;
    DiscreteVariable<Vecd> *dv_div_a_coupling_;
    DiscreteVariable<Matd> *dv_graddiv_a_block_;
};

/** Contact neighbor contribution to block-Jacobi diagonal; accumulates onto owner fields. */
template <typename... Parameters>
class AphiComputeBlockJacobiDiagonalCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiComputeBlockJacobiDiagonalCK(Contact<Parameters...> &contact_relation,
                                              const AphiMaterialNames &material_names, Real omega,
                                              const AphiLhsAssemblyOptions &options,
                                              const AphiBlockJacobiDiagonalNames &diag_names = AphiBlockJacobiDiagonalNames{},
                                              Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg>
    explicit AphiComputeBlockJacobiDiagonalCK(
        DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg, ThirdArg> parameters)
        : AphiComputeBlockJacobiDiagonalCK(parameters.identifier_, std::get<0>(parameters.others_),
                                           std::get<1>(parameters.others_), std::get<2>(parameters.others_)){};
    virtual ~AphiComputeBlockJacobiDiagonalCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
        Real *contact_Vol_;
        Real *sigma_;
        Real *nu_;
        Real *contact_sigma_;
        Real *contact_nu_;
        Real *laplace_a_diag_;
        Real *laplace_phi_diag_;
        Vecd *grad_phi_coupling_;
        Vecd *div_a_coupling_;
        Matd *graddiv_a_block_;
        bool use_a_divergence_penalty_;
        AphiContactADivergencePenaltyStencilMode contact_a_divergence_penalty_stencil_;
        Real omega_;
        AphiLhsTermFlags terms_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    AphiLhsAssemblyOptions options_;
    AphiBlockJacobiDiagonalNames diag_names_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    DiscreteVariable<Real> *dv_laplace_a_diag_;
    DiscreteVariable<Real> *dv_laplace_phi_diag_;
    DiscreteVariable<Vecd> *dv_grad_phi_coupling_;
    DiscreteVariable<Vecd> *dv_div_a_coupling_;
    DiscreteVariable<Matd> *dv_graddiv_a_block_;
    StdVec<DiscreteVariable<Real> *> dv_contact_sigma_;
    StdVec<DiscreteVariable<Real> *> dv_contact_nu_;
};

/** Apply block-Jacobi inverse: dst = M^{-1} src. A block uses local 2x2 real/imag inverse when reaction is on. */
class AphiApplyBlockJacobiInverseCK : public LocalDynamics
{
  public:
    explicit AphiApplyBlockJacobiInverseCK(
        SPHBody &sph_body, const AphiBlockNames &src_names, const AphiBlockNames &dst_names,
        const AphiMaterialNames &material_names, Real omega, const AphiLhsAssemblyOptions &options,
        AphiBlockJacobiPreconditionerKind preconditioner_kind = AphiBlockJacobiPreconditionerKind::Decoupled,
        const AphiBlockJacobiDiagonalNames &diag_names = AphiBlockJacobiDiagonalNames{},
        const AphiBlockJacobiDiagnosticNames *diagnostic_names = nullptr);
    virtual ~AphiApplyBlockJacobiInverseCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real omega_;
        AphiLhsTermFlags terms_;
        bool use_phi_gauge_penalty_;
        Real phi_gauge_penalty_;
        bool use_a_divergence_penalty_;
        Real a_divergence_penalty_;
        AphiBlockJacobiPreconditionerKind preconditioner_kind_;
        Real *sigma_;
        Real *laplace_a_diag_;
        Real *laplace_phi_diag_;
        Vecd *grad_phi_coupling_;
        Vecd *div_a_coupling_;
        Matd *graddiv_a_block_;
        Vecd *src_a_real_;
        Vecd *src_a_imag_;
        Real *src_phi_real_;
        Real *src_phi_imag_;
        Vecd *dst_a_real_;
        Vecd *dst_a_imag_;
        Real *dst_phi_real_;
        Real *dst_phi_imag_;
        Real *fallback_flag_;
        Real *min_pivot_;
    };

  protected:
    Real omega_;
    AphiLhsAssemblyOptions options_;
    AphiBlockJacobiPreconditionerKind preconditioner_kind_;
    bool record_diagnostics_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_laplace_a_diag_;
    DiscreteVariable<Real> *dv_laplace_phi_diag_;
    DiscreteVariable<Vecd> *dv_grad_phi_coupling_;
    DiscreteVariable<Vecd> *dv_div_a_coupling_;
    DiscreteVariable<Matd> *dv_graddiv_a_block_;
    DiscreteVariable<Vecd> *dv_src_a_real_;
    DiscreteVariable<Vecd> *dv_src_a_imag_;
    DiscreteVariable<Real> *dv_src_phi_real_;
    DiscreteVariable<Real> *dv_src_phi_imag_;
    DiscreteVariable<Vecd> *dv_dst_a_real_;
    DiscreteVariable<Vecd> *dv_dst_a_imag_;
    DiscreteVariable<Real> *dv_dst_phi_real_;
    DiscreteVariable<Real> *dv_dst_phi_imag_;
    DiscreteVariable<Real> *dv_fallback_flag_;
    DiscreteVariable<Real> *dv_min_pivot_;
};

/** Inner + Contact block-Jacobi diagonal assembly (Inner assigns, Contact accumulates). */
template <class ExecutionPolicy>
class AphiComputeBlockJacobiContactDynamicsBundle
{
  public:
    AphiComputeBlockJacobiContactDynamicsBundle(SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
                                                const AphiMaterialNames &material_names, Real omega,
                                                const AphiLhsAssemblyOptions &options,
                                                const AphiBlockJacobiDiagonalNames &diag_names = AphiBlockJacobiDiagonalNames{},
                                                Real pair_weight_regularization = Real(0.01));

    void exec();

  protected:
    InteractionDynamicsCK<ExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_inner_;
    InteractionDynamicsCK<ExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Contact<>>> compute_contact_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BLOCK_JACOBI_PRECONDITIONER_CK_H
