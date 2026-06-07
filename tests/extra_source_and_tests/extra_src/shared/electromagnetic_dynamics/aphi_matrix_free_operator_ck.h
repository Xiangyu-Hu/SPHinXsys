#ifndef APHI_MATRIX_FREE_OPERATOR_CK_H
#define APHI_MATRIX_FREE_OPERATOR_CK_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_contact_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/diagnostics/aphi_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_block_zero_ck.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "interaction_algorithms_ck.h"
#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... RelationTypes>
class AphiApplyCK;

/**
 * Fused matrix-free K(input) on Inner<> without global intermediate Lap/Grad/Div fields.
 * Current implementation uses one InteractKernel and multiple optional neighbor loops.
 * A later optimization may merge them into one loop.
 * Supports only PairwiseUncorrected grad/div modes; caller must zero output before exec
 * (or use AphiApplyDynamicsBundle). Diagnostic B-corrected modes are only in debug assemble.
 */
template <typename... Parameters>
class AphiApplyCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit AphiApplyCK(Inner<Parameters...> &inner_relation, const AphiBlockNames &input_block,
                         const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
                         const AphiLhsAssemblyOptions &options, Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg, typename FourthArg, typename FifthArg,
              typename SixthArg>
    explicit AphiApplyCK(
        DynamicsArgs<Inner<Parameters...>, FirstArg, SecondArg, ThirdArg, FourthArg, FifthArg, SixthArg> parameters)
        : AphiApplyCK(parameters.identifier_, std::get<0>(parameters.others_), std::get<1>(parameters.others_),
                      std::get<2>(parameters.others_), std::get<3>(parameters.others_),
                      std::get<4>(parameters.others_), std::get<5>(parameters.others_)){};
    virtual ~AphiApplyCK() = default;

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
        Vecd *in_a_real_;
        Vecd *in_a_imag_;
        Real *in_phi_real_;
        Real *in_phi_imag_;
        Vecd *out_a_real_;
        Vecd *out_a_imag_;
        Real *out_phi_real_;
        Real *out_phi_imag_;
        AphiLhsTermFlags terms_;
        Real omega_;
        bool use_phi_gauge_penalty_;
        Real phi_gauge_penalty_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    AphiLhsAssemblyOptions options_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Vecd> *dv_in_a_real_;
    DiscreteVariable<Vecd> *dv_in_a_imag_;
    DiscreteVariable<Real> *dv_in_phi_real_;
    DiscreteVariable<Real> *dv_in_phi_imag_;
    DiscreteVariable<Vecd> *dv_out_a_real_;
    DiscreteVariable<Vecd> *dv_out_a_imag_;
    DiscreteVariable<Real> *dv_out_phi_real_;
    DiscreteVariable<Real> *dv_out_phi_imag_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
};

/** Contact contribution to fused K(X); accumulates onto output block (Inner assigns first). */
template <typename... Parameters>
class AphiApplyCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiApplyCK(Contact<Parameters...> &contact_relation, const AphiBlockNames &input_block,
                         const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
                         const AphiLhsAssemblyOptions &options, Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg, typename FourthArg, typename FifthArg,
              typename SixthArg>
    explicit AphiApplyCK(
        DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg, ThirdArg, FourthArg, FifthArg, SixthArg> parameters)
        : AphiApplyCK(parameters.identifier_, std::get<0>(parameters.others_), std::get<1>(parameters.others_),
                      std::get<2>(parameters.others_), std::get<3>(parameters.others_),
                      std::get<4>(parameters.others_), std::get<5>(parameters.others_)){};
    virtual ~AphiApplyCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        Real *sigma_;
        Real *nu_;
        Real *contact_sigma_;
        Real *contact_nu_;
        Real *in_phi_real_;
        Real *in_phi_imag_;
        Vecd *in_a_real_;
        Vecd *in_a_imag_;
        Real *contact_in_phi_real_;
        Real *contact_in_phi_imag_;
        Vecd *contact_in_a_real_;
        Vecd *contact_in_a_imag_;
        Vecd *out_a_real_;
        Vecd *out_a_imag_;
        Real *out_phi_real_;
        Real *out_phi_imag_;
        AphiLhsTermFlags terms_;
        Real omega_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    AphiLhsAssemblyOptions options_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Vecd> *dv_in_a_real_;
    DiscreteVariable<Vecd> *dv_in_a_imag_;
    DiscreteVariable<Real> *dv_in_phi_real_;
    DiscreteVariable<Real> *dv_in_phi_imag_;
    DiscreteVariable<Vecd> *dv_out_a_real_;
    DiscreteVariable<Vecd> *dv_out_a_imag_;
    DiscreteVariable<Real> *dv_out_phi_real_;
    DiscreteVariable<Real> *dv_out_phi_imag_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    StdVec<DiscreteVariable<Real> *> dv_contact_sigma_;
    StdVec<DiscreteVariable<Real> *> dv_contact_nu_;
    StdVec<DiscreteVariable<Real> *> dv_contact_in_phi_real_;
    StdVec<DiscreteVariable<Real> *> dv_contact_in_phi_imag_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_in_a_real_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_in_a_imag_;
};

/** Zero output block then fused apply: output = K(input). Optionally appends grad(div A) penalty. */
template <class ExecutionPolicy>
class AphiApplyDynamicsBundle
{
  public:
    AphiApplyDynamicsBundle(SPHBody &sph_body, Inner<> &inner_relation, const AphiBlockNames &input_block,
                            const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
                            const AphiLhsAssemblyOptions &options, Real pair_weight_regularization = Real(0.01));

    void exec();

  protected:
    AphiLhsAssemblyOptions options_;
    AphiBlockNames input_block_;
    AphiBlockNames output_block_;
    SPHBody &sph_body_;
    Inner<> &inner_relation_;
    StateDynamics<ExecutionPolicy, AphiZeroBlockCK> zero_output_;
    InteractionDynamicsCK<ExecutionPolicy, AphiApplyCK<Inner<>>> apply_;
    AphiPairwiseADivergencePenaltyPipelineBundle<ExecutionPolicy> a_divergence_penalty_pairwise_;
    AphiADivergencePenaltyPipelineBundle<ExecutionPolicy> a_divergence_penalty_b_corrected_;
};

/** Zero output then Inner + Contact fused apply for multi-body Contact<> relations. */
template <class ExecutionPolicy>
class AphiApplyContactDynamicsBundle
{
  public:
    AphiApplyContactDynamicsBundle(SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
                                   const AphiBlockNames &input_block, const AphiBlockNames &output_block,
                                   const AphiMaterialNames &material_names, Real omega,
                                   const AphiLhsAssemblyOptions &options, Real pair_weight_regularization = Real(0.01));

    void exec();

  protected:
    AphiLhsAssemblyOptions options_;
    StateDynamics<ExecutionPolicy, AphiZeroBlockCK> zero_output_;
    InteractionDynamicsCK<ExecutionPolicy, AphiApplyCK<Inner<>>> apply_inner_;
    InteractionDynamicsCK<ExecutionPolicy, AphiApplyCK<Contact<>>> apply_contact_;
    AphiContactPairwiseADivergencePenaltyPipelineBundle<ExecutionPolicy> a_divergence_penalty_inner_contact_;
    AphiPairwiseADivergencePenaltyPipelineBundle<ExecutionPolicy> a_divergence_penalty_inner_only_;
};

/**
 * Contact block-diagonal K_ii(z_i): Inner full stencil + Contact diagonal-only cross-body terms.
 * Used for single-body GMRES Arnoldi; does not read neighbor body Krylov input blocks.
 */
template <typename... RelationTypes>
class AphiApplyContactBlockDiagonalCK;

template <typename... Parameters>
class AphiApplyContactBlockDiagonalCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit AphiApplyContactBlockDiagonalCK(Contact<Parameters...> &contact_relation, const AphiBlockNames &input_block,
                                             const AphiBlockNames &output_block, const AphiMaterialNames &material_names,
                                             Real omega, const AphiLhsAssemblyOptions &options,
                                             Real pair_weight_regularization = Real(0.01));
    template <typename FirstArg, typename SecondArg, typename ThirdArg, typename FourthArg, typename FifthArg,
              typename SixthArg>
    explicit AphiApplyContactBlockDiagonalCK(
        DynamicsArgs<Contact<Parameters...>, FirstArg, SecondArg, ThirdArg, FourthArg, FifthArg, SixthArg> parameters)
        : AphiApplyContactBlockDiagonalCK(parameters.identifier_, std::get<0>(parameters.others_),
                                        std::get<1>(parameters.others_), std::get<2>(parameters.others_),
                                        std::get<3>(parameters.others_), std::get<4>(parameters.others_),
                                        std::get<5>(parameters.others_)){};
    virtual ~AphiApplyContactBlockDiagonalCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        Real *sigma_;
        Real *nu_;
        Real *contact_sigma_;
        Real *contact_nu_;
        Real *in_phi_real_;
        Real *in_phi_imag_;
        Vecd *in_a_real_;
        Vecd *in_a_imag_;
        Vecd *out_a_real_;
        Vecd *out_a_imag_;
        Real *out_phi_real_;
        Real *out_phi_imag_;
        AphiLhsTermFlags terms_;
        Real omega_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    AphiLhsAssemblyOptions options_;
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Vecd> *dv_in_a_real_;
    DiscreteVariable<Vecd> *dv_in_a_imag_;
    DiscreteVariable<Real> *dv_in_phi_real_;
    DiscreteVariable<Real> *dv_in_phi_imag_;
    DiscreteVariable<Vecd> *dv_out_a_real_;
    DiscreteVariable<Vecd> *dv_out_a_imag_;
    DiscreteVariable<Real> *dv_out_phi_real_;
    DiscreteVariable<Real> *dv_out_phi_imag_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_nu_;
    StdVec<DiscreteVariable<Real> *> dv_contact_sigma_;
    StdVec<DiscreteVariable<Real> *> dv_contact_nu_;
};

/** Inner full + Contact block-diagonal apply for single-body Krylov matvec. */
template <class ExecutionPolicy>
class AphiApplyContactBlockDiagonalDynamicsBundle
{
  public:
    AphiApplyContactBlockDiagonalDynamicsBundle(SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
                                                const AphiBlockNames &input_block, const AphiBlockNames &output_block,
                                                const AphiMaterialNames &material_names, Real omega,
                                                const AphiLhsAssemblyOptions &options,
                                                Real pair_weight_regularization = Real(0.01));

    void exec();

  protected:
    StateDynamics<ExecutionPolicy, AphiZeroBlockCK> zero_output_;
    InteractionDynamicsCK<ExecutionPolicy, AphiApplyCK<Inner<>>> apply_inner_;
    InteractionDynamicsCK<ExecutionPolicy, AphiApplyContactBlockDiagonalCK<Contact<>>> apply_contact_block_diagonal_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_MATRIX_FREE_OPERATOR_CK_H
