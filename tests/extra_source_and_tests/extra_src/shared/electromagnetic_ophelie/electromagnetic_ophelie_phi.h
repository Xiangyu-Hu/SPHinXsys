#ifndef ELECTROMAGNETIC_OPHELIE_PHI_H
#define ELECTROMAGNETIC_OPHELIE_PHI_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "electromagnetic_ophelie_phi_gradient.h"
#include "electromagnetic_ophelie_phi_rhs_diagnostics.h"
#include "interaction_algorithms_ck.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

class ZeroOphelieScalarFieldCK : public LocalDynamics
{
  public:
    ZeroOphelieScalarFieldCK(SPHBody &sph_body, const std::string &field_name)
        : LocalDynamics(sph_body), dv_field_(particles_->template getVariableByName<Real>(field_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : field_(encloser.dv_field_->DelegatedData(ex_policy))
        {
        }
        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            field_[index_i] = 0.0;
        }

      protected:
        Real *field_;
    };

  protected:
    DiscreteVariable<Real> *dv_field_;
};

/** output = sigma * input (Vecd fields). */
class OphelieScaleVecdFieldBySigmaCK : public LocalDynamics
{
  public:
    OphelieScaleVecdFieldBySigmaCK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
                                   const std::string &input_field, const std::string &output_field)
        : LocalDynamics(sph_body), dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
          dv_input_(particles_->template getVariableByName<Vecd>(input_field)),
          dv_output_(particles_->template getVariableByName<Vecd>(output_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              input_(encloser.dv_input_->DelegatedData(ex_policy)),
              output_(encloser.dv_output_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            output_[index_i] = sigma_[index_i] * input_[index_i];
        }

      protected:
        Real *sigma_;
        Vecd *input_;
        Vecd *output_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_input_;
    DiscreteVariable<Vecd> *dv_output_;
};

class OphelieScaleScalarFieldCK : public LocalDynamics
{
  public:
    OphelieScaleScalarFieldCK(SPHBody &sph_body, const std::string &field_name, Real scale_factor)
        : LocalDynamics(sph_body), scale_factor_(scale_factor),
          dv_field_(particles_->template getVariableByName<Real>(field_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : scale_factor_(encloser.scale_factor_), field_(encloser.dv_field_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            field_[index_i] *= scale_factor_;
        }

      protected:
        Real scale_factor_;
        Real *field_;
    };

  protected:
    Real scale_factor_;
    DiscreteVariable<Real> *dv_field_;
};

template <typename... RelationTypes>
class ComputeOpheliePhiImagRhsFromASrcCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOpheliePhiImagRhsFromASrcCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOpheliePhiImagRhsFromASrcCK(RelationType<Parameters...> &relation, Real omega,
                                       const OphelieGlassFieldNames &names)
        : BaseInteraction(relation), omega_(omega),
          dv_a_src_real_(this->particles_->template getVariableByName<Vecd>(names.a_src_real)),
          dv_phi_rhs_imag_(this->particles_->template getVariableByName<Real>(names.phi_rhs_imag)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
              phi_rhs_imag_(encloser.dv_phi_rhs_imag_->DelegatedData(ex_policy)), omega_(encloser.omega_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real sigma_i = sigma_[index_i];
            const Vecd a_re_i = a_src_real_[index_i];
            Real rhs_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
                const Real sigma_ij = harmonicMean(sigma_i, sigma_[index_j]);
                // Pairwise flux; for constant sigma matches +omega*sigma*div(A) (opposite sign to DivSigmaA RHS).
                rhs_i -= omega_ * sigma_ij * g_ij.dot(a_re_i - a_src_real_[index_j]);
            }
            phi_rhs_imag_[index_i] += rhs_i;
        }

      protected:
        Real *Vol_;
        Real *sigma_;
        Vecd *a_src_real_;
        Real *phi_rhs_imag_;
        Real omega_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Real> *dv_phi_rhs_imag_;
    DiscreteVariable<Real> *dv_sigma_;
};

template <typename... Parameters>
class ComputeOpheliePhiImagRhsFromASrcCK<Inner<Parameters...>>
    : public ComputeOpheliePhiImagRhsFromASrcCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOpheliePhiImagRhsFromASrcCK(Inner<Parameters...> &inner_relation, Real omega,
                                       const OphelieGlassFieldNames &names)
        : ComputeOpheliePhiImagRhsFromASrcCK<Base, Inner<Parameters...>>(inner_relation, omega, names)
    {
    }
};

template <typename... RelationTypes>
class ComputeOphelieScalarPhiGradientCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieScalarPhiGradientCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieScalarPhiGradientCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names)
        : BaseInteraction(relation),
          dv_phi_imag_(this->particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_grad_phi_imag_(this->particles_->template getVariableByName<Vecd>(names.grad_phi_imag))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_imag_[index_i];
            Vecd grad_phi = Vecd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
                grad_phi += g_ij * (phi_i - phi_imag_[index_j]);
            }
            grad_phi_imag_[index_i] = grad_phi;
        }

      protected:
        Real *Vol_;
        Real *phi_imag_;
        Vecd *grad_phi_imag_;
    };

  protected:
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
};

template <typename... Parameters>
class ComputeOphelieScalarPhiGradientCK<Inner<Parameters...>>
    : public ComputeOphelieScalarPhiGradientCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieScalarPhiGradientCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names)
        : ComputeOphelieScalarPhiGradientCK<Base, Inner<Parameters...>>(inner_relation, names)
    {
    }
};

/** Add gauge penalty term onto Laplace lhs: lhs_phi += penalty * phi. */
class OpheliePhiImagGaugePenaltyToLhsCK : public LocalDynamics
{
  public:
    OpheliePhiImagGaugePenaltyToLhsCK(SolidBody &sph_body, const OphelieGlassFieldNames &names, Real phi_gauge_penalty)
        : LocalDynamics(sph_body), phi_gauge_penalty_(phi_gauge_penalty),
          dv_phi_imag_(particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_phi_lhs_imag_(particles_->template getVariableByName<Real>(names.phi_lhs_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : phi_gauge_penalty_(encloser.phi_gauge_penalty_),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              phi_lhs_imag_(encloser.dv_phi_lhs_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            phi_lhs_imag_[index_i] += phi_gauge_penalty_ * phi_imag_[index_i];
        }

      protected:
        Real phi_gauge_penalty_;
        Real *phi_imag_;
        Real *phi_lhs_imag_;
    };

  protected:
    Real phi_gauge_penalty_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Real> *dv_phi_lhs_imag_;
};

class OphelieJacobiPhiImagUpdateCK : public LocalDynamics
{
  public:
    OphelieJacobiPhiImagUpdateCK(SolidBody &sph_body, const OphelieGlassFieldNames &names, Real phi_gauge_penalty,
                                 Real jacobi_relaxation)
        : LocalDynamics(sph_body), phi_gauge_penalty_(phi_gauge_penalty), jacobi_relaxation_(jacobi_relaxation),
          dv_phi_imag_(particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_phi_rhs_imag_(particles_->template getVariableByName<Real>(names.phi_rhs_imag)),
          dv_phi_lhs_imag_(particles_->template getVariableByName<Real>(names.phi_lhs_imag)),
          dv_phi_diag_(particles_->template getVariableByName<Real>(names.phi_laplace_diag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : phi_gauge_penalty_(encloser.phi_gauge_penalty_), jacobi_relaxation_(encloser.jacobi_relaxation_),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              phi_rhs_imag_(encloser.dv_phi_rhs_imag_->DelegatedData(ex_policy)),
              phi_lhs_imag_(encloser.dv_phi_lhs_imag_->DelegatedData(ex_policy)),
              phi_diag_(encloser.dv_phi_diag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real diag_i = phi_diag_[index_i];
            const Real denom = diag_i + phi_gauge_penalty_ + TinyReal;
            const Real phi_i = phi_imag_[index_i];
            const Real numerator = phi_rhs_imag_[index_i] + diag_i * phi_i - phi_lhs_imag_[index_i];
            const Real phi_jacobi = numerator / denom;
            phi_imag_[index_i] = (Real(1.0) - jacobi_relaxation_) * phi_i + jacobi_relaxation_ * phi_jacobi;
        }

      protected:
        Real phi_gauge_penalty_;
        Real jacobi_relaxation_;
        Real *phi_imag_;
        Real *phi_rhs_imag_;
        Real *phi_lhs_imag_;
        Real *phi_diag_;
    };

  protected:
    Real phi_gauge_penalty_;
    Real jacobi_relaxation_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Real> *dv_phi_rhs_imag_;
    DiscreteVariable<Real> *dv_phi_lhs_imag_;
    DiscreteVariable<Real> *dv_phi_diag_;
};

class ComputeOphelieEJQWithPhiCK : public LocalDynamics
{
  public:
    ComputeOphelieEJQWithPhiCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, const OphelieParameters &params);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real omega_;
        Real *sigma_;
        Vecd *a_src_real_;
        Vecd *grad_phi_imag_;
        Vecd *e_real_;
        Vecd *e_imag_;
        Vecd *j_real_;
        Vecd *j_imag_;
        Real *joule_heat_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    DiscreteVariable<Vecd> *dv_e_real_;
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Real> *dv_joule_heat_;
};

template <class ExecutionPolicy>
inline void setupOpheliePhiImagRhsFromASrcDivConsistent(SolidBody &glass_body, Inner<> &inner,
                                                        const OphelieGlassFieldNames &names,
                                                        const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_rhs(glass_body, names.phi_rhs_imag);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdFieldBySigmaCK> scale_sigma_a(
        glass_body, names, names.a_src_real, names.j_imag);
    StateDynamics<ExecutionPolicy, OphelieScaleScalarFieldCK> scale_rhs(
        glass_body, names.phi_rhs_imag, -params.omega());

    zero_rhs.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (opheliePhiUseCorrectedDivergence(params))
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params);
    }
    scale_sigma_a.exec();
    execOphelieVecdDivergence<ExecutionPolicy>(inner, names, params, names.j_imag, names.phi_rhs_imag);
    scale_rhs.exec();
}

template <class ExecutionPolicy>
inline void setupOpheliePhiImagRhsFromASrcLegacyFlux(SolidBody &glass_body, Inner<> &inner,
                                                     const OphelieGlassFieldNames &names,
                                                     const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_rhs(glass_body, names.phi_rhs_imag);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);

    zero_rhs.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        const OphelieEdgeFluxComponent imag_component = makeOphelieEdgeFluxImagComponent(names, params);
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Inner<>>> compute_rhs(
            inner, names, imag_component, params.omega(), params.pair_weight_regularization_);
        compute_rhs.exec();
        return;
    }
    InteractionDynamicsCK<ExecutionPolicy, ComputeOpheliePhiImagRhsFromASrcCK<Inner<>>> compute_rhs(
        inner, params.omega(), names);
    compute_rhs.exec();
}

template <class ExecutionPolicy>
inline void setupOpheliePhiImagRhsFromASrc(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                           const OphelieParameters &params)
{
    if (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::DivSigmaA)
    {
        setupOpheliePhiImagRhsFromASrcDivConsistent<ExecutionPolicy>(glass_body, inner, names, params);
        return;
    }
    setupOpheliePhiImagRhsFromASrcLegacyFlux<ExecutionPolicy>(glass_body, inner, names, params);
}

template <class ExecutionPolicy>
inline void setupOpheliePhiImagRhsProblem(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                          const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, names.phi_imag);
    zero_phi.exec();
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
}

template <class ExecutionPolicy>
inline void applyOpheliePhiImagDivSigmaGradLhsOperator(SolidBody &glass_body, Inner<> &inner,
                                                       const OphelieGlassFieldNames &names,
                                                       const OphelieParameters &params)
{
    (void)params;
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_lhs(glass_body, names.phi_lhs_imag);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdFieldBySigmaCK> scale_sigma_grad(
        glass_body, names, names.grad_phi_imag, names.j_imag);
    StateDynamics<ExecutionPolicy, OpheliePhiImagGaugePenaltyToLhsCK> apply_gauge_penalty(
        glass_body, names, params.phi_gauge_penalty_);

    update_cell_linked_list.exec();
    update_inner_relation.exec();
    zero_lhs.exec();
    if (opheliePhiUseCorrectedGradient(params))
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params);
        execOphelieScalarPhiGradientCorrectedOnly<ExecutionPolicy>(inner, names);
    }
    else
    {
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(inner,
                                                                                                             names);
        compute_grad_phi.exec();
    }
    applyOpheliePhiBoundaryGradNeumannProjectionDynamics<ExecutionPolicy>(glass_body, names, params, true);
    scale_sigma_grad.exec();
    execOphelieVecdDivergence<ExecutionPolicy>(inner, names, params, names.j_imag, names.phi_lhs_imag);
    apply_gauge_penalty.exec();
}

template <class ExecutionPolicy>
inline void applyOpheliePhiImagLegacyPairwiseLhsOperator(SolidBody &glass_body, Inner<> &inner,
                                                         const OphelieGlassFieldNames &names,
                                                         const OphelieParameters &params)
{
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>> apply_laplace(
        inner, names.phi_imag, names.sigma, names.phi_lhs_imag, params.pair_weight_regularization_);
    StateDynamics<ExecutionPolicy, OpheliePhiImagGaugePenaltyToLhsCK> apply_gauge_penalty(
        glass_body, names, params.phi_gauge_penalty_);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_lhs(glass_body, names.phi_lhs_imag);

    zero_lhs.exec();
    apply_laplace.exec();
    apply_gauge_penalty.exec();
}

template <class ExecutionPolicy>
inline void computeOpheliePhiGMRESPreconditionerDiagonal(SolidBody &glass_body, Inner<> &inner,
                                                          const OphelieGlassFieldNames &names,
                                                          const OphelieParameters &params)
{
    if (params.phi_lhs_operator_kind_ != OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        computeOpheliePhiOperatorDiagonal<ExecutionPolicy>(glass_body, inner, names, params);
        return;
    }

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> div_diag(n, Real(0));
    StdVec<Real> legacy_diag(n, Real(0));

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_div(glass_body, names.phi_laplace_diag);
    InteractionDynamicsCK<ExecutionPolicy, OphelieDivSigmaGradDiagonalCK<Inner<>>> compute_div_diag(
        inner, names.sigma, names.phi_laplace_diag);
    zero_div.exec();
    compute_div_diag.exec();
    hostReadScalarField(particles, names.phi_laplace_diag, div_diag.data(), n);

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_legacy(glass_body, names.phi_laplace_diag);
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_legacy_diag(
        inner, names.sigma, names.phi_laplace_diag, params.pair_weight_regularization_);
    zero_legacy.exec();
    compute_legacy_diag.exec();
    hostReadScalarField(particles, names.phi_laplace_diag, legacy_diag.data(), n);

    for (size_t i = 0; i < n; ++i)
    {
        div_diag[i] = std::max(div_diag[i], legacy_diag[i]);
    }
    hostAssignScalarField(particles, names.phi_laplace_diag, div_diag.data(), n);
}

template <class ExecutionPolicy>
inline void computeOpheliePhiOperatorDiagonal(SolidBody &glass_body, Inner<> &inner,
                                              const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_diag(glass_body, names.phi_laplace_diag);
    zero_diag.exec();
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        InteractionDynamicsCK<ExecutionPolicy, OphelieDivSigmaGradDiagonalCK<Inner<>>> compute_diag(
            inner, names.sigma, names.phi_laplace_diag);
        compute_diag.exec();
        return;
    }
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_diag(
        inner, names.sigma, names.phi_laplace_diag, params.pair_weight_regularization_);
    compute_diag.exec();
}

template <class ExecutionPolicy>
inline void applyOpheliePhiImagLhsOperator(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                           const OphelieParameters &params)
{
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        applyOpheliePhiImagDivSigmaGradLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        return;
    }
    applyOpheliePhiImagLegacyPairwiseLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
}

template <class ExecutionPolicy>
inline Real evaluateOpheliePhiImagRelativeResidual(SolidBody &glass_body, Inner<> &inner,
                                                   const OphelieGlassFieldNames &names,
                                                   const OphelieParameters &params)
{
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    BaseParticles &particles = glass_body.getBaseParticles();
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    const size_t n = particles.TotalRealParticles();
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    Real rhs_max = 0.0;
    Real residual_max = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real residual = std::abs(rhs[i] - lhs[i]);
        rhs_max = std::max(rhs_max, std::abs(rhs[i]));
        residual_max = std::max(residual_max, residual);
    }
    return residual_max / (rhs_max + TinyReal);
}

template <class ExecutionPolicy>
inline Real evaluateOpheliePhiImagJacobiRelativeResidual(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                         Real phi_gauge_penalty)
{
    syncVariableToHost<Real>(particles, names.phi_imag);
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    const size_t n = particles.TotalRealParticles();
    const Real *phi = particles.getVariableDataByName<Real>(names.phi_imag);
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    Real rhs_max = 0.0;
    Real residual_max = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real residual = std::abs(rhs[i] - lhs[i] - phi_gauge_penalty * phi[i]);
        rhs_max = std::max(rhs_max, std::abs(rhs[i]));
        residual_max = std::max(residual_max, residual);
    }
    return residual_max / (rhs_max + TinyReal);
}

template <class ExecutionPolicy>
inline Real solvePhiImagJacobiWithCurrentRhs(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                             const OphelieParameters &params)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, names.phi_imag);
    zero_phi.exec();
    logOpheliePhiRhsFingerprint("jacobi_entry", computeOpheliePhiRhsFingerprint(particles, names.phi_rhs_imag, n));

    computeOpheliePhiOperatorDiagonal<ExecutionPolicy>(glass_body, inner, names, params);
    StateDynamics<ExecutionPolicy, OphelieJacobiPhiImagUpdateCK> jacobi_update(
        glass_body, names, params.phi_gauge_penalty_, params.phi_jacobi_relaxation_);

    Real max_residual = -1.0;
    for (size_t iter = 0; iter < params.phi_jacobi_max_iterations_; ++iter)
    {
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        jacobi_update.exec();
        if ((iter + 1) % 10 == 0 || iter + 1 == params.phi_jacobi_max_iterations_)
        {
            applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
            max_residual = evaluateOpheliePhiImagJacobiRelativeResidual<ExecutionPolicy>(
                glass_body.getBaseParticles(), names, params.phi_gauge_penalty_);
            if (max_residual < params.phi_jacobi_tolerance_)
            {
                break;
            }
            syncVariableToDevice<Real>(glass_body.getBaseParticles(), names.phi_imag);
        }
    }
    if (max_residual < 0.0)
    {
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        max_residual = evaluateOpheliePhiImagJacobiRelativeResidual<ExecutionPolicy>(
            glass_body.getBaseParticles(), names, params.phi_gauge_penalty_);
    }
    return max_residual;
}

template <class ExecutionPolicy>
inline Real solvePhiImagJacobi(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                               const OphelieParameters &params)
{
    setupOpheliePhiImagRhsProblem<ExecutionPolicy>(glass_body, inner, names, params);
    return solvePhiImagJacobiWithCurrentRhs<ExecutionPolicy>(glass_body, inner, names, params);
}

/** Jacobi smoothing with fixed RHS (call after GMRES/PCG). */
template <class ExecutionPolicy>
inline Real refinePhiImagJacobiPostSolve(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                         const OphelieParameters &params)
{
    if (params.phi_post_jacobi_refinement_iterations_ == 0)
    {
        return evaluateOpheliePhiImagRelativeResidual<ExecutionPolicy>(glass_body, inner, names, params);
    }

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_diag(glass_body, names.phi_laplace_diag);
    zero_diag.exec();
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_diag(
            inner, names.sigma, names.phi_laplace_diag, params.pair_weight_regularization_);
        compute_diag.exec();
    }
    else
    {
        computeOpheliePhiOperatorDiagonal<ExecutionPolicy>(glass_body, inner, names, params);
    }
    StateDynamics<ExecutionPolicy, OphelieJacobiPhiImagUpdateCK> jacobi_update(
        glass_body, names, params.phi_gauge_penalty_, params.phi_post_jacobi_relaxation_);

    Real max_residual = -1.0;
    for (size_t iter = 0; iter < params.phi_post_jacobi_refinement_iterations_; ++iter)
    {
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        jacobi_update.exec();
        if ((iter + 1) % 25 == 0 || iter + 1 == params.phi_post_jacobi_refinement_iterations_)
        {
            applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
            max_residual = evaluateOpheliePhiImagJacobiRelativeResidual<ExecutionPolicy>(
                glass_body.getBaseParticles(), names, params.phi_gauge_penalty_);
            if (!std::isfinite(max_residual) || max_residual < params.phi_jacobi_tolerance_)
            {
                break;
            }
            syncVariableToDevice<Real>(glass_body.getBaseParticles(), names.phi_imag);
        }
    }
    if (max_residual < 0.0)
    {
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        max_residual = evaluateOpheliePhiImagJacobiRelativeResidual<ExecutionPolicy>(
            glass_body.getBaseParticles(), names, params.phi_gauge_penalty_);
    }
    return max_residual;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#include "electromagnetic_ophelie_phi.hpp"
#endif // ELECTROMAGNETIC_OPHELIE_PHI_H
