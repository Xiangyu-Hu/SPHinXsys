/**
 * @file electromagnetic_ophelie_phi_variational.h
 * @brief Variational φ operator K=G*(σG) (historically Doc214 E2): K=G^*(σG), b=−G^*(σ a_sign ω A) with exact discrete adjoint of SPH G.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_PHI_VARIATIONAL_H
#define ELECTROMAGNETIC_OPHELIE_PHI_VARIATIONAL_H

#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_phi_gradient.h"
#include "interaction_algorithms_ck.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** out = scale * in (Vecd). */
class OphelieScaleVecdFieldCK : public LocalDynamics
{
  public:
    OphelieScaleVecdFieldCK(SPHBody &sph_body, const std::string &field_name, Real scale_factor)
        : LocalDynamics(sph_body), scale_factor_(scale_factor),
          dv_field_(particles_->template getVariableByName<Vecd>(field_name))
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
        Vecd *field_;
    };

  protected:
    Real scale_factor_;
    DiscreteVariable<Vecd> *dv_field_;
};

class ZeroOphelieScalarFieldForE2CK : public LocalDynamics
{
  public:
    ZeroOphelieScalarFieldForE2CK(SPHBody &sph_body, const std::string &field_name)
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

class OphelieScaleScalarFieldForE2CK : public LocalDynamics
{
  public:
    OphelieScaleScalarFieldForE2CK(SPHBody &sph_body, const std::string &field_name, Real scale_factor)
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

class OphelieScaleVecdBySigmaForE2CK : public LocalDynamics
{
  public:
    OphelieScaleVecdBySigmaForE2CK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
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

class OpheliePhiGaugePenaltyForE2CK : public LocalDynamics
{
  public:
    OpheliePhiGaugePenaltyForE2CK(SPHBody &sph_body, const std::string &phi_field, const std::string &phi_lhs_field,
                                  Real gauge_penalty)
        : LocalDynamics(sph_body), gauge_penalty_(gauge_penalty),
          dv_phi_(particles_->template getVariableByName<Real>(phi_field)),
          dv_lhs_(particles_->template getVariableByName<Real>(phi_lhs_field))
    {
    }
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : gauge_penalty_(encloser.gauge_penalty_), phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
              lhs_(encloser.dv_lhs_->DelegatedData(ex_policy))
        {
        }
        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            lhs_[index_i] += gauge_penalty_ * phi_[index_i];
        }

      protected:
        Real gauge_penalty_;
        Real *phi_;
        Real *lhs_;
    };

  protected:
    Real gauge_penalty_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Real> *dv_lhs_;
};

/**
 * Gφ_i = Σ_j g_ij (φ_i − φ_j).
 * Doc214 G2: optional linear correction C_i on e_ij (same as ScalarPhiGradientCorrected).
 */
template <typename... RelationTypes>
class ComputeOphelieE2ScalarGradientCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieE2ScalarGradientCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieE2ScalarGradientCK(RelationType<Parameters...> &relation, const std::string &phi_field,
                                     const std::string &grad_field, const OphelieGlassFieldNames &names,
                                     bool use_linear_correction = false)
        : BaseInteraction(relation), use_linear_correction_(use_linear_correction),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(phi_field)),
          dv_grad_(this->particles_->template getVariableByName<Vecd>(grad_field)),
          dv_linear_correction_(this->particles_->template getVariableByName<Matd>(names.phi_grad_linear_correction))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
              grad_(encloser.dv_grad_->DelegatedData(ex_policy)),
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy)),
              use_linear_correction_(encloser.use_linear_correction_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            Vecd grad_phi = Vecd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Vecd e = use_linear_correction_ ? (linear_correction_[index_i] * this->e_ij(index_i, index_j))
                                                       : this->e_ij(index_i, index_j);
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, e);
                grad_phi += g_ij * (phi_i - phi_[index_j]);
            }
            grad_[index_i] = grad_phi;
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *grad_;
        Matd *linear_correction_;
        bool use_linear_correction_;
    };

  protected:
    bool use_linear_correction_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_grad_;
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... Parameters>
class ComputeOphelieE2ScalarGradientCK<Inner<Parameters...>>
    : public ComputeOphelieE2ScalarGradientCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieE2ScalarGradientCK(Inner<Parameters...> &inner_relation, const std::string &phi_field,
                                     const std::string &grad_field, const OphelieGlassFieldNames &names,
                                     bool use_linear_correction = false)
        : ComputeOphelieE2ScalarGradientCK<Base, Inner<Parameters...>>(inner_relation, phi_field, grad_field, names,
                                                                         use_linear_correction)
    {
    }
};

template <typename... RelationTypes>
class ComputeOphelieAdjointGradientDotCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieAdjointGradientDotCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieAdjointGradientDotCK(RelationType<Parameters...> &relation, const std::string &vec_field,
                                       const std::string &scalar_out_field, const OphelieGlassFieldNames &names,
                                       bool use_linear_correction = false)
        : BaseInteraction(relation), use_linear_correction_(use_linear_correction),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_vec_(this->particles_->template getVariableByName<Vecd>(vec_field)),
          dv_out_(this->particles_->template getVariableByName<Real>(scalar_out_field)),
          dv_linear_correction_(this->particles_->template getVariableByName<Matd>(names.phi_grad_linear_correction))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              vec_(encloser.dv_vec_->DelegatedData(ex_policy)),
              out_(encloser.dv_out_->DelegatedData(ex_policy)),
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy)),
              use_linear_correction_(encloser.use_linear_correction_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            // Doc215 E1: true discrete adjoint of
            //   (Gφ)_i = Σ_j g_ij(C_i,Vol_j) (φ_i − φ_j)
            // via gather of the scatter transpose:
            //   for each directed row i→j:  out_i += g_ij·v_i ; out_j -= g_ij·v_i
            // which at particle i is Σ_j [ g_ij·v_i − g_ji·v_j ].
            // Doc214 reused C_i for both ends (≡ Hold), wrong when C_i≠C_j.
            // Uncorrected + symmetric kernel coincides with the old
            // g_ij·(v_i+(Vol_i/Vol_j)v_j) form; corrected G must use C_j on g_ji.
            const Vecd v_i = vec_[index_i];
            const Real Vol_i = Vol_[index_i];
            Real acc = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Real Vol_j = Vol_[index_j];
                if (Vol_j <= TinyReal || Vol_i <= TinyReal)
                {
                    continue;
                }
                const Real dW_ij = this->dW_ij(index_i, index_j);
                const Vecd e_ij = this->e_ij(index_i, index_j);
                const Vecd e_i = use_linear_correction_ ? (linear_correction_[index_i] * e_ij) : e_ij;
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ij * Vol_j, e_i);

                const Real dW_ji = this->dW_ij(index_j, index_i);
                const Vecd e_ji = this->e_ij(index_j, index_i);
                const Vecd e_j = use_linear_correction_ ? (linear_correction_[index_j] * e_ji) : e_ji;
                const Vecd g_ji = pairwiseGradientWeightUncorrected(dW_ji * Vol_i, e_j);

                acc += g_ij.dot(v_i) - g_ji.dot(vec_[index_j]);
            }
            out_[index_i] += acc;
        }

      protected:
        Real *Vol_;
        Vecd *vec_;
        Real *out_;
        Matd *linear_correction_;
        bool use_linear_correction_;
    };

  protected:
    bool use_linear_correction_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_vec_;
    DiscreteVariable<Real> *dv_out_;
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... Parameters>
class ComputeOphelieAdjointGradientDotCK<Inner<Parameters...>>
    : public ComputeOphelieAdjointGradientDotCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieAdjointGradientDotCK(Inner<Parameters...> &inner_relation, const std::string &vec_field,
                                       const std::string &scalar_out_field, const OphelieGlassFieldNames &names,
                                       bool use_linear_correction = false)
        : ComputeOphelieAdjointGradientDotCK<Base, Inner<Parameters...>>(inner_relation, vec_field, scalar_out_field,
                                                                           names, use_linear_correction)
    {
    }
};

template <class ExecutionPolicy>
inline void execOphelieAdjointGradientDot(Inner<> &inner, const OphelieGlassFieldNames &names,
                                          const std::string &vec_field, const std::string &scalar_out_field,
                                          bool use_linear_correction = false)
{
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieAdjointGradientDotCK<Inner<>>> op(
        inner, vec_field, scalar_out_field, names, use_linear_correction);
    op.exec();
}

inline bool ophelieDoc214E2UseCorrectedG(const OphelieParameters &params)
{
    return opheliePhiUseCorrectedGradient(params);
}

/** Doc214 G2b: use paired D_c as G^* with G^* ≈ −D_c (compatible package). */
inline bool ophelieDoc214E2UseCompatibleAdjoint(const OphelieParameters &params)
{
    return opheliePhiUseCorrectedDivergence(params);
}

/** b = −G^* (σ (a_sign ω A)). Compatible: G^*=−D_c ⇒ b = D_c(σ a_sign ω A). */
template <class ExecutionPolicy>
inline void setupOpheliePhiRhsDoc214E2Variational(SolidBody &glass_body, Inner<> &inner,
                                                  const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                                  const std::string &a_field, const std::string &sigma_a_scratch,
                                                  const std::string &phi_rhs_field, Real a_sign, Real omega)
{
    const bool use_corr = ophelieDoc214E2UseCorrectedG(params);
    const bool use_dc = ophelieDoc214E2UseCompatibleAdjoint(params);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldForE2CK> zero_rhs(glass_body, phi_rhs_field);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdBySigmaForE2CK> scale_sigma_a(glass_body, names, a_field,
                                                                                   sigma_a_scratch);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdFieldCK> scale_v(glass_body, sigma_a_scratch, a_sign * omega);
    StateDynamics<ExecutionPolicy, OphelieScaleScalarFieldForE2CK> neg_rhs(glass_body, phi_rhs_field, Real(-1));

    zero_rhs.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (use_corr)
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params, Real(1.0e-6), true);
    }
    scale_sigma_a.exec();
    scale_v.exec();
    if (use_dc)
    {
        // b = D_c(σ a_sign ω A)  (= −G^* … when G^* = −D_c)
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCorrectedCK<Inner<>>> compute_div(
            inner, sigma_a_scratch, phi_rhs_field, names);
        compute_div.exec();
    }
    else
    {
        execOphelieAdjointGradientDot<ExecutionPolicy>(inner, names, sigma_a_scratch, phi_rhs_field, use_corr);
        neg_rhs.exec();
    }
}

/**
 * Same-G recon: E = -(Gφ + a_sign ω A), J=σE, Q=½σ|E|² into edge-recon field names
 * (so downstream sync/primary copy stays unchanged).
 */
template <class ExecutionPolicy>
inline void execOphelieDoc214E2ParticleGReconForComponent(SolidBody &glass_body, Inner<> &inner,
                                                          const OphelieGlassFieldNames &names,
                                                          const OphelieParameters &params, const std::string &phi_field,
                                                          const std::string &grad_scratch,
                                                          const std::string &active_a_field, const std::string &e_field,
                                                          const std::string &j_field, const std::string &q_field,
                                                          Real a_sign, Real omega)
{
    const bool use_corr = ophelieDoc214E2UseCorrectedG(params);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (use_corr)
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params, Real(1.0e-6), false);
    }
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieE2ScalarGradientCK<Inner<>>> compute_grad(
        inner, phi_field, grad_scratch, names, use_corr);
    compute_grad.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, grad_scratch);
    syncVariableToHost<Vecd>(particles, active_a_field);
    syncVariableToHost<Real>(particles, names.sigma);
    Vecd *grad = particles.getVariableDataByName<Vecd>(grad_scratch);
    Vecd *a = particles.getVariableDataByName<Vecd>(active_a_field);
    Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    Vecd *e = particles.getVariableDataByName<Vecd>(e_field);
    Vecd *j = particles.getVariableDataByName<Vecd>(j_field);
    Real *q = particles.getVariableDataByName<Real>(q_field);
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd e_i = -(grad[i] + a_sign * omega * a[i]);
        e[i] = e_i;
        j[i] = sigma[i] * e_i;
        q[i] = Real(0.5) * sigma[i] * e_i.squaredNorm();
    }
    syncVariableToDevice<Vecd>(particles, e_field);
    syncVariableToDevice<Vecd>(particles, j_field);
    syncVariableToDevice<Real>(particles, q_field);
}

template <class ExecutionPolicy>
inline void execOphelieDoc214E2ParticleGRecon(SolidBody &glass_body, Inner<> &inner,
                                              const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                              bool skip_real = false)
{
    std::cout << "[ophelie][doc214][e2] recon=particle_G E=-(G phi + a_sign omega A) (same G as solve)"
              << " corrected_G=" << (ophelieDoc214E2UseCorrectedG(params) ? 1 : 0)
              << " compatible_adj=" << (ophelieDoc214E2UseCompatibleAdjoint(params) ? 1 : 0) << std::endl;
    execOphelieDoc214E2ParticleGReconForComponent<ExecutionPolicy>(
        glass_body, inner, names, params, names.phi_imag, names.grad_phi_imag, names.a_src_real,
        names.e_edge_recon_imag, names.j_edge_recon_imag, names.joule_heat_edge_recon_imag, Real(1), params.omega());
    if (!skip_real && params.edge_flux_complex_)
    {
        execOphelieDoc214E2ParticleGReconForComponent<ExecutionPolicy>(
            glass_body, inner, names, params, names.phi_real, names.grad_phi_real, names.a_src_imag,
            names.e_edge_recon_real, names.j_edge_recon_real, names.joule_heat_edge_recon_real, Real(-1),
            params.omega());
        // complex Q used by cylinder compare
        BaseParticles &particles = glass_body.getBaseParticles();
        const size_t n = particles.TotalRealParticles();
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_imag);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_real);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        const Real *qi = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);
        const Real *qr = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_real);
        Real *qc = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        for (size_t i = 0; i < n; ++i)
        {
            qc[i] = qi[i] + qr[i];
        }
        syncVariableToDevice<Real>(particles, names.joule_heat_edge_recon_complex);
    }
}

/** Kφ = G^* (σ G φ) + gauge. Compatible: G^*=−D_c ⇒ Kφ = −D_c(σ G_c φ). */
template <class ExecutionPolicy>
inline void applyOpheliePhiLhsDoc214E2Variational(SolidBody &glass_body, Inner<> &inner,
                                                  const OphelieGlassFieldNames &names, const OphelieParameters &params,
                                                  const std::string &phi_field, const std::string &grad_scratch,
                                                  const std::string &sigma_grad_scratch, const std::string &phi_lhs_field,
                                                  Real gauge_penalty)
{
    const bool use_corr = ophelieDoc214E2UseCorrectedG(params);
    const bool use_dc = ophelieDoc214E2UseCompatibleAdjoint(params);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldForE2CK> zero_lhs(glass_body, phi_lhs_field);
    StateDynamics<ExecutionPolicy, OphelieScaleVecdBySigmaForE2CK> scale_sigma_grad(glass_body, names, grad_scratch,
                                                                                      sigma_grad_scratch);
    StateDynamics<ExecutionPolicy, OpheliePhiGaugePenaltyForE2CK> apply_gauge(glass_body, phi_field, phi_lhs_field,
                                                                               gauge_penalty);
    StateDynamics<ExecutionPolicy, OphelieScaleScalarFieldForE2CK> neg_lhs(glass_body, phi_lhs_field, Real(-1));

    update_cell_linked_list.exec();
    update_inner_relation.exec();
    zero_lhs.exec();
    if (use_corr)
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params, Real(1.0e-6), false);
    }
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieE2ScalarGradientCK<Inner<>>> compute_grad(
        inner, phi_field, grad_scratch, names, use_corr);
    compute_grad.exec();
    scale_sigma_grad.exec();
    if (use_dc)
    {
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCorrectedCK<Inner<>>> compute_div(
            inner, sigma_grad_scratch, phi_lhs_field, names);
        compute_div.exec();
        neg_lhs.exec(); // K = −D_c(σ G_c φ)
    }
    else
    {
        execOphelieAdjointGradientDot<ExecutionPolicy>(inner, names, sigma_grad_scratch, phi_lhs_field, use_corr);
    }
    apply_gauge.exec();
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_VARIATIONAL_H
