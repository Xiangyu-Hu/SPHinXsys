#ifndef ELECTROMAGNETIC_OPHELIE_PHI_GRADIENT_H
#define ELECTROMAGNETIC_OPHELIE_PHI_GRADIENT_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_parameters.h"
#include "interaction_ck.h"
#include "vector_functions.h"

#include <cmath>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** G_c: corrected gradient (--phi-compatible-correction or diagnostic --phi-gradient-correction). */
inline bool opheliePhiUseCorrectedGradient(const OphelieParameters &params)
{
    return params.phi_compatible_correction_ || params.phi_gradient_correction_;
}

/** D_c: corrected divergence (only with full compatible package). */
inline bool opheliePhiUseCorrectedDivergence(const OphelieParameters &params)
{
    return params.phi_compatible_correction_;
}

struct OpheliePhiGradCorrectionMatrixStats
{
    Real det_min = 0.0;
    Real det_max = 0.0;
    size_t det_negative_count = 0;
    size_t det_small_count = 0;
    size_t fallback_identity_count = 0;
    Real condition_proxy_max = 0.0;
};

inline OpheliePhiGradCorrectionMatrixStats computeOpheliePhiGradCorrectionMatrixStatsFromB(
    BaseParticles &particles, const std::string &field_name, size_t total_real_particles, Real det_min_threshold)
{
    syncVariableToHost<Matd>(particles, field_name);
    const Matd *B = particles.getVariableDataByName<Matd>(field_name);
    OpheliePhiGradCorrectionMatrixStats stats;
    if (total_real_particles == 0)
    {
        return stats;
    }

    stats.det_min = B[0].determinant();
    stats.det_max = stats.det_min;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real det = B[i].determinant();
        if (!std::isfinite(det))
        {
            stats.fallback_identity_count++;
            continue;
        }
        stats.det_min = std::min(stats.det_min, det);
        stats.det_max = std::max(stats.det_max, det);
        if (det <= Real(0))
        {
            stats.det_negative_count++;
        }
        if (std::abs(det) < det_min_threshold)
        {
            stats.det_small_count++;
        }
        const Matd BtB = B[i].transpose() * B[i];
        stats.condition_proxy_max =
            std::max(stats.condition_proxy_max, BtB.norm() / (std::abs(det) + TinyReal));
        if (det <= Real(0) || std::abs(det) < det_min_threshold)
        {
            stats.fallback_identity_count++;
        }
    }
    return stats;
}

inline void logOpheliePhiGradCorrectionMatrixStats(const OpheliePhiGradCorrectionMatrixStats &stats)
{
    std::cout << "[ophelie] phi_grad_B_stats: det_min=" << stats.det_min << " det_max=" << stats.det_max
              << " det_negative_count=" << stats.det_negative_count << " det_small_count=" << stats.det_small_count
              << " fallback_identity_count=" << stats.fallback_identity_count
              << " condition_proxy_max=" << stats.condition_proxy_max << std::endl;
}

template <typename... RelationTypes>
class ComputeOpheliePhiGradLinearCorrectionMatrixCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOpheliePhiGradLinearCorrectionMatrixCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOpheliePhiGradLinearCorrectionMatrixCK(RelationType<Parameters...> &relation,
                                                  const OphelieGlassFieldNames &names)
        : BaseInteraction(relation),
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
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            Matd local_configuration = Matd::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd gradW_ijV_j =
                    this->dW_ij(index_i, index_j) * Vol_[index_j] * this->e_ij(index_i, index_j);
                local_configuration -= this->vec_r_ij(index_i, index_j) * gradW_ijV_j.transpose();
            }
            linear_correction_[index_i] = local_configuration;
        }

      protected:
        Real *Vol_;
        Matd *linear_correction_;
    };

  protected:
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... Parameters>
class ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<Parameters...>>
    : public ComputeOpheliePhiGradLinearCorrectionMatrixCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOpheliePhiGradLinearCorrectionMatrixCK(Inner<Parameters...> &inner_relation,
                                                  const OphelieGlassFieldNames &names)
        : ComputeOpheliePhiGradLinearCorrectionMatrixCK<Base, Inner<Parameters...>>(inner_relation, names)
    {
    }
};

class InvertOpheliePhiGradLinearCorrectionMatrixCK : public LocalDynamics
{
  public:
    InvertOpheliePhiGradLinearCorrectionMatrixCK(RealBody &sph_body, const OphelieGlassFieldNames &names,
                                                 Real det_min_threshold = Real(1.0e-6))
        : LocalDynamics(sph_body), det_min_threshold_(det_min_threshold),
          dv_linear_correction_(particles_->template getVariableByName<Matd>(names.phi_grad_linear_correction))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : det_min_threshold_(encloser.det_min_threshold_),
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real det = linear_correction_[index_i].determinant();
            if (!std::isfinite(det) || det <= Real(0) || std::abs(det) < det_min_threshold_)
            {
                linear_correction_[index_i] = Matd::Identity();
                return;
            }
            linear_correction_[index_i] = inverseTikhonov(linear_correction_[index_i], SqrtEps);
        }

      protected:
        Real det_min_threshold_;
        Matd *linear_correction_;
    };

  protected:
    Real det_min_threshold_;
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... RelationTypes>
class ComputeOphelieScalarPhiGradientCorrectedCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieScalarPhiGradientCorrectedCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieScalarPhiGradientCorrectedCK(RelationType<Parameters...> &relation,
                                               const OphelieGlassFieldNames &names)
        : BaseInteraction(relation),
          dv_phi_imag_(this->particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_grad_phi_imag_(this->particles_->template getVariableByName<Vecd>(names.grad_phi_imag)),
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
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy)),
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy))
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
                const Vecd corrected_e_ij = linear_correction_[index_i] * this->e_ij(index_i, index_j);
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, corrected_e_ij);
                grad_phi += g_ij * (phi_i - phi_imag_[index_j]);
            }
            grad_phi_imag_[index_i] = grad_phi;
        }

      protected:
        Real *Vol_;
        Real *phi_imag_;
        Vecd *grad_phi_imag_;
        Matd *linear_correction_;
    };

  protected:
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... Parameters>
class ComputeOphelieScalarPhiGradientCorrectedCK<Inner<Parameters...>>
    : public ComputeOphelieScalarPhiGradientCorrectedCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieScalarPhiGradientCorrectedCK(Inner<Parameters...> &inner_relation,
                                               const OphelieGlassFieldNames &names)
        : ComputeOphelieScalarPhiGradientCorrectedCK<Base, Inner<Parameters...>>(inner_relation, names)
    {
    }
};

/** Uncorrected SPH div(F)_i = sum_j (F_j - F_i) · grad_i W_ij. */
template <typename... RelationTypes>
class ComputeOphelieVecdDivergenceCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieVecdDivergenceCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieVecdDivergenceCK(RelationType<Parameters...> &relation, const std::string &vec_field_name,
                                   const std::string &div_field_name)
        : BaseInteraction(relation), dv_vec_field_(this->particles_->template getVariableByName<Vecd>(vec_field_name)),
          dv_div_field_(this->particles_->template getVariableByName<Real>(div_field_name))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              vec_field_(encloser.dv_vec_field_->DelegatedData(ex_policy)),
              div_field_(encloser.dv_div_field_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd vec_i = vec_field_[index_i];
            Real div_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd grad_W_ij = pairwiseGradientWeightUncorrected(this->dW_ij(index_i, index_j) * Vol_[index_j],
                                                                         this->e_ij(index_i, index_j));
                div_i += (vec_field_[index_j] - vec_i).dot(grad_W_ij);
            }
            div_field_[index_i] = div_i;
        }

      protected:
        Real *Vol_;
        Vecd *vec_field_;
        Real *div_field_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_vec_field_;
    DiscreteVariable<Real> *dv_div_field_;
};

template <typename... Parameters>
class ComputeOphelieVecdDivergenceCK<Inner<Parameters...>>
    : public ComputeOphelieVecdDivergenceCK<Base, Inner<Parameters...>>
{
  public:
    explicit ComputeOphelieVecdDivergenceCK(Inner<Parameters...> &relation, const std::string &vec_field_name,
                                            const std::string &div_field_name)
        : ComputeOphelieVecdDivergenceCK<Base, Inner<Parameters...>>(relation, vec_field_name, div_field_name)
    {
    }
};

/** D_c(F): compatible divergence using the same C_i as G_c. */
template <typename... RelationTypes>
class ComputeOphelieVecdDivergenceCorrectedCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieVecdDivergenceCorrectedCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieVecdDivergenceCorrectedCK(RelationType<Parameters...> &relation, const std::string &vec_field_name,
                                            const std::string &div_field_name, const OphelieGlassFieldNames &names)
        : BaseInteraction(relation), dv_vec_field_(this->particles_->template getVariableByName<Vecd>(vec_field_name)),
          dv_div_field_(this->particles_->template getVariableByName<Real>(div_field_name)),
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
              vec_field_(encloser.dv_vec_field_->DelegatedData(ex_policy)),
              div_field_(encloser.dv_div_field_->DelegatedData(ex_policy)),
              linear_correction_(encloser.dv_linear_correction_->DelegatedData(ex_policy))
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd vec_i = vec_field_[index_i];
            Real div_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Vecd corrected_e_ij = linear_correction_[index_i] * this->e_ij(index_i, index_j);
                const Vecd grad_W_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, corrected_e_ij);
                div_i += (vec_field_[index_j] - vec_i).dot(grad_W_ij);
            }
            div_field_[index_i] = div_i;
        }

      protected:
        Real *Vol_;
        Vecd *vec_field_;
        Real *div_field_;
        Matd *linear_correction_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_vec_field_;
    DiscreteVariable<Real> *dv_div_field_;
    DiscreteVariable<Matd> *dv_linear_correction_;
};

template <typename... Parameters>
class ComputeOphelieVecdDivergenceCorrectedCK<Inner<Parameters...>>
    : public ComputeOphelieVecdDivergenceCorrectedCK<Base, Inner<Parameters...>>
{
  public:
    explicit ComputeOphelieVecdDivergenceCorrectedCK(Inner<Parameters...> &relation, const std::string &vec_field_name,
                                                     const std::string &div_field_name,
                                                     const OphelieGlassFieldNames &names)
        : ComputeOphelieVecdDivergenceCorrectedCK<Base, Inner<Parameters...>>(relation, vec_field_name, div_field_name,
                                                                                names)
    {
    }
};

template <class ExecutionPolicy>
inline OpheliePhiGradCorrectionMatrixStats execOpheliePhiGradCorrectionMatrixPrep(
    RealBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    Real det_min_threshold = Real(1.0e-6), bool log_stats = false)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    InteractionDynamicsCK<ExecutionPolicy, ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<>>> compute_B(inner,
                                                                                                             names);
    compute_B.exec();
    const OpheliePhiGradCorrectionMatrixStats stats =
        computeOpheliePhiGradCorrectionMatrixStatsFromB(particles, names.phi_grad_linear_correction, n,
                                                        det_min_threshold);
    StateDynamics<ExecutionPolicy, InvertOpheliePhiGradLinearCorrectionMatrixCK> invert_B(glass_body, names,
                                                                                            det_min_threshold);
    invert_B.exec();
    if (log_stats)
    {
        logOpheliePhiGradCorrectionMatrixStats(stats);
    }
    return stats;
}

template <class ExecutionPolicy>
inline void execOphelieScalarPhiGradientCorrectedOnly(Inner<> &inner, const OphelieGlassFieldNames &names)
{
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCorrectedCK<Inner<>>> compute_grad_phi(inner,
                                                                                                                 names);
    compute_grad_phi.exec();
}

template <class ExecutionPolicy>
inline void execOphelieVecdDivergenceUncorrected(Inner<> &inner, const std::string &vec_field_name,
                                                 const std::string &div_field_name)
{
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCK<Inner<>>> compute_div(inner, vec_field_name,
                                                                                              div_field_name);
    compute_div.exec();
}

template <class ExecutionPolicy>
inline void execOphelieVecdDivergenceCorrected(RealBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                               const std::string &vec_field_name, const std::string &div_field_name,
                                               Real det_min_threshold = Real(1.0e-6))
{
    execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, OphelieParameters{},
                                                            det_min_threshold, false);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCorrectedCK<Inner<>>> compute_div(
        inner, vec_field_name, div_field_name, names);
    compute_div.exec();
}

template <class ExecutionPolicy>
inline void execOphelieVecdDivergence(Inner<> &inner, const OphelieGlassFieldNames &names,
                                      const OphelieParameters &params, const std::string &vec_field_name,
                                      const std::string &div_field_name)
{
    if (opheliePhiUseCorrectedDivergence(params))
    {
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieVecdDivergenceCorrectedCK<Inner<>>> compute_div(
            inner, vec_field_name, div_field_name, names);
        compute_div.exec();
        return;
    }
    execOphelieVecdDivergenceUncorrected<ExecutionPolicy>(inner, vec_field_name, div_field_name);
}

template <class ExecutionPolicy, typename GradCK, typename CorrectedGradCK, typename CorrectionMatrixCK>
inline void execOphelieScalarPhiGradient(RealBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                         const OphelieParameters &params)
{
    if (opheliePhiUseCorrectedGradient(params))
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, inner, names, params);
        execOphelieScalarPhiGradientCorrectedOnly<ExecutionPolicy>(inner, names);
        return;
    }
    InteractionDynamicsCK<ExecutionPolicy, GradCK> compute_grad_phi(inner, names);
    compute_grad_phi.exec();
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_GRADIENT_H
