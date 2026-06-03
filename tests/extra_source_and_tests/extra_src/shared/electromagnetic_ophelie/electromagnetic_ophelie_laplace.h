#ifndef ELECTROMAGNETIC_OPHELIE_LAPLACE_H
#define ELECTROMAGNETIC_OPHELIE_LAPLACE_H

#include "interaction_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real harmonicMean(Real lhs, Real rhs) { return Real(2) * lhs * rhs / (lhs + rhs + TinyReal); }

inline Real pairwiseNegativeLaplaceWeight(Real dW_ijV_j, Real distance, Real distance_sq,
                                          Real pair_weight_regularization, Real reference_smoothing_length)
{
    const Real h_sq = reference_smoothing_length * reference_smoothing_length;
    return -Real(2) * distance * dW_ijV_j / (distance_sq + pair_weight_regularization * h_sq + TinyReal);
}

inline Vecd pairwiseGradientWeightUncorrected(Real dW_ijV_j, const Vecd &e_ij) { return -dW_ijV_j * e_ij; }

template <typename... RelationTypes>
class OpheliePairwiseLaplaceCK;

template <template <typename...> class RelationType, typename... Parameters>
class OpheliePairwiseLaplaceCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    explicit OpheliePairwiseLaplaceCK(RelationType<Parameters...> &relation, const std::string &input_name,
                                      const std::string &coefficient_name, const std::string &output_name,
                                      Real pair_weight_regularization = Real(0.01))
        : BaseInteraction(relation), pair_weight_regularization_(pair_weight_regularization),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_input_(this->particles_->template getVariableByName<Real>(input_name)),
          dv_output_(this->particles_->template getVariableByName<Real>(output_name)),
          dv_coefficient_(this->particles_->template getVariableByName<Real>(coefficient_name))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              input_(encloser.dv_input_->DelegatedData(ex_policy)),
              output_(encloser.dv_output_->DelegatedData(ex_policy)),
              coefficient_(encloser.dv_coefficient_->DelegatedData(ex_policy)),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = input_[index_i];
            const Real coefficient_i = coefficient_[index_i];
            Real laplace_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real harmonic_coefficient = harmonicMean(coefficient_i, coefficient_[index_j]);
                const Real pair_weight = harmonic_coefficient * pairwiseNegativeLaplaceWeight(
                                                                  dW_ijV_j, distance, distance_sq,
                                                                  pair_weight_regularization_, reference_smoothing_length_);
                laplace_i += pair_weight * (phi_i - input_[index_j]);
            }
            output_[index_i] = laplace_i;
        }

      protected:
        Real *Vol_;
        Real *input_;
        Real *output_;
        Real *coefficient_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_input_;
    DiscreteVariable<Real> *dv_output_;
    DiscreteVariable<Real> *dv_coefficient_;
};

template <typename... Parameters>
class OpheliePairwiseLaplaceCK<Inner<Parameters...>>
    : public OpheliePairwiseLaplaceCK<Base, Inner<Parameters...>>
{
  public:
    explicit OpheliePairwiseLaplaceCK(Inner<Parameters...> &inner_relation, const std::string &input_name,
                                    const std::string &coefficient_name, const std::string &output_name,
                                    Real pair_weight_regularization = Real(0.01))
        : OpheliePairwiseLaplaceCK<Base, Inner<Parameters...>>(inner_relation, input_name, coefficient_name, output_name,
                                                                 pair_weight_regularization)
    {
    }
};

template <typename... RelationTypes>
class OpheliePairwiseLaplaceDiagonalCK;

template <template <typename...> class RelationType, typename... Parameters>
class OpheliePairwiseLaplaceDiagonalCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    explicit OpheliePairwiseLaplaceDiagonalCK(RelationType<Parameters...> &relation, const std::string &coefficient_name,
                                              const std::string &diagonal_name, Real pair_weight_regularization = Real(0.01))
        : BaseInteraction(relation), pair_weight_regularization_(pair_weight_regularization),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_diagonal_(this->particles_->template getVariableByName<Real>(diagonal_name)),
          dv_coefficient_(this->particles_->template getVariableByName<Real>(coefficient_name))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              diagonal_(encloser.dv_diagonal_->DelegatedData(ex_policy)),
              coefficient_(encloser.dv_coefficient_->DelegatedData(ex_policy)),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real coefficient_i = coefficient_[index_i];
            Real diag_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real harmonic_coefficient = harmonicMean(coefficient_i, coefficient_[index_j]);
                diag_i += harmonic_coefficient * pairwiseNegativeLaplaceWeight(
                                                    dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                    reference_smoothing_length_);
            }
            diagonal_[index_i] = diag_i;
        }

      protected:
        Real *Vol_;
        Real *diagonal_;
        Real *coefficient_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_diagonal_;
    DiscreteVariable<Real> *dv_coefficient_;
};

template <typename... Parameters>
class OpheliePairwiseLaplaceDiagonalCK<Inner<Parameters...>>
    : public OpheliePairwiseLaplaceDiagonalCK<Base, Inner<Parameters...>>
{
  public:
    explicit OpheliePairwiseLaplaceDiagonalCK(Inner<Parameters...> &inner_relation, const std::string &coefficient_name,
                                            const std::string &diagonal_name, Real pair_weight_regularization = Real(0.01))
        : OpheliePairwiseLaplaceDiagonalCK<Base, Inner<Parameters...>>(inner_relation, coefficient_name, diagonal_name,
                                                                       pair_weight_regularization)
    {
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_LAPLACE_H
