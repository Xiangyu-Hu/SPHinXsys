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

/**
 * Contact (e.g. fluid–wall) Laplace contribution. Accumulates into the source body's output
 * field (Inner must run first and assign). Contact body must register the same input/coefficient
 * variable names. Use k=0 on adiabatic wall faces to suppress heat exchange.
 */
template <typename... Parameters>
class OpheliePairwiseLaplaceCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit OpheliePairwiseLaplaceCK(Contact<Parameters...> &contact_relation, const std::string &input_name,
                                      const std::string &coefficient_name, const std::string &output_name,
                                      Real pair_weight_regularization = Real(0.01),
                                      bool dirichlet_ghost = false)
        : BaseInteraction(contact_relation), pair_weight_regularization_(pair_weight_regularization),
          dirichlet_ghost_(dirichlet_ghost),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_input_(this->particles_->template getVariableByName<Real>(input_name)),
          dv_output_(this->particles_->template getVariableByName<Real>(output_name)),
          dv_coefficient_(this->particles_->template getVariableByName<Real>(coefficient_name))
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            dv_contact_input_.push_back(
                this->contact_particles_[k]->template getVariableByName<Real>(input_name));
            dv_contact_coefficient_.push_back(
                this->contact_particles_[k]->template getVariableByName<Real>(coefficient_name));
        }
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
              input_(encloser.dv_input_->DelegatedData(ex_policy)),
              contact_input_(encloser.dv_contact_input_[contact_index]->DelegatedData(ex_policy)),
              output_(encloser.dv_output_->DelegatedData(ex_policy)),
              coefficient_(encloser.dv_coefficient_->DelegatedData(ex_policy)),
              contact_coefficient_(encloser.dv_contact_coefficient_[contact_index]->DelegatedData(ex_policy)),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              dirichlet_ghost_(encloser.dirichlet_ghost_),
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
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
                const Real harmonic_coefficient = harmonicMean(coefficient_i, contact_coefficient_[index_j]);
                const Real pair_weight = harmonic_coefficient * pairwiseNegativeLaplaceWeight(
                                                                  dW_ijV_j, distance, distance_sq,
                                                                  pair_weight_regularization_, reference_smoothing_length_);
                // Dirichlet face at the geometric wall: dummy φ_j = 2 φ_bc − φ_i (Adami).
                // Stored contact_input is φ_bc. Face-value φ_j=φ_bc halves (φ_i−φ_j) on a linear field.
                const Real phi_j = dirichlet_ghost_ ? (Real(2) * contact_input_[index_j] - phi_i)
                                                    : contact_input_[index_j];
                laplace_i += pair_weight * (phi_i - phi_j);
            }
            output_[index_i] += laplace_i;
        }

      protected:
        Real *Vol_;
        Real *contact_Vol_;
        Real *input_;
        Real *contact_input_;
        Real *output_;
        Real *coefficient_;
        Real *contact_coefficient_;
        Real pair_weight_regularization_;
        bool dirichlet_ghost_;
        Real reference_smoothing_length_;
    };

  protected:
    Real pair_weight_regularization_;
    bool dirichlet_ghost_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_input_;
    DiscreteVariable<Real> *dv_output_;
    DiscreteVariable<Real> *dv_coefficient_;
    StdVec<DiscreteVariable<Real> *> dv_contact_input_;
    StdVec<DiscreteVariable<Real> *> dv_contact_coefficient_;
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

/** Jacobi/GMRES preconditioner diagonal estimate for div(sigma grad phi) with uncorrected grad/div kernels. */
template <typename... RelationTypes>
class OphelieDivSigmaGradDiagonalCK;

template <template <typename...> class RelationType, typename... Parameters>
class OphelieDivSigmaGradDiagonalCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    explicit OphelieDivSigmaGradDiagonalCK(RelationType<Parameters...> &relation, const std::string &coefficient_name,
                                           const std::string &diagonal_name)
        : BaseInteraction(relation),
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
              coefficient_(encloser.dv_coefficient_->DelegatedData(ex_policy))
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
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
                const Real harmonic_coefficient = harmonicMean(coefficient_i, coefficient_[index_j]);
                diag_i += harmonic_coefficient * g_ij.squaredNorm();
            }
            diagonal_[index_i] = diag_i;
        }

      protected:
        Real *Vol_;
        Real *diagonal_;
        Real *coefficient_;
    };

  protected:
    DiscreteVariable<Real> *dv_diagonal_;
    DiscreteVariable<Real> *dv_coefficient_;
};

template <typename... Parameters>
class OphelieDivSigmaGradDiagonalCK<Inner<Parameters...>>
    : public OphelieDivSigmaGradDiagonalCK<Base, Inner<Parameters...>>
{
  public:
    explicit OphelieDivSigmaGradDiagonalCK(Inner<Parameters...> &inner_relation, const std::string &coefficient_name,
                                           const std::string &diagonal_name)
        : OphelieDivSigmaGradDiagonalCK<Base, Inner<Parameters...>>(inner_relation, coefficient_name, diagonal_name)
    {
    }
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_LAPLACE_H
