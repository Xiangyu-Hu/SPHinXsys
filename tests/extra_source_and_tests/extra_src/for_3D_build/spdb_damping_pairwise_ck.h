/**
 * @file spdb_damping_pairwise_ck.h
 * @brief CK/device pairwise velocity damping aligned with host
 *        DampingPairwiseInner + InteractionSplit (test_3d_sliding / CPU SPDB).
 * @details Host algorithm is Gauss–Seidel with cell splitting (updates i and j).
 *          This CK port uses a Jacobi sweep (update i only from frozen neighbors),
 *          which is parallel-safe on SYCL and has the same pairwise exchange form
 *          when both particles run the same step. Not artificial bulk viscosity;
 *          eta is the FixedDampingRate used by DampingPairwiseInner.
 *          Extra-test only; does not modify src/.
 * @warning Parallel Jacobi on SYCL raced in F3 elastic balance (vmax~300).
 *          Current programmed_temperature_sycl uses VelocityDecayCK instead.
 */
#ifndef SPDB_DAMPING_PAIRWISE_CK_H
#define SPDB_DAMPING_PAIRWISE_CK_H

#include "interaction_ck.hpp"

namespace SPH
{

template <typename...>
class DampingPairwiseCK;

template <typename... Parameters>
class DampingPairwiseCK<Inner<Parameters...>> : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit DampingPairwiseCK(DynamicsIdentifier &identifier, Real eta)
        : BaseInteraction(identifier),
          eta_(eta),
          dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
          dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
    virtual ~DampingPairwiseCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseInteraction::InteractKernel(ex_policy, encloser),
              eta_(encloser.eta_),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}

        void interact(size_t index_i, Real dt = 0.0)
        {
            const Real Vol_i = Vol_[index_i];
            const Real capacity_i = mass_[index_i]; // FixedDampingRate with c=1
            const Vecd vel_i = vel_[index_i];
            Vecd delta = Vecd::Zero();

            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd rij = this->vec_r_ij(index_i, index_j);
                const Real r_ij = rij.norm();
                if (r_ij <= Eps)
                    continue;

                const Real capacity_j = mass_[index_j];
                // Host Neighborhood::dW_ij_ is |dW/dr|; CK dW_ij may be signed.
                const Real dW = ABS(this->dW_ij(index_i, index_j));
                const Real parameter_b = eta_ * dW * Vol_i * Vol_[index_j] * dt / r_ij;
                const Real denom =
                    capacity_i * capacity_j - parameter_b * (capacity_i + capacity_j);
                if (ABS(denom) <= Eps)
                    continue;

                const Vecd pair_difference = vel_i - vel_[index_j];
                const Vecd increment = parameter_b * pair_difference / denom;
                delta += increment * capacity_j;
            }
            vel_[index_i] = vel_i + delta;
        }

      protected:
        Real eta_;
        Vecd *vel_;
        Real *mass_;
        Real *Vol_;
    };

  protected:
    Real eta_;
    DiscreteVariable<Vecd> *dv_vel_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Real> *dv_Vol_;
};

} // namespace SPH

#endif
