/**
 * @file spdb_velocity_decay_ck.h
 * @brief Simple CK velocity decay for quasi-static relaxation on device.
 * @details Library has no CK pairwise damping. This is a demo-side substitute
 *          so SYCL acoustic loops can calm without host DampingPairwiseInner.
 */
#ifndef SPDB_VELOCITY_DECAY_CK_H
#define SPDB_VELOCITY_DECAY_CK_H

#include "base_local_dynamics.h"

namespace SPH
{

template <class DynamicsIdentifier>
class VelocityDecayCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    VelocityDecayCK(DynamicsIdentifier &identifier, Real decay_per_s)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
          decay_per_s_(decay_per_s) {}
    virtual ~VelocityDecayCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : vel_(encloser.dv_vel_->DelegatedData(ex_policy)), decay_per_s_(encloser.decay_per_s_) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            const Real factor = Real(1) / (Real(1) + decay_per_s_ * dt);
            vel_[index_i] *= factor;
        }

      protected:
        Vecd *vel_;
        Real decay_per_s_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_vel_;
    Real decay_per_s_;
};

} // namespace SPH

#endif
