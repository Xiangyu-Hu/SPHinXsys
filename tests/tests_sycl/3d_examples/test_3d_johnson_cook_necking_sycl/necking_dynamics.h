#ifndef NECKING_DYNAMICS_CK_H
#define NECKING_DYNAMICS_CK_H

#include "necking_case.h"

namespace SPH
{
/** Prescribed grips share exactly the classic example's displacement and ramp. */
class PrescribedNeckingGripsCK : public LocalDynamics
{
  public:
    PrescribedNeckingGripsCK(SPHBody &body, GripMotion &motion)
        : LocalDynamics(body), motion_(motion),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
          dv_pos0_(particles_->getVariableByName<Vecd>("InitialPosition")),
          sv_motion_(particles_->addUniqueSingleVariable<GripMotion>("GripMotion")) {}

    void setupDynamics(Real = 0) override { sv_motion_->setValue(motion_); }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
              motion_(encloser.sv_motion_->DelegatedData(ex_policy)) {}

        void update(size_t i, Real = 0)
        {
            if (math::fabs(pos0_[i][2]) < NeckingRod::grip_position)
                return;
            const Real side = pos0_[i][2] > 0 ? Real(1) : Real(-1);
            vel_[i] = Vecd(0, 0, side * motion_->velocity);
            if (motion_->constrain_position)
            {
                pos_[i] = pos0_[i];
                pos_[i][2] += side * motion_->displacement;
            }
        }

      protected:
        Vecd *pos_, *vel_, *pos0_;
        GripMotion *motion_;
    };

  protected:
    GripMotion &motion_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_, *dv_pos0_;
    SingleVariable<GripMotion> *sv_motion_;
};

/** Reduce on the execution device; transfer only two boundary reactions per step. */
class NeckingGripReactionCK : public LocalDynamicsReduce<ReduceSum<Vec2d>>
{
  public:
    explicit NeckingGripReactionCK(SPHBody &body)
        : LocalDynamicsReduce<ReduceSum<Vec2d>>(body),
          dv_force_(particles_->getVariableByName<Vecd>("Force")),
          dv_pos0_(particles_->getVariableByName<Vecd>("InitialPosition")),
          dv_mass_(particles_->getVariableByName<Real>("Mass")) {}

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : force_(encloser.dv_force_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
              mass_(encloser.dv_mass_->DelegatedData(ex_policy)) {}

        Vec2d reduce(size_t i, Real acceleration)
        {
            if (pos0_[i][2] >= NeckingRod::grip_position)
                return Vec2d(mass_[i] * acceleration - force_[i][2], 0);
            if (pos0_[i][2] <= -NeckingRod::grip_position)
                return Vec2d(0, -mass_[i] * acceleration - force_[i][2]);
            return Vec2d::Zero();
        }

      protected:
        Vecd *force_, *pos0_;
        Real *mass_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_force_, *dv_pos0_;
    DiscreteVariable<Real> *dv_mass_;
};
} // namespace SPH
#endif // NECKING_DYNAMICS_CK_H
