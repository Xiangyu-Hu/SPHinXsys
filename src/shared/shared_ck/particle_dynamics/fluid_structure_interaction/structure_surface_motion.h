/**
 * @file    structure_surface_motion.h
 * @brief   Surface motion of a structure as seen by the surrounding fluid.
 *          The position is recorded before the structure sub loop and the
 *          averaged velocity and acceleration over that interval are recovered
 *          afterwards, so the fluid sees a consistent moving boundary.
 * @author  Pruthvik Arasikere Mallikarjuna and Xiangyu Hu
 */

#ifndef STRUCTURE_SURFACE_MOTION_H
#define STRUCTURE_SURFACE_MOTION_H

#include "base_local_dynamics.h"
#include "sph_system.h"

namespace SPH
{
namespace FSI
{
//----------------------------------------------------------------------
//      Record the position before the structure sub loop.
//----------------------------------------------------------------------
class InitializeDisplacementCK : public LocalDynamics
{
  public:
    explicit InitializeDisplacementCK(RealBody &solid_body)
        : LocalDynamics(solid_body),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_pos_temp_(particles_->registerStateVariable<Vecd>("TemporaryPosition")) {}

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, InitializeDisplacementCK &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos_temp_(encloser.dv_pos_temp_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0) { pos_temp_[index_i] = pos_[index_i]; }

      protected:
        Vecd *pos_, *pos_temp_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos_temp_;
};
//----------------------------------------------------------------------
//      Averaged velocity and acceleration over the last coupling interval.
//----------------------------------------------------------------------
class UpdateAverageVelocityAndAccelerationCK : public LocalDynamics
{
  public:
    explicit UpdateAverageVelocityAndAccelerationCK(RealBody &solid_body)
        : LocalDynamics(solid_body),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_pos_temp_(particles_->getVariableByName<Vecd>("TemporaryPosition")),
          dv_vel_ave_(particles_->registerStateVariable<Vecd>("AverageVelocity")),
          dv_acc_ave_(particles_->registerStateVariable<Vecd>("AverageAcceleration"))
    {
        // These drive the fluid-side FSI force directly, so restoring the
        // exact restart values is required; not persisting them would leave
        // the coupling at zero right after resume.
        particles_->addEvolvingVariable<Vecd>("AverageVelocity");
        particles_->addEvolvingVariable<Vecd>("AverageAcceleration");
    }

    struct UpdateKernel
    {
        template <typename ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, UpdateAverageVelocityAndAccelerationCK &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos_temp_(encloser.dv_pos_temp_->DelegatedData(ex_policy)),
              vel_ave_(encloser.dv_vel_ave_->DelegatedData(ex_policy)),
              acc_ave_(encloser.dv_acc_ave_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            Vecd updated_vel_ave = (pos_[index_i] - pos_temp_[index_i]) / (dt + Eps);
            acc_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
            vel_ave_[index_i] = updated_vel_ave;
        }

      protected:
        Vecd *pos_, *pos_temp_, *vel_ave_, *acc_ave_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos_temp_, *dv_vel_ave_, *dv_acc_ave_;
};
} // namespace FSI
} // namespace SPH

#endif // STRUCTURE_SURFACE_MOTION_H