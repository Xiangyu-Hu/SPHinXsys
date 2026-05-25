/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * ...license header...
 * ------------------------------------------------------------------------- */
/**
 * @file    fluid_boundary_ck.h
 * @brief   CK-pattern GPU-compatible free-stream velocity correction.
 *          Parallel to fluid_boundary.h::FreeStreamVelocityCorrection but
 *          using DiscreteVariable/SingularVariable DelegatedData so the kernel
 *          runs correctly on device (SYCL) as well as on CPU.
 * @author  Pruthvik Arasikere Mallikarjuna and Xiangyu Hu
 */

#pragma once

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{

/**
 * @class FreeStreamVelocityCorrectionCK
 * @brief GPU-compatible (CK pattern) version of FreeStreamVelocityCorrection.
 *        Corrects stream-wise velocity of surface particles toward the far-field
 *        target, weighted by density ratio. Runs on device via StateDynamics.
 */
template <typename TargetVelocity>
class FreeStreamVelocityCorrectionCK : public LocalDynamics
{
  public:
    explicit FreeStreamVelocityCorrectionCK(SPHBody &sph_body,
                                            const Transform &transform = Transform())
        : LocalDynamics(sph_body),
          transform_(transform),
          rho0_(DynamicCast<Fluid>(this, sph_body_->getBaseMaterial()).ReferenceDensity()),
          sv_physical_time_(&sph_system_->svPhysicalTime()),
          dv_rho_sum_(particles_->template getVariableByName<Real>("DensitySummation")),
          dv_pos_(particles_->template getVariableByName<Vecd>("Position")),
          dv_vel_(particles_->template getVariableByName<Vecd>("Velocity")),
          dv_indicator_(particles_->template getVariableByName<int>("Indicator")),
          target_velocity(*this) {}

    virtual ~FreeStreamVelocityCorrectionCK() {}

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : transform_(encloser.transform_),
              rho0_(encloser.rho0_),
              physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
              rho_sum_(encloser.dv_rho_sum_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
              target_velocity(encloser.target_velocity) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            if (indicator_[index_i] == 1)
            {
                Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
                Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
                Real frame_u_stream = frame_velocity[0];
                Real u_free = target_velocity(frame_position, frame_velocity,
                                              *physical_time_)[0];
                frame_velocity[0] = u_free + (frame_u_stream - u_free) *
                                                 SMIN(rho_sum_[index_i], rho0_) / rho0_;
                vel_[index_i] = transform_.xformFrameVecToBase(frame_velocity);
            }
        }

      protected:
        Transform transform_;
        Real rho0_;
        Real *physical_time_;
        Real *rho_sum_;
        Vecd *pos_, *vel_;
        int *indicator_;
        TargetVelocity target_velocity;
    };

  protected:
    Transform transform_;
    Real rho0_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Real> *dv_rho_sum_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Vecd> *dv_vel_;
    DiscreteVariable<int> *dv_indicator_;
    TargetVelocity target_velocity;
};

} // namespace fluid_dynamics
} // namespace SPH