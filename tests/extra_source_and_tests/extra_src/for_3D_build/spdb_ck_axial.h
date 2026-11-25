/**
 * @file spdb_ck_axial.h
 * @brief CK axial support and nominal end force for the SPDB cylinder.
 * @details Same pattern as FixConstraintCK and GravityForceCK. The force scale
 * is a host SingleVariable read through its pointer, so a pseudo-time ramp
 * is visible on later exec calls without rebuilding the cached kernel.
 */
#ifndef SPDB_CK_AXIAL_H
#define SPDB_CK_AXIAL_H

#include "base_local_dynamics.h"
#include "sphinxsys_variable.h"

namespace SPH
{
template <class DynamicsIdentifier>
class AxialZSupportCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    explicit AxialZSupportCK(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
          dv_pos0_(this->particles_->template registerStateVariableFrom<Vecd>("InitialPosition", "Position")),
          dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity"))
    {
        this->particles_->template addEvolvingVariable<Vecd>(dv_pos0_);
    }
    virtual ~AxialZSupportCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            pos_[index_i][2] = pos0_[index_i][2];
            vel_[index_i][2] = 0.0;
        }

      protected:
        Vecd *pos_, *pos0_, *vel_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos0_, *dv_vel_;
};

/** Writes ForcePrior on a body part. z component is force_per_particle * scale. */
template <class DynamicsIdentifier>
class NominalAxialForceCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    NominalAxialForceCK(DynamicsIdentifier &identifier, Real force_per_particle, SingleVariable<Real> &scale)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_force_prior_(this->particles_->template getVariableByName<Vecd>("ForcePrior")),
          sv_scale_(&scale), force_per_particle_(force_per_particle) {}
    virtual ~NominalAxialForceCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
              scale_(encloser.sv_scale_->DelegatedData(ex_policy)),
              force_per_particle_(encloser.force_per_particle_) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            force_prior_[index_i][0] = 0.0;
            force_prior_[index_i][1] = 0.0;
            force_prior_[index_i][2] = force_per_particle_ * (*scale_);
        }

      protected:
        Vecd *force_prior_;
        Real *scale_;
        Real force_per_particle_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_force_prior_;
    SingleVariable<Real> *sv_scale_;
    Real force_per_particle_;
};
} // namespace SPH

#endif
