/**
 * @file electromagnetic_ophelie_boussinesq.h
 * @brief Stage 4.1 Boussinesq buoyancy force: F = -m β (T − T0) g.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H
#define ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H

#include "force_prior_ck.h"
#include "force_prior_ck.hpp"
#include "external_force.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Boussinesq buoyancy as ForcePrior companion to GravityForceCK.
 * F_b = −mass · β · (T − T_ref) · g_vector
 * with g_vector = Gravity::InducedAcceleration (usually (0,0,−g)).
 */
class BoussinesqBuoyancyForceCK : public LocalDynamics, public ForcePriorCK
{
  public:
    BoussinesqBuoyancyForceCK(SPHBody &sph_body, const Gravity &gravity, Real beta, Real t_ref,
                              const std::string &temperature_field = "Temperature")
        : LocalDynamics(sph_body), ForcePriorCK(this->particles_, "BoussinesqBuoyancyForceCK"), gravity_(gravity),
          beta_(beta), t_ref_(t_ref),
          sv_physical_time_(&sph_system_->svPhysicalTime()),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_mass_(particles_->getVariableByName<Real>("Mass")),
          dv_temperature_(particles_->getVariableByName<Real>(temperature_field))
    {
    }

    class UpdateKernel : public ForcePriorCK::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : ForcePriorCK::UpdateKernel(ex_policy, encloser), gravity_(encloser.gravity_), beta_(encloser.beta_),
              t_ref_(encloser.t_ref_), physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)), mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd g_acc = gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
            this->current_force_[index_i] = -mass_[index_i] * beta_ * (temperature_[index_i] - t_ref_) * g_acc;
            ForcePriorCK::UpdateKernel::update(index_i, dt);
        }

      protected:
        Gravity gravity_;
        Real beta_;
        Real t_ref_;
        Real *physical_time_;
        Vecd *pos_;
        Real *mass_;
        Real *temperature_;
    };

  protected:
    Gravity gravity_;
    Real beta_;
    Real t_ref_;
    SingleVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Real> *dv_temperature_;
};

class OphelieVelocityMaxReduceCK : public BaseLocalDynamicsReduce<ReduceMax, SPHBody>
{
  public:
    explicit OphelieVelocityMaxReduceCK(SPHBody &sph_body)
        : BaseLocalDynamicsReduce<ReduceMax, SPHBody>(sph_body),
          dv_vel_(particles_->template getVariableByName<Vecd>("Velocity"))
    {
        quantity_name_ = "OphelieVelocityMax";
    }

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : vel_(encloser.dv_vel_->DelegatedData(ex_policy))
        {
        }

        Real reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return vel_[index_i].norm();
        }

      protected:
        Vecd *vel_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_vel_;
};

class OphelieTemperatureMeanReduceCK : public BaseLocalDynamicsReduce<ReduceSum<std::pair<Real, Real>>, SPHBody>
{
  public:
    using ReduceReturnType = std::pair<Real, Real>;
    using BaseDynamicsType = BaseLocalDynamicsReduce<ReduceSum<ReduceReturnType>, SPHBody>;

    explicit OphelieTemperatureMeanReduceCK(SPHBody &sph_body,
                                            const std::string &temperature_field = "Temperature")
        : BaseDynamicsType(sph_body),
          dv_temperature_(particles_->template getVariableByName<Real>(temperature_field)),
          dv_vol_(particles_->template getVariableByName<Real>("VolumetricMeasure"))
    {
        quantity_name_ = "OphelieTemperatureMean";
    }

    class FinishDynamics
    {
      public:
        using OutputType = Real;
        template <class EncloserType>
        FinishDynamics(EncloserType &)
        {
        }
        Real Result(const ReduceReturnType &reduced_value)
        {
            return reduced_value.first / (reduced_value.second + TinyReal);
        }
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              vol_(encloser.dv_vol_->DelegatedData(ex_policy))
        {
        }

        ReduceReturnType reduce(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            return ReduceReturnType(temperature_[index_i] * vol_[index_i], vol_[index_i]);
        }

      protected:
        Real *temperature_;
        Real *vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_vol_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_BOUSSINESQ_H
