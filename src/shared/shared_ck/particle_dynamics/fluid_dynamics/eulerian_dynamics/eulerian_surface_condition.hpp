#ifndef EULERIAN_SURFACE_CONDITION_HPP
#define EULERIAN_SURFACE_CONDITION_HPP

#include "eulerian_surface_condition.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
template <typename... Args>
EulerianSurfaceCondition<KernelCorrectionType, ConditionType>::
    EulerianSurfaceCondition(BodyPartByParticle &box_part, Args &&...args)
    : BaseLocalDynamics<BodyPartByParticle>(box_part),
      BaseStateCondition(this->particles_),
      kernel_correction_method_(this->particles_),
      condition_(std::forward<Args>(args)...),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      dv_n_(nullptr), dv_surface_area_(nullptr)
{
    Vecd *kernel_gradient_integral =
        this->particles_->template getVariableDataByName<Vecd>("KernelGradientIntegral");
    Real *Vol = this->particles_->template getVariableDataByName<Real>("VolumetricMeasure");
    dv_n_ = this->particles_->template registerStateVariable<Vecd>(
        "NormalDirection",
        [&](size_t index)
        {
            return (-kernel_gradient_integral[index] * Vol[index]).normalized();
        });
    dv_surface_area_ = this->particles_->template registerStateVariable<Real>(
        "SurfaceArea",
        [&](size_t index)
        {
            return (kernel_gradient_integral[index].norm() * Vol[index]).norm();
        });
}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
template <class ExecutionPolicy, class EncloserType>
EulerianSurfaceCondition<KernelCorrectionType, ConditionType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseStateCondition::ComputingKernel(ex_policy, encloser),
      correction_kernel_(ex_policy, encloser.kernel_correction_method_),
      eos_(encloser.fluid_), condition_(encloser.condition_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      dmom_dt_(encloser.dv_dmom_dt_->DelegatedData(ex_policy)),
      surface_area_(encloser.dv_surface_area_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      dmass_dt_(encloser.dv_dmass_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
void EulerianSurfaceCondition<KernelCorrectionType, ConditionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real p = condition_.getPressure(p_[index_i], *physical_time_);
    Vecd vel = condition_.getVelocity(vel_[index_i], *physical_time_);
    Real rho = eos_.DensityFromPressure(p);
    FluidStateIn state_boundary = State(rho, vel, p);
    FluidStateOut state_star = this->riemann_solver_.InterfaceState(state_i, state_boundary, -n_[index_i]);
    dmass_dt_[index_i] -= 2.0 * surface_area[index_i] * (state_star.rho_ * state_star.vel_).dot(n_[index_i]);
    Matd convect_flux = state_star.rho_ * state_star.vel_ * state_star.vel_.transpose();
    dmom_dt_[index_i] -= 2.0 * surface_area[index_i] * (convect_flux + state_star.p_ * Matd::Identity()) * n_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_SURFACE_CONDITION_HPP
