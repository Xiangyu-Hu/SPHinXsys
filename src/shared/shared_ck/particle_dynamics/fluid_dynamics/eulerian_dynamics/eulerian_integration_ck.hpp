#ifndef EULERIAN_INTEGRATION_CK_HPP
#define EULERIAN_INTEGRATION_CK_HPP

#include "eulerian_integration_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class BaseInteractionType>
template <class BaseRelationType>
EulerianIntegrationCK<BaseInteractionType>::EulerianIntegrationCK(BaseRelationType &base_relation)
    : BaseInteractionType(base_relation),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template registerStateVariable<Real>("Pressure")),
      dv_dmass_dt_(this->particles_->template registerStateVariable<Real>("MassChangeRate")),
      dv_vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")),
      dv_mom_(this->particles_->template registerStateVariable<Vecd>("Momentum")),
      dv_dmom_dt_(this->particles_->template registerStateVariable<Vecd>("MomentumChangeRate"))
{
    //----------------------------------------------------------------------
    //		add evolving variables
    //----------------------------------------------------------------------
    this->particles_->template addEvolvingVariable<Real>("Density");
    this->particles_->template addEvolvingVariable<Real>("Pressure");
    this->particles_->template addEvolvingVariable<Vecd>("Velocity");
    //----------------------------------------------------------------------
    //		add output particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToWrite<Real>("Pressure");
    this->particles_->template addVariableToWrite<Vecd>("Velocity");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
template <class BaseRelationType>
EulerianIntegrationCK<Inner<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    EulerianIntegrationCK(BaseRelationType &base_relation)
    : BaseInteraction(base_relation),
      kernel_correction_method_(this->particles_),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, this->sph_body_.getBaseMaterial())),
      riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
EulerianIntegrationCK<Inner<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_method_),
      riemann_solver_(encloser.riemann_solver_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      dmass_dt_(encloser.dv_dmass_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      mom_(encloser.dv_mom_->DelegatedData(ex_policy)),
      dmom_dt_(encloser.dv_dmom_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
void EulerianIntegrationCK<Inner<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real mass_rate = 0.0;
    Vecd mom_rate = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        FluidStateIn state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStateOut state_star = riemann_solver_.InterfaceState(state_i, state_j, e_ij);

        auto correction_ij = correction_(index_i) + correction_(index_j);
        mass_rate -= (correction_ij * state_star.vel_).dot(e_ij) * state_star.rho_ * dW_ijV_j;
        Matd convect_flux = state_star.rho_ * state_star.vel_ * state_star.vel_.transpose();
        mom_rate -= (convect_flux + state_star.p_ * Matd::Identity()) * correction_ij * e_ij * dW_ijV_j;
    }
    dmass_dt_[index_i] += mass_rate * Vol_[index_i];
    dmom_dt_[index_i] += mom_rate * Vol_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
template <class BaseRelationType>
EulerianIntegrationCK<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    EulerianIntegrationCK(BaseRelationType &base_relation)
    : BaseInteraction(base_relation), Interaction<Wall>(base_relation),
      kernel_correction_method_(this->particles_),
      riemann_solver_(DynamicCast<WeaklyCompressibleFluid>(this, this->sph_body_.getBaseMaterial()))
{
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
EulerianIntegrationCK<Contact<Boundary, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_method_),
      riemann_solver_(encloser.riemann_solver_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      dmass_dt_(encloser.dv_dmass_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      mom_(encloser.dv_mom_->DelegatedData(ex_policy)),
      dmom_dt_(encloser.dv_dmom_dt_->DelegatedData(ex_policy)),
      wall_n_(encloser.dv_wall_n_->DelegatedData(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
void EulerianIntegrationCK<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    FluidStateIn state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real mass_rate = 0.0;
    Vecd mom_rate = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        Vecd vel_j_in_wall = 2.0 * wall_vel_ave_[index_j] - state_i.vel_;
        Real p_j_in_wall = state_i.p_;
        Real rho_in_wall = state_i.rho_;

        FluidStateIn state_j(rho_in_wall, vel_j_in_wall, p_j_in_wall);
        FluidStateOut state_star = this->riemann_solver_.InterfaceState(state_i, state_j, wall_n_[index_j]);
        mass_rate -= 2.0 * correction_(index_i) * (state_star.rho_ * state_star.vel_).dot(e_ij) * dW_ijV_j;
        Matd convect_flux = state_star.rho_ * state_star.vel_ * state_star.vel_.transpose();
        mom_rate -= 2.0 * (convect_flux + state_star.p_ * Matd::Identity()) * correction_(index_i) * e_ij * dW_ijV_j;
    }
    dmass_dt_[index_i] += mass_rate * Vol_[index_i];
    dmom_dt_[index_i] += mom_rate * Vol_[index_i];
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
class EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    EulerianIntegrationCK(BaseRelationType &base_relation)
    : BaseDynamicsType(base_relation) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      dmass_dt_(encloser.dv_dmass_dt_->DelegatedData(ex_policy)),
      dmom_dt_(encloser.dv_dmom_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    dmass_dt_[index_i] = 0.0;
    dmom_dt_[index_i] = Vecd::Zero();
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      eos_(encloser.fluid),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      dmass_dt_(encloser.dv_dmass_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      mom_(encloser.dv_mom_->DelegatedData(ex_policy)),
      dmom_dt_(encloser.dv_dmom_dt_->DelegatedData(ex_policy))
      //=================================================================================================//
      template <template <typename...> class RelationType, class... InteractionParameters>
      void EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
          UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];
    p_[index_i] = eos_.getPressure(rho_[index_i]);
    mom_[index_i] += dmom_dt_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_INTEGRATION_CK_HPP