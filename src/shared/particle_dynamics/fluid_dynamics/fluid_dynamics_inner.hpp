/**
 * @file 	fluid_dynamics_inner.hpp
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details 	We consider here weakly compressible fluids. The algorithms may be
 * 			different for free surface flow and the one without free surface.
 * @author	Chi Zhang and Xiangyu Hu
 */
#pragma once

#include "fluid_dynamics_inner.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseDensitySummation<DataDelegationType>::BaseDensitySummation(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_),
      rho0_(this->sph_body_.base_material_->ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_body_.sph_adaptation_->LatticeNumberDensity())
{
    this->particles_->registerVariable(rho_sum_, "DensitySummation");
}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseViscousAcceleration<DataDelegationType>::BaseViscousAcceleration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), vel_(this->particles_->vel_),
      acc_prior_(this->particles_->acc_prior_),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseIntegration<DataDelegationType>::BaseIntegration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      fluid_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial())),
      rho_(this->particles_->rho_),
      p_(*this->particles_->template getVariableByName<Real>("Pressure")),
      drho_dt_(*this->particles_->template registerSharedVariable<Real>("DensityChangeRate")),
      pos_(this->particles_->pos_), vel_(this->particles_->vel_),
      acc_(this->particles_->acc_), acc_prior_(this->particles_->acc_prior_) {}
//=================================================================================================//
template <class RiemannSolverType, class PressureType>
BaseIntegration1stHalfInner<RiemannSolverType, PressureType>::
    BaseIntegration1stHalfInner(BaseInnerRelation &inner_relation)
    : BaseIntegration<FluidDataInner>(inner_relation),
      pressure_(fluid_, particles_), riemann_solver_(fluid_, fluid_)
{
    /**
     *	register sortable particle data
     */
    particles_->registerSortableVariable<Vecd>("Position");
    particles_->registerSortableVariable<Vecd>("Velocity");
    particles_->registerSortableVariable<Real>("MassiveMeasure");
    particles_->registerSortableVariable<Real>("Density");
    particles_->registerSortableVariable<Real>("Pressure");
    particles_->registerSortableVariable<Real>("VolumetricMeasure");
    //----------------------------------------------------------------------
    //		add restart output particle data
    //----------------------------------------------------------------------
    particles_->addVariableToRestart<Real>("Pressure");
}
//=================================================================================================//
template <class RiemannSolverType, class PressureType>
void BaseIntegration1stHalfInner<RiemannSolverType, PressureType>::initialization(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    pressure_.update(index_i);
}
//=================================================================================================//
template <class RiemannSolverType, class PressureType>
void BaseIntegration1stHalfInner<RiemannSolverType, PressureType>::update(size_t index_i, Real dt)
{
    vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
}
//=================================================================================================//
template <class RiemannSolverType, class PressureType>
void BaseIntegration1stHalfInner<RiemannSolverType, PressureType>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        acceleration -= 2.0 * pressure_.atInterface(index_i, index_j) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
    }
    acc_[index_i] += acceleration / rho_[index_i];
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration2ndHalfInner<RiemannSolverType>::
    BaseIntegration2ndHalfInner(BaseInnerRelation &inner_relation)
    : BaseIntegration<FluidDataInner>(inner_relation), riemann_solver_(this->fluid_, this->fluid_),
      Vol_(particles_->Vol_), mass_(particles_->mass_) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfInner<RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfInner<RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    Vol_[index_i] = mass_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfInner<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    acc_[index_i] = p_dissipation / rho_[index_i];
};
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//