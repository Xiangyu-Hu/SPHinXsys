#pragma once

#include "fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseIntegration<DataDelegationType>::BaseIntegration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      fluid_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial())),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_),
      p_(*this->particles_->template getVariableByName<Real>("Pressure")),
      drho_dt_(*this->particles_->template registerSharedVariable<Real>("DensityChangeRate")),
      pos_(this->particles_->pos_), vel_(this->particles_->vel_),
      force_(this->particles_->force_), force_prior_(this->particles_->force_prior_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>::
    Integration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration<FluidDataInner>(inner_relation),
      correction_(particles_), riemann_solver_(fluid_, fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
    //----------------------------------------------------------------------
    //		register sortable particle data
    //----------------------------------------------------------------------
    particles_->registerSortableVariable<Vecd>("Position");
    particles_->registerSortableVariable<Vecd>("Velocity");
    particles_->registerSortableVariable<Real>("Mass");
    particles_->registerSortableVariable<Vecd>("ForcePrior");
    particles_->registerSortableVariable<Vecd>("Force");
    particles_->registerSortableVariable<Real>("DensityChangeRate");
    particles_->registerSortableVariable<Real>("Density");
    particles_->registerSortableVariable<Real>("Pressure");
    particles_->registerSortableVariable<Real>("VolumetricMeasure");
    //----------------------------------------------------------------------
    //		add restart output particle data
    //----------------------------------------------------------------------
    particles_->addVariableToRestart<Real>("Pressure");
    particles_->addVariableToRestart<Real>("DensityChangeRate");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>::initialization(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        force -= mass_[index_i] * (p_[index_i] * correction_(index_i) + p_[index_j] * correction_(index_j)) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
    }
    force_[index_i] += force / rho_[index_i];
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>::
    Integration1stHalf(BaseContactRelation &wall_contact_relation)
    : BaseIntegrationWithWall(wall_contact_relation),
      correction_(particles_), riemann_solver_(fluid_, fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &force_ave_k = *(wall_force_ave_[k]);
        StdLargeVec<Real> &wall_mass_k = *(wall_mass_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];

            Real face_wall_external_acceleration = (force_prior_[index_i] / mass_[index_i] - force_ave_k[index_j] / wall_mass_k[index_j]).dot(-e_ij);
            Real p_in_wall = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            force -= mass_[index_i] * (p_[index_i] + p_in_wall) * correction_(index_i) * dW_ijV_j * e_ij;
            rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    force_[index_i] += force / rho_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
Integration1stHalf<Contact<Wall, Extended>, RiemannSolverType, KernelCorrectionType>::
    Integration1stHalf(BaseContactRelation &wall_contact_relation, Real penalty_strength)
    : Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>(wall_contact_relation),
      penalty_strength_(penalty_strength) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Contact<Wall, Extended>, RiemannSolverType, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>::interaction(index_i, dt);

    Real rho_i = this->rho_[index_i];
    Real penalty_pressure = this->p_[index_i];
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real particle_spacing_j1 = 1.0 / this->contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real particle_spacing_ratio2 = 1.0 / (this->sph_body_.sph_adaptation_->ReferenceSpacing() * particle_spacing_j1);
        particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Vecd &n_j = n_k[index_j];

            /** penalty method to prevent particle running into boundary */
            Real projection = e_ij.dot(n_j);
            Real delta = 2.0 * projection * r_ij * particle_spacing_j1;
            Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
            // penalty must be positive so that the penalty force is pointed to fluid inner domain
            Real penalty = penalty_strength_ * beta * fabs(projection * penalty_pressure);

            // penalty force induced acceleration
            force -= 2.0 * this->mass_[index_i] * penalty * n_j * dW_ijV_j / rho_i;
        }
    }
    this->force_[index_i] += force;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
Integration1stHalf<Contact<>, RiemannSolverType, KernelCorrectionType>::
    Integration1stHalf(BaseContactRelation &contact_relation)
    : BaseIntegration<FluidContactData>(contact_relation),
      correction_(this->particles_)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_corrections_.push_back(KernelCorrectionType(this->contact_particles_[k]));
        Fluid &contact_fluid = DynamicCast<Fluid>(this, this->contact_particles_[k]->getBaseMaterial());
        riemann_solvers_.push_back(RiemannSolverType(this->fluid_, contact_fluid));
        contact_p_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Pressure"));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType>
void Integration1stHalf<Contact<>, RiemannSolverType, KernelCorrectionType>::
    interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &p_k = *(this->contact_p_[k]);
        KernelCorrectionType &correction_k = contact_corrections_[k];
        RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

            force -= this->mass_[index_i] * riemann_solver_k.AverageP(this->p_[index_i] * correction_k(index_i), p_k[index_j] * correction_(index_j)) *
                     2.0 * e_ij * dW_ijV_j;
            rho_dissipation += riemann_solver_k.DissipativeUJump(this->p_[index_i] - p_k[index_j]) * dW_ijV_j;
        }
    }
    this->force_[index_i] += force / this->rho_[index_i];
    this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
Integration2ndHalf<Inner<>, RiemannSolverType>::
    Integration2ndHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration<FluidDataInner>(inner_relation), riemann_solver_(this->fluid_, this->fluid_),
      Vol_(particles_->Vol_), mass_(particles_->mass_) {}
//=================================================================================================//
template <class RiemannSolverType>
void Integration2ndHalf<Inner<>, RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType>
void Integration2ndHalf<Inner<>, RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    Vol_[index_i] = mass_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void Integration2ndHalf<Inner<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
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
        p_dissipation += mass_[index_i] * riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation / rho_[index_i];
};
//=================================================================================================//
template <class RiemannSolverType>
Integration2ndHalf<Contact<Wall>, RiemannSolverType>::
    Integration2ndHalf(BaseContactRelation &wall_contact_relation)
    : BaseIntegrationWithWall(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType>
void Integration2ndHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_[index_i];
            density_change_rate += (vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += mass_[index_i] * riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
        }
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
Integration2ndHalf<Contact<>, RiemannSolverType>::
    Integration2ndHalf(BaseContactRelation &contact_relation)
    : BaseIntegration<FluidContactData>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Fluid &contact_fluid = DynamicCast<Fluid>(this, contact_particles_[k]->getBaseMaterial());
        riemann_solvers_.push_back(RiemannSolverType(fluid_, contact_fluid));
        contact_vel_.push_back(contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void Integration2ndHalf<Contact<>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        RiemannSolverType &riemann_solver_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ijV_j_[n];

            Vecd vel_ave = riemann_solver_k.AverageV(this->vel_[index_i], vel_k[index_j]);
            density_change_rate += 2.0 * (this->vel_[index_i] - vel_ave).dot(e_ij) * dW_ijV_j;
            Real u_jump = (this->vel_[index_i] - vel_k[index_j]).dot(e_ij);
            p_dissipation += this->mass_[index_i] * riemann_solver_k.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
        }
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->force_[index_i] += p_dissipation / this->rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
