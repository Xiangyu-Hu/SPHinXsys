#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BaseDensitySummationContact::BaseDensitySummationContact(BaseContactRelation &contact_relation)
    : BaseDensitySummation<FluidContactData>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        contact_mass_.push_back(&(this->contact_particles_[k]->mass_));
    }
}
//=================================================================================================//
Real BaseDensitySummationContact::ContactSummation(size_t index_i)
{
    Real sigma(0.0);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &contact_mass_k = *(this->contact_mass_[k]);
        Real contact_inv_rho0_k = contact_inv_rho0_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            sigma += contact_neighborhood.W_ij_[n] * contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
        }
    }
    return sigma;
};
//=================================================================================================//
void DensitySummationContact::interaction(size_t index_i, Real dt)
{
    Real sigma = BaseDensitySummationContact::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
}
//=================================================================================================//
DensitySummationContactAdaptive::
    DensitySummationContactAdaptive(BaseContactRelation &contact_relation)
    : BaseDensitySummationContact(contact_relation),
      sph_adaptation_(*sph_body_.sph_adaptation_),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void DensitySummationContactAdaptive::interaction(size_t index_i, Real dt)
{
    Real sigma = BaseDensitySummationContact::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i] /
                         sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
//=================================================================================================//
void ViscousAccelerationWithWall::interaction(size_t index_i, Real dt)
{
    Real rho_i = this->rho_[index_i];
    const Vecd &vel_i = this->vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
            acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
        }
    }

    this->acc_prior_[index_i] += acceleration;
}
//=================================================================================================//
Oldroyd_BIntegration1stHalfWithWall::
    Oldroyd_BIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation)
    : Integration1stHalfWithWallDissipativeRiemann(wall_contact_relation),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")){};
//=================================================================================================//
void Oldroyd_BIntegration1stHalfWithWall::interaction(size_t index_i, Real dt)
{
    Integration1stHalfWithWallDissipativeRiemann::interaction(index_i, dt);

    Real rho_i = rho_[index_i];
    Matd tau_i = tau_[index_i];

    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            /** stress boundary condition. */
            acceleration += 2.0 * tau_i * nablaW_ijV_j / rho_i;
        }
    }

    acc_[index_i] += acceleration;
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalfWithWall::
    Oldroyd_BIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation)
    : Integration2ndHalfWithWallDissipativeRiemann(wall_contact_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalfWithWall::interaction(size_t index_i, Real dt)
{
    Integration2ndHalfWithWallDissipativeRiemann::interaction(index_i, dt);

    Vecd vel_i = vel_[index_i];
    Matd tau_i = tau_[index_i];

    Matd stress_rate = Matd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];

            Matd velocity_gradient = -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
            stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
                           (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
        }
    }
    dtau_dt_[index_i] += stress_rate;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
