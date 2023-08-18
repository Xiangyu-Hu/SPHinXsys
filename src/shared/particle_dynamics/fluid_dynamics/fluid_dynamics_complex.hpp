/**
 * @file 	fluid_dynamics_complex.hpp
 * @brief 	Here, we define the algorithm classes for complex fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Chi Zhang and Xiangyu Hu
 */
#pragma once

#include "fluid_dynamics_complex.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
void DensitySummationComplex::
    interaction(size_t index_i, Real dt)
{
    BaseDensitySummationComplex<DensitySummationInner>::interaction(index_i, dt);
    Real sigma = BaseDensitySummationComplex<DensitySummationInner>::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
}
//=================================================================================================//
void DensitySummationComplexAdaptive::
    interaction(size_t index_i, Real dt)
{
    BaseDensitySummationComplex<DensitySummationInnerAdaptive>::interaction(index_i, dt);
    Real sigma = BaseDensitySummationComplex<DensitySummationInnerAdaptive>::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i] /
                         sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
//=================================================================================================//
void TransportVelocityCorrectionComplex::
    interaction(size_t index_i, Real dt)
{
    TransportVelocityCorrectionInner::interaction(index_i, dt);

    Vecd acceleration_trans = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

            // acceleration for transport velocity
            acceleration_trans -= 2.0 * nablaW_ijV_j;
        }
    }

    /** correcting particle position */
    if (surface_indicator_[index_i] == 0)
        pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
}
//=================================================================================================//
void TransportVelocityCorrectionComplexAdaptive::
    interaction(size_t index_i, Real dt)
{
    TransportVelocityCorrectionInnerAdaptive::interaction(index_i, dt);

    Vecd acceleration_trans = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

            // acceleration for transport velocity
            acceleration_trans -= 2.0 * nablaW_ijV_j;
        }
    }

    /** correcting particle position */
    if (surface_indicator_[index_i] == 0)
    {
        Real inv_h_ratio = 1.0 / sph_adaptation_.SmoothingLengthRatio(index_i);
        pos_[index_i] += coefficient_ * smoothing_length_sqr_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalfWithWall::
    interaction(size_t index_i, Real dt)
{
    BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>::interaction(index_i, dt);

    Real rho_i = rho_[index_i];
    Matd tau_i = tau_[index_i];

    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
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
void Oldroyd_BIntegration2ndHalfWithWall::
    interaction(size_t index_i, Real dt)
{
    BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>::interaction(index_i, dt);

    Vecd vel_i = vel_[index_i];
    Matd tau_i = tau_[index_i];

    Matd stress_rate = Matd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
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
template <class BaseIntegrationType>
template <class BaseBodyRelationType, typename... Args>
InteractionWithWall<BaseIntegrationType>::
    InteractionWithWall(BaseContactRelation &wall_contact_relation,
                        BaseBodyRelationType &base_body_relation, Args &&...args)
    : BaseIntegrationType(base_body_relation, std::forward<Args>(args)...),
      FluidWallData(wall_contact_relation)
{
    if (&base_body_relation.getSPHBody() != &wall_contact_relation.getSPHBody())
    {
        std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    for (size_t k = 0; k != FluidWallData::contact_particles_.size(); ++k)
    {
        Real rho0_k = FluidWallData::contact_bodies_[k]->base_material_->ReferenceDensity();
        wall_inv_rho0_.push_back(1.0 / rho0_k);
        wall_mass_.push_back(&(FluidWallData::contact_particles_[k]->mass_));
        wall_vel_ave_.push_back(FluidWallData::contact_particles_[k]->AverageVelocity());
        wall_acc_ave_.push_back(FluidWallData::contact_particles_[k]->AverageAcceleration());
        wall_n_.push_back(&(FluidWallData::contact_particles_[k]->n_));
    }
}
//=================================================================================================//
template <class DensitySummationInnerType>
template <typename... Args>
BaseDensitySummationComplex<DensitySummationInnerType>::
    BaseDensitySummationComplex(Args &&...args)
    : BaseInteractionComplex<DensitySummationInnerType, FluidContactData>(std::forward<Args>(args)...)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        contact_mass_.push_back(&(this->contact_particles_[k]->mass_));
    }
};
//=================================================================================================//
template <class DensitySummationInnerType>
Real BaseDensitySummationComplex<DensitySummationInnerType>::ContactSummation(size_t index_i)
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
template <class ViscousAccelerationInnerType>
void BaseViscousAccelerationWithWall<ViscousAccelerationInnerType>::
    interaction(size_t index_i, Real dt)
{
    ViscousAccelerationInnerType::interaction(index_i, dt);

    Real rho_i = this->rho_[index_i];
    const Vecd &vel_i = this->vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        Neighborhood &contact_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
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
template <class BaseIntegration1stHalfType>
void BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::
    interaction(size_t index_i, Real dt)
{
    BaseIntegration1stHalfType::interaction(index_i, dt);

    Vecd acc_prior_i = computeNonConservativeAcceleration(index_i);

    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &acc_ave_k = *(this->wall_acc_ave_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];

            Real face_wall_external_acceleration = (acc_prior_i - acc_ave_k[index_j]).dot(-e_ij);
            Real p_in_wall = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            acceleration -= (this->p_[index_i] + p_in_wall) * dW_ijV_j * e_ij;
            rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    this->acc_[index_i] += acceleration / this->rho_[index_i];
    this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class BaseIntegration1stHalfType>
Vecd BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::computeNonConservativeAcceleration(size_t index_i)
{
    return this->acc_prior_[index_i];
}
//=================================================================================================//
template <class BaseIntegration1stHalfType>
void BaseExtendIntegration1stHalfWithWall<BaseIntegration1stHalfType>::initialization(size_t index_i, Real dt)
{
    BaseIntegration1stHalfType::initialization(index_i, dt);
    non_cnsrv_acc_[index_i] = Vecd::Zero();
}
//=================================================================================================//
template <class BaseIntegration1stHalfType>
void BaseExtendIntegration1stHalfWithWall<BaseIntegration1stHalfType>::
    interaction(size_t index_i, Real dt)
{
    BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::interaction(index_i, dt);

    Real rho_i = this->rho_[index_i];
    Real penalty_pressure = this->p_[index_i];
    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        Real particle_spacing_j1 = 1.0 / FluidWallData::contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real particle_spacing_ratio2 = 1.0 / (this->sph_body_.sph_adaptation_->ReferenceSpacing() * particle_spacing_j1);
        particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
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
            acceleration -= 2.0 * penalty * n_j * dW_ijV_j / rho_i;
        }
    }
    this->acc_[index_i] += acceleration;
}
//=================================================================================================//
template <class BaseIntegration1stHalfType>
Vecd BaseExtendIntegration1stHalfWithWall<BaseIntegration1stHalfType>::
    computeNonConservativeAcceleration(size_t index_i)
{
    Vecd acceleration = BaseIntegration1stHalfType::computeNonConservativeAcceleration(index_i);
    non_cnsrv_acc_[index_i] = acceleration;
    return acceleration;
}
//=================================================================================================//
template <class BaseIntegration2ndHalfType>
void BaseIntegration2ndHalfWithWall<BaseIntegration2ndHalfType>::
    interaction(size_t index_i, Real dt)
{
    BaseIntegration2ndHalfType::interaction(index_i, dt);

    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            density_change_rate += (this->vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (this->vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
        }
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->acc_[index_i] += p_dissipation / this->rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//