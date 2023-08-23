#pragma once
#include "continuum_dynamics_complex.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
namespace continuum_dynamics
{
//=================================================================================================//
//==========================BaseShearStressRelaxation1stHalfWithWall================================//
//=================================================================================================//
template <class BaseShearStressRelaxation1stHalfType>
void BaseShearStressRelaxation1stHalfWithWall<BaseShearStressRelaxation1stHalfType>::interaction(size_t index_i, Real dt)
{
    BaseShearStressRelaxation1stHalfType::interaction(index_i, dt);
    Matd shear_stress_i = this->shear_stress_[index_i];
    Real rho_i = this->rho_[index_i];
    Real rho_in_wall = this->continuum_.getDensity();
    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k) // There may be several wall bodies.
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        StdLargeVec<Vecd> &acc_ave_k = *(this->wall_acc_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];

            Matd shear_stress_in_wall = shear_stress_i;
            acceleration += this->rho_[index_j] * (shear_stress_i / (rho_i * rho_i) + shear_stress_in_wall / (rho_in_wall * rho_in_wall)) * nablaW_ijV_j;
        }
    }
    this->acc_shear_[index_i] += acceleration;
}
//=================================================================================================//
//==========================BaseShearStressRelaxation2ndHalfWithWall================================//
//=================================================================================================//
template <class BaseShearStressRelaxation2ndHalfType>
void BaseShearStressRelaxation2ndHalfWithWall<BaseShearStressRelaxation2ndHalfType>::interaction(size_t index_i, Real dt)
{
    BaseShearStressRelaxation2ndHalfType::interaction(index_i, dt);

    Vecd vel_i = this->vel_[index_i];
    Matd velocity_gradient = Matd::Zero();
    Real rho_in_wall = this->continuum_.getDensity();
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {

        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_i;

            Matd velocity_gradient_ij = -(vel_i - vel_in_wall) * nablaW_ijV_j.transpose();
            velocity_gradient += velocity_gradient_ij;
        }
    }
    this->velocity_gradient_[index_i] += velocity_gradient;
}

template <class BaseStressRelaxation1stHalfType>
Vecd BaseStressRelaxation1stHalfWithWall<BaseStressRelaxation1stHalfType>::computeNonConservativeAcceleration(size_t index_i)
{
    return this->acc_prior_[index_i];
}

template <class BaseStressRelaxation1stHalfType>
void BaseStressRelaxation1stHalfWithWall<BaseStressRelaxation1stHalfType>::interaction(size_t index_i, Real dt)
{
    BaseStressRelaxation1stHalfType::interaction(index_i, dt);

    Vecd acc_prior_i = computeNonConservativeAcceleration(index_i);
    Vecd acceleration = acc_prior_i;
    Real rho_dissipation(0);

    Matd stress_tensor_i = this->reduceTensor(this->stress_tensor_3D_[index_i]);

    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd>& acc_ave_k = *(this->wall_acc_ave_[k]);
        Neighborhood& wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd& e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];

            Real face_wall_external_acceleration = (acc_prior_i - acc_ave_k[index_j]).dot(-e_ij);
            Real p_in_wall = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            acceleration += 2 * stress_tensor_i * wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    this->acc_[index_i] += acceleration / this->rho_[index_i];
    this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class BaseStressRelaxation2ndHalfType>
void BaseStressRelaxation2ndHalfWithWall<BaseStressRelaxation2ndHalfType>::interaction(size_t index_i, Real dt)
{
    BaseStressRelaxation2ndHalfType::interaction(index_i, dt);

    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();

    Vecd vel_i = this->vel_[index_i];
    Matd velocity_gradient = Matd::Zero();

    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
        StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);

        Neighborhood& wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd& e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            Vecd vel_in_wall = 1.9 * vel_ave_k[index_j] - 0.9 *vel_i;
            Matd velocity_gradient_ij = -(vel_i - vel_in_wall) * nablaW_ijV_j.transpose();
            velocity_gradient += velocity_gradient_ij;

            density_change_rate += (this->vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Vecd u_jump = 2.0 * (this->vel_[index_i] - vel_ave_k[index_j]);
            p_dissipation += this->riemann_solver_.DissipativePJumpExtra(u_jump, n_k[index_j]) * dW_ijV_j;
        }
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->velocity_gradient_[index_i] += velocity_gradient;
    this->acc_[index_i] += p_dissipation / this->rho_[index_i];
}

template <class BaseStressDiffusionType>
void BaseStressDiffusionWithWall<BaseStressDiffusionType>::interaction(size_t index_i, Real dt)
{
    BaseStressDiffusionType::interaction(index_i, dt);

    Real rho_in_wall = this->plastic_continuum_.getDensity();
    Real gravity = abs(this->acc_prior_[index_i](1, 0));
    Real density = rho_in_wall;
    Mat3d diffusion_stress_rate_ = Mat3d::Zero();
    Mat3d diffusion_stress_ = Mat3d::Zero();
    Mat3d stress_tensor_i = this->stress_tensor_3D_[index_i];
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        Neighborhood& wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Real y_ij = this->pos_[index_i](1, 0) - this->pos_[index_j](1, 0);
            //stress boundary condition
            Mat3d stress_tensor_j = stress_tensor_i;
            diffusion_stress_ = stress_tensor_i - stress_tensor_j;
            diffusion_stress_(0, 0) = diffusion_stress_(0, 0) - (1 - sin(this->fai_)) * density * gravity * y_ij;
            diffusion_stress_(1, 1) = diffusion_stress_(1, 1) - density * gravity * y_ij;
            diffusion_stress_(2, 2) = diffusion_stress_(2, 2) - (1 - sin(this->fai_)) * density * gravity * y_ij;
            diffusion_stress_rate_ += 2 * this->zeta_ * this->smoothing_length_ * this->sound_speed_ * diffusion_stress_ * r_ij * dW_ijV_j / (r_ij * r_ij + 0.01 * this->smoothing_length_);
        }
    }
    this->stress_rate_3D_[index_i] += diffusion_stress_rate_;
}

} // namespace continuum_dynamics
} // namespace SPH