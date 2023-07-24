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
} // namespace continuum_dynamics
} // namespace SPH