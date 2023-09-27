#include "viscous_dynamics.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
    }

    acc_prior_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
void AngularConservativeViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];

        /** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
         * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
        Real v_r_ij = (vel_[index_i] - vel_[index_j]).dot(r_ij * e_ij);
        Real eta_ij = 8.0 * mu_ * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
        acceleration += eta_ij * inner_neighborhood.dW_ijV_j_[n] * e_ij;
    }

    acc_prior_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
void ViscousWallBoundary::interaction(size_t index_i, Real dt)
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
VorticityInner::VorticityInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      vel_(particles_->vel_)
{
    particles_->registerVariable(vorticity_, "VorticityInner");
    particles_->addVariableToWrite<AngularVecd>("VorticityInner");
}
//=================================================================================================//
void VorticityInner::interaction(size_t index_i, Real dt)
{
    AngularVecd vorticity = ZeroData<AngularVecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Vecd vel_diff = vel_[index_i] - vel_[index_j];
        vorticity += getCrossProduct(vel_diff, inner_neighborhood.e_ij_[n]) * inner_neighborhood.dW_ijV_j_[n];
    }

    vorticity_[index_i] = vorticity;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
