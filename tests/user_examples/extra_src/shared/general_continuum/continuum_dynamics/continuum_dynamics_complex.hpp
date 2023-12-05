#pragma once
#include "continuum_dynamics_complex.h"
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
namespace continuum_dynamics
{
//=================================================================================================//
template <class BaseIntegrationType>
template <class BaseBodyRelationType, typename... Args>
InteractionWithWall<BaseIntegrationType>::
    InteractionWithWall(BaseBodyRelationType &base_body_relation,
                        BaseContactRelation &wall_contact_relation, Args &&...args)
    : BaseIntegrationType(base_body_relation, std::forward<Args>(args)...),
      FSIContactData(wall_contact_relation)
{
    if (&base_body_relation.getSPHBody() != &wall_contact_relation.getSPHBody())
    {
        std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    for (size_t k = 0; k != FSIContactData::contact_particles_.size(); ++k)
    {
        Real rho0_k = FSIContactData::contact_bodies_[k]->base_material_->ReferenceDensity();
        wall_inv_rho0_.push_back(1.0 / rho0_k);
        wall_mass_.push_back(&(FSIContactData::contact_particles_[k]->mass_));
        wall_vel_ave_.push_back(FSIContactData::contact_particles_[k]->AverageVelocity());
        wall_force_ave_.push_back(FSIContactData::contact_particles_[k]->AverageForce());
        wall_n_.push_back(&(FSIContactData::contact_particles_[k]->n_));
    }
}
//=================================================================================================//
template <class BaseStressRelaxation1stHalfType>
Vecd BaseStressRelaxation1stHalfWithWall<BaseStressRelaxation1stHalfType>::computeNonConservativeForce(size_t index_i)
{
    return this->force_prior_[index_i];
}

template <class BaseStressRelaxation1stHalfType>
void BaseStressRelaxation1stHalfWithWall<BaseStressRelaxation1stHalfType>::interaction(size_t index_i, Real dt)
{
    BaseStressRelaxation1stHalfType::interaction(index_i, dt);

    Vecd force_prior_i = computeNonConservativeForce(index_i);
    Vecd force = force_prior_i;
    Real rho_dissipation(0);
    Matd stress_tensor_i = this->reduceTensor(this->stress_tensor_3D_[index_i]);
    for (size_t k = 0; k < FSIContactData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd>& force_ave_k = *(this->wall_force_ave_[k]);
        StdLargeVec<Real>& wall_mass_k = *(this->wall_mass_[k]);
        Neighborhood& wall_neighborhood = (*FSIContactData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd& e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Real face_wall_external_acceleration = (force_prior_i / this->mass_[index_i] - force_ave_k[index_j] / wall_mass_k[index_j]).dot(-e_ij);
            Real p_in_wall = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            force += 2 * this->mass_[index_i] * stress_tensor_i * wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
            rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    this->force_[index_i] += force / this->rho_[index_i];
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
    for (size_t k = 0; k < FSIContactData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
        StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
        Neighborhood& wall_neighborhood = (*FSIContactData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd& e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            density_change_rate += (this->vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (this->vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += this->mass_[index_i] * this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
            velocity_gradient -= (vel_i - vel_in_wall) * dW_ijV_j * e_ij.transpose();
        }
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->force_[index_i] += p_dissipation / this->rho_[index_i];
    this->velocity_gradient_[index_i] += velocity_gradient;
}
} // namespace continuum_dynamics
} // namespace SPH
