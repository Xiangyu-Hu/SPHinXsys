#include "viscous_dynamics.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void ViscousAcceleration<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        force += 2.0 * mass_[index_i] * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
    }

    force_prior_[index_i] += force / rho_[index_i];
}
//=================================================================================================//
void ViscousAcceleration<AngularConservative<Inner<>>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
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
        force += eta_ij * mass_[index_i] * inner_neighborhood.dW_ijV_j_[n] * e_ij;
    }

    force_prior_[index_i] += force / rho_[index_i];
}
//=================================================================================================//
void ViscousAcceleration<ContactWall<>>::interaction(size_t index_i, Real dt)
{
    Real rho_i = this->rho_[index_i];
    const Vecd &vel_i = this->vel_[index_i];

    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            Vecd vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
            force += 2.0 * this->mu_ * this->mass_[index_i] * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
        }
    }

    this->force_prior_[index_i] += force;
}
//=================================================================================================//
ViscousAcceleration<Contact<>>::ViscousAcceleration(BaseContactRelation &contact_relation)
    : ViscousAcceleration<FluidContactData>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real mu_k = DynamicCast<Fluid>(this, &contact_particles_[k]->getBaseMaterial())->ReferenceViscosity();
        contact_mu_.push_back(Real(2) * (mu_ * mu_k) / (mu_ + mu_k));
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
    }
}
//=================================================================================================//
void ViscousAcceleration<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real contact_mu_k = this->contact_mu_[k];
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd vel_derivative = (this->vel_[index_i] - vel_k[index_j]) /
                                  (contact_neighborhood.r_ij_[n] + 0.01 * this->smoothing_length_);
            force += 2.0 * this->mass_[index_i] * contact_mu_k * vel_derivative * contact_neighborhood.dW_ijV_j_[n];
        }
    }
    force_prior_[index_i] += force / this->rho_[index_i];
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
