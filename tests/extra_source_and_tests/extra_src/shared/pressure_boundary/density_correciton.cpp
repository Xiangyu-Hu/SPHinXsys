#include "density_correciton.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void DensitySummationPressure<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = W0_;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma += inner_neighborhood.W_ij_[n];

    rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
}
//=================================================================================================//
DensitySummationPressure<Contact<Base>>::DensitySummationPressure(BaseContactRelation &contact_relation)
    : DensitySummationPressure<Base, DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
    }
}
//=================================================================================================//
Real DensitySummationPressure<Contact<Base>>::ContactSummation(size_t index_i)
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
void DensitySummationPressure<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = DensitySummationPressure<Contact<Base>>::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
}
} // namespace fluid_dynamics
} // namespace SPH
