#include "density_summation.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void DensitySummation<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = W0_;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma += inner_neighborhood.W_ij_[n];

    rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
}
//=================================================================================================//
void DensitySummation<Inner<>>::update(size_t index_i, Real dt)
{
    rho_[index_i] = rho_sum_[index_i];
    Vol_[index_i] = mass_[index_i] / rho_[index_i];
}
//=================================================================================================//
DensitySummation<Inner<AdaptiveSmoothingLength>>::DensitySummation(BaseInnerRelation &inner_relation)
    : DensitySummation<Inner<Base>>(inner_relation),
      sph_adaptation_(getSPHAdaptation()),
      kernel_(*sph_adaptation_.getKernel()),
      h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void DensitySummation<Inner<AdaptiveSmoothingLength>>::update(size_t index_i, Real dt)
{
    rho_[index_i] = rho_sum_[index_i];
    Vol_[index_i] = mass_[index_i] / rho_[index_i];
}
//=================================================================================================//
void DensitySummation<Inner<AdaptiveSmoothingLength>>::interaction(size_t index_i, Real dt)
{
    Real sigma_i = mass_[index_i] * kernel_.W0(h_ratio_[index_i], ZeroVecd);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma_i += inner_neighborhood.W_ij_[n] * mass_[inner_neighborhood.j_[n]];

    rho_sum_[index_i] = sigma_i * rho0_ * inv_sigma0_ / mass_[index_i] /
                        sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
//=================================================================================================//
DensitySummation<Contact<Base>>::DensitySummation(BaseContactRelation &contact_relation)
    : DensitySummation<Base, DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
    }
}
//=================================================================================================//
Real DensitySummation<Contact<Base>>::ContactSummation(size_t index_i)
{
    Real sigma(0.0);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *contact_mass_k = this->contact_mass_[k];
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
void DensitySummation<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = DensitySummation<Contact<Base>>::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
}
//=================================================================================================//
DensitySummation<Contact<AdaptiveSmoothingLength>>::
    DensitySummation(BaseContactRelation &contact_relation)
    : DensitySummation<Contact<Base>>(contact_relation),
      sph_adaptation_(getSPHAdaptation()),
      h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void DensitySummation<Contact<AdaptiveSmoothingLength>>::interaction(size_t index_i, Real dt)
{
    Real sigma = DensitySummation<Contact<Base>>::ContactSummation(index_i);
    rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i] /
                         sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
