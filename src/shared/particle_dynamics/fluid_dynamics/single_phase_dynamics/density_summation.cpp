#include "density_summation.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
void BaseDensitySummationInner::update(size_t index_i, Real dt)
{
    rho_[index_i] = rho_sum_[index_i];
}
//=================================================================================================//
DensitySummationInner::DensitySummationInner(BaseInnerRelation &inner_relation)
    : BaseDensitySummationInner(inner_relation),
      W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//
void DensitySummationInner::interaction(size_t index_i, Real dt)
{
    Real sigma = W0_;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma += inner_neighborhood.W_ij_[n];

    rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
}
//=================================================================================================//
DensitySummationInnerAdaptive::DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation)
    : BaseDensitySummationInner(inner_relation),
      sph_adaptation_(*sph_body_.sph_adaptation_),
      kernel_(*sph_adaptation_.getKernel()),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void DensitySummationInnerAdaptive::interaction(size_t index_i, Real dt)
{
    Real sigma_i = mass_[index_i] * kernel_.W0(h_ratio_[index_i], ZeroVecd);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma_i += inner_neighborhood.W_ij_[n] * mass_[inner_neighborhood.j_[n]];

    rho_sum_[index_i] = sigma_i * rho0_ * inv_sigma0_ / mass_[index_i] /
                        sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
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
} // namespace fluid_dynamics
} // namespace SPH
