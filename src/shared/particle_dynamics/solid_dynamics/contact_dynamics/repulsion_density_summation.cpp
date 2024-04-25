#include "repulsion_density_summation.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
RepulsionDensitySummation<Inner<>>::
    RepulsionDensitySummation(SelfSurfaceContactRelation &self_contact_relation)
    : RepulsionDensitySummation<Base, SolidDataInner>(self_contact_relation, "SelfRepulsionDensity"),
      mass_(*particles_->getVariableByName<Real>("Mass"))
{
    Real dp_1 = self_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    offset_W_ij_ = self_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(dp_1, ZeroVecd);
}
//=================================================================================================//
void RepulsionDensitySummation<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Real corrected_W_ij = std::max(inner_neighborhood.W_ij_[n] - offset_W_ij_, Real(0));
        sigma += corrected_W_ij * mass_[inner_neighborhood.j_[n]];
    }
    repulsion_density_[index_i] = sigma;
}
//=================================================================================================//
RepulsionDensitySummation<Contact<>>::
    RepulsionDensitySummation(SurfaceContactRelation &solid_body_contact_relation)
    : RepulsionDensitySummation<Base, ContactDynamicsData>(solid_body_contact_relation, "RepulsionDensity"),
      mass_(*particles_->getVariableByName<Real>("Mass")),
      offset_W_ij_(StdVec<Real>(contact_configuration_.size(), 0.0))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(contact_particles_[k]->getVariableByName<Real>("Mass"));
    }

    // we modify the default formulation by an offset, so that exactly touching bodies produce 0 initial force
    // subtract summation of the kernel function of 2 particles at 1 particle distance, and if the result is negative, we take 0
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    // dp1, dp2 half reference spacing
    Real dp_1 = solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real dp_2 = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real distance = 0.5 * dp_1 + 0.5 * dp_2;
        offset_W_ij_[k] = solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(distance, ZeroVecd);
    }
}
//=================================================================================================//
void RepulsionDensitySummation<Contact<>>::interaction(size_t index_i, Real dt)
{
    /** Contact interaction. */
    Real sigma = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], Real(0));
            sigma += corrected_W_ij * contact_mass_k[contact_neighborhood.j_[n]];
        }
    }
    repulsion_density_[index_i] = sigma;
};
//=================================================================================================//
ShellContactDensity::ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation)
    : RepulsionDensitySummation<Base, ContactDynamicsData>(solid_body_contact_relation, "RepulsionDensity"),
      solid_(DynamicCast<Solid>(this, sph_body_.getBaseMaterial())),
      kernel_(solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()),
      particle_spacing_(solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real dp_k = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real average_spacing_k = 0.5 * particle_spacing_ + 0.5 * dp_k;
        Real h_ratio_k = particle_spacing_ / average_spacing_k;
        offset_W_ij_.push_back(kernel_->W(h_ratio_k, average_spacing_k, ZeroVecd));

        Real contact_max(0.0);
        for (int l = 0; l != 3; ++l)
        {
            Real temp = three_gaussian_points_[l] * average_spacing_k * 0.5 + average_spacing_k * 0.5;
            Real contact_temp = 2.0 * (kernel_->W(h_ratio_k, temp, ZeroVecd) - offset_W_ij_[k]) *
                                average_spacing_k * 0.5 * three_gaussian_weights_[l];
            contact_max += Dimensions == 2 ? contact_temp : contact_temp * Pi * temp;
        }
        /** a calibration factor to avoid particle penetration into shell structure */
        calibration_factor_.push_back(solid_.ReferenceDensity() / (contact_max + Eps));

        contact_Vol_.push_back(contact_particles_[k]->getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void ShellContactDensity::interaction(size_t index_i, Real dt)
{
    /** shell contact interaction. */
    Real sigma = 0.0;
    Real contact_density_i = 0.0;

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &contact_Vol_k = *(contact_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], Real(0));
            sigma += corrected_W_ij * contact_Vol_k[contact_neighborhood.j_[n]];
        }
        constexpr Real heuristic_limiter = 0.1;
        // With heuristic_limiter, the maximum contact pressure is heuristic_limiter * K (Bulk modulus).
        // The contact pressure applied to fewer particles than on solids, yielding high acceleration locally,
        // which is one source of instability. Thus, we add a heuristic_limiter
        // to maintain enough contact pressure to prevent penetration while also maintaining stability.
        contact_density_i += heuristic_limiter * sigma * calibration_factor_[k];
    }
    repulsion_density_[index_i] = contact_density_i;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
