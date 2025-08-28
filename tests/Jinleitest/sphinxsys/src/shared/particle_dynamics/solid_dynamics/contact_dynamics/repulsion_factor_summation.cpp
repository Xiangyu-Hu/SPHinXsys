#include "repulsion_factor_summation.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
RepulsionFactorSummation<Inner<>>::
    RepulsionFactorSummation(SelfSurfaceContactRelation &self_contact_relation)
    : RepulsionFactorSummation<Base, DataDelegateInner>(self_contact_relation, "SelfRepulsionFactor")
{
    Real dp_1 = self_contact_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing();
    offset_W_ij_ = self_contact_relation.getSPHBody().getSPHAdaptation().getKernel()->W(dp_1, ZeroVecd);
}
//=================================================================================================//
void RepulsionFactorSummation<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real sigma = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Real corrected_W_ij = std::max(inner_neighborhood.W_ij_[n] - offset_W_ij_, Real(0));
        sigma += corrected_W_ij * particles_->ParticleVolume(inner_neighborhood.j_[n]);
    }
    repulsion_factor_[index_i] = sigma;
}
//=================================================================================================//
RepulsionFactorSummation<Contact<>>::
    RepulsionFactorSummation(SurfaceContactRelation &solid_body_contact_relation)
    : RepulsionFactorSummation<Base, DataDelegateContact>(solid_body_contact_relation, "RepulsionFactor") {}
//=================================================================================================//
void RepulsionFactorSummation<Contact<>>::interaction(size_t index_i, Real dt)
{
    /** Contact interaction. */
    Real sigma = 0.0;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            sigma += contact_neighborhood.W_ij_[n] * contact_particles_[k]->ParticleVolume(contact_neighborhood.j_[n]);
        }
    }
    repulsion_factor_[index_i] = sigma;
};
//=================================================================================================//
ShellContactFactor::ShellContactFactor(ShellSurfaceContactRelation &solid_body_contact_relation)
    : RepulsionFactorSummation<Base, DataDelegateContact>(solid_body_contact_relation, "RepulsionFactor"),
      solid_(DynamicCast<Solid>(this, sph_body_.getBaseMaterial())),
      kernel_(solid_body_contact_relation.getSPHBody().getSPHAdaptation().getKernel()),
      particle_spacing_(solid_body_contact_relation.getSPHBody().getSPHAdaptation().ReferenceSpacing())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real dp_k = solid_body_contact_relation.contact_bodies_[k]->getSPHAdaptation().ReferenceSpacing();
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
        calibration_factor_.push_back(1.0 / (contact_max + Eps));

        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void ShellContactFactor::interaction(size_t index_i, Real dt)
{
    /** shell contact interaction. */
    Real sigma = 0.0;
    Real contact_density_i = 0.0;

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *contact_Vol_k = contact_Vol_[k];
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
    repulsion_factor_[index_i] = contact_density_i;
}
//=================================================================================================//
ShellSelfContactFactorSummation::ShellSelfContactFactorSummation(ShellSelfContactRelation &self_contact_relation)
    : RepulsionFactorSummation<Base, DataDelegateInner>(self_contact_relation, "SelfRepulsionFactor")
{
}
//=================================================================================================//
void ShellSelfContactFactorSummation::interaction(size_t index_i, Real dt)
{
    Real sigma = 0.0;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma += inner_neighborhood.W_ij_[n] * particles_->ParticleVolume(inner_neighborhood.j_[n]);
    repulsion_factor_[index_i] = sigma;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
