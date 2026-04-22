#include "base_body_relation.h"

#include "adaptation.h"
#include "base_particle_dynamics.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
SearchDepthContact::SearchDepthContact(SPHBody &sph_body, Mesh &target_mesh)
    : search_depth_(1)
{
    Real inv_grid_spacing_ = 1.0 / target_mesh.GridSpacing();
    Kernel *kernel_ = sph_body.getSPHAdaptation().getKernel();
    search_depth_ = 1 + (int)floor(kernel_->CutOffRadius() * inv_grid_spacing_);
}
//=================================================================================================//
SearchDepthAdaptive::SearchDepthAdaptive(SPHBody &sph_body, Mesh &target_mesh)
    : inv_grid_spacing_(1.0 / target_mesh.GridSpacing()),
      kernel_(sph_body.getSPHAdaptation().getKernel()),
      h_ratio_(sph_body.getBaseParticles().getVariableDataByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
SearchDepthAdaptiveContact::SearchDepthAdaptiveContact(SPHBody &sph_body, Mesh &target_mesh)
    : inv_grid_spacing_(1.0 / target_mesh.GridSpacing()),
      sph_adaptation_(sph_body.getSPHAdaptation()),
      kernel_(*sph_adaptation_.getKernel()) {};
//=================================================================================================//
int SearchDepthAdaptiveContact::operator()(size_t particle_index) const
{
    return 1 + (int)floor(
                   kernel_.CutOffRadius(sph_adaptation_.SmoothingLengthRatio(particle_index)) *
                   inv_grid_spacing_);
}
//=================================================================================================//
RealBodyVector BodyPartsToRealBodies(BodyPartVector body_parts)
{
    RealBodyVector real_bodies;

    for (size_t k = 0; k != body_parts.size(); ++k)
    {
        real_bodies.push_back(DynamicCast<RealBody>(body_parts[k], &body_parts[k]->getSPHBody()));
    }
    return real_bodies;
}
//=================================================================================================//
SPHRelation::SPHRelation(SPHBody &sph_body)
    : sph_body_(sph_body),
      base_particles_(sph_body.getBaseParticles()) {}
//=================================================================================================//
BaseInnerRelation::BaseInnerRelation(RealBody &real_body)
    : SPHRelation(real_body), real_body_(&real_body)
{
    subscribeToBody();
    inner_configuration_.resize(base_particles_.ParticlesBound(), Neighborhood());
}
//=================================================================================================//
void BaseInnerRelation::resetNeighborhoodCurrentSize()
{
    tbb::parallel_for(
        IndexRange(0, base_particles_.TotalRealParticles()),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                inner_configuration_[num].current_size_ = 0;
            }
        },
        ap);
}
//=================================================================================================//
BaseContactRelation::BaseContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
    : SPHRelation(sph_body), contact_bodies_(contact_sph_bodies)
{
    subscribeToBody();
    contact_configuration_.resize(contact_bodies_.size());
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        const std::string name = contact_bodies_[k]->getName();
        contact_particles_.push_back(&contact_bodies_[k]->getBaseParticles());
        contact_adaptations_.push_back(&contact_bodies_[k]->getSPHAdaptation());
        contact_configuration_[k].resize(base_particles_.ParticlesBound(), Neighborhood());
    }
}
//=================================================================================================//
void BaseContactRelation::resetNeighborhoodCurrentSize()
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        tbb::parallel_for(
            IndexRange(0, base_particles_.TotalRealParticles()),
            [&](const IndexRange &r)
            {
                for (size_t num = r.begin(); num != r.end(); ++num)
                {
                    contact_configuration_[k][num].current_size_ = 0;
                }
            },
            ap);
    }
}
//=================================================================================================//
} // namespace SPH
