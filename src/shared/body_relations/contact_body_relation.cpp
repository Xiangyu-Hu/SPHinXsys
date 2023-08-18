#include "contact_body_relation.h"
#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

namespace SPH
{
//=================================================================================================//
ContactRelation::ContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies)
    : ContactRelationCrossResolution(sph_body, contact_bodies)
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_contact_neighbors_.push_back(
            neighbor_builder_contact_ptrs_keeper_.createPtr<NeighborBuilderContact>(
                sph_body_, *contact_bodies_[k]));
    }
}
//=================================================================================================//
void ContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        target_cell_linked_lists_[k]->searchNeighborsByParticles(
            sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
SurfaceContactRelation::SurfaceContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies)
    : ContactRelationCrossResolution(sph_body, contact_bodies),
      body_surface_layer_(shape_surface_ptr_keeper_.createPtr<BodySurfaceLayer>(sph_body)),
      body_part_particles_(body_surface_layer_->body_part_particles_)
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_contact_neighbors_.push_back(
            neighbor_builder_contact_ptrs_keeper_
                .createPtr<NeighborBuilderSurfaceContact>(sph_body_, *contact_bodies_[k]));
    }
}
//=================================================================================================//
void SurfaceContactRelation::resetNeighborhoodCurrentSize()
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        particle_for(execution::ParallelPolicy(), body_part_particles_,
                     [&](size_t index_i)
                     {
                         contact_configuration_[k][index_i].current_size_ = 0;
                     });
    }
}
//=================================================================================================//
void SurfaceContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        target_cell_linked_lists_[k]->searchNeighborsByParticles(
            *body_surface_layer_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
ContactRelationToBodyPart::
    ContactRelationToBodyPart(SPHBody &sph_body, BodyPartVector contact_body_parts_)
    : ContactRelationCrossResolution(sph_body, contact_body_parts_)
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_part_contact_neighbors_.push_back(
            neighbor_builder_contact_ptrs_keeper_
                .createPtr<NeighborBuilderContactBodyPart>(sph_body_, *contact_body_parts_[k]));
    }
}
//=================================================================================================//
void ContactRelationToBodyPart::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        target_cell_linked_lists_[k]->searchNeighborsByParticles(
            sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_part_contact_neighbors_[k]);
    }
}
//=================================================================================================//
AdaptiveContactRelation::AdaptiveContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
    : BaseContactRelation(sph_body, contact_sph_bodies)
{
    cell_linked_list_levels_.resize(contact_bodies_.size());
    get_multi_level_search_range_.resize(contact_bodies_.size());
    get_contact_neighbors_adaptive_.resize(contact_bodies_.size());

    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        cell_linked_list_levels_[k] = contact_bodies_[k]->getCellLinkedList().CellLinkedListLevels();

        for (size_t l = 0; l != cell_linked_list_levels_[k].size(); ++l)
        {
            get_multi_level_search_range_[k].push_back(
                adaptive_search_depth_ptr_vector_keeper_
                    .createPtr<SearchDepthAdaptiveContact>(
                        sph_body_, cell_linked_list_levels_[k][l]));

            get_contact_neighbors_adaptive_[k].push_back(
                neighbor_builder_contact_adaptive_ptr_vector_keeper_
                    .createPtr<NeighborBuilderContactAdaptive>(
                        sph_body, *contact_sph_bodies[k]));
        }
    }
    resizeConfiguration();
}
//=================================================================================================//
void AdaptiveContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        for (size_t l = 0; l != cell_linked_list_levels_[k].size(); ++l)
        {
            cell_linked_list_levels_[k][l]->searchNeighborsByParticles(
                sph_body_, contact_configuration_[k],
                *get_multi_level_search_range_[k][l], *get_contact_neighbors_adaptive_[k][l]);
        }
    }
}
//=================================================================================================//
} // namespace SPH
