#include "contact_body_relation.h"
#include "all_particles.h"
#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

#include <numeric>

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
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
ShellSurfaceContactRelation::ShellSurfaceContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies)
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
void ShellSurfaceContactRelation::resetNeighborhoodCurrentSize()
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
void ShellSurfaceContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, *body_surface_layer_, contact_configuration_[k],
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
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_part_contact_neighbors_[k]);
    }
}
//=================================================================================================//
AdaptiveContactRelation::AdaptiveContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
    : BaseContactRelation(sph_body, contact_sph_bodies)
{
    cell_linked_lists_.resize(contact_bodies_.size());
    get_multi_level_search_range_.resize(contact_bodies_.size());
    get_contact_neighbors_adaptive_.resize(contact_bodies_.size());

    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        cell_linked_lists_[k] = &contact_bodies_[k]->getCellLinkedList();
        StdVec<Mesh *> &meshes = cell_linked_lists_[k]->getMeshes();
        for (size_t l = 0; l != meshes.size(); ++l)
        {
            get_multi_level_search_range_[k].push_back(
                adaptive_search_depth_ptr_vector_keeper_
                    .createPtr<SearchDepthAdaptiveContact>(sph_body_, *meshes[l]));

            get_contact_neighbors_adaptive_[k].push_back(
                neighbor_builder_contact_adaptive_ptr_vector_keeper_
                    .createPtr<NeighborBuilderContactAdaptive>(sph_body, *contact_sph_bodies[k]));
        }
    }
}
//=================================================================================================//
void AdaptiveContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        StdVec<Mesh *> &meshes = cell_linked_lists_[k]->getMeshes();
        for (size_t l = 0; l != meshes.size(); ++l)
        {
            cell_linked_lists_[k]->searchNeighborsByMesh(
                *meshes[l], sph_body_, contact_configuration_[k],
                *get_multi_level_search_range_[k][l], *get_contact_neighbors_adaptive_[k][l]);
        }
    }
}
//=================================================================================================//
ContactRelationFromShellToFluid::ContactRelationFromShellToFluid(SPHBody &sph_body, RealBodyVector contact_bodies,
                                                                 const StdVec<bool> &normal_corrections)
    : ContactRelationCrossResolution(sph_body, contact_bodies)
{
    if (contact_bodies.size() != normal_corrections.size())
    {
        throw std::runtime_error("ContactRelationFromShell: sizes of normal_corrections and contact_bodies are different!");
    }

    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_shell_contact_neighbors_.push_back(
            neighbor_builder_contact_from_shell_ptrs_keeper_.createPtr<NeighborBuilderContactFromShellToFluid>(
                sph_body_, *contact_bodies_[k], normal_corrections[k]));
    }
}
//=================================================================================================//
void ContactRelationFromShellToFluid::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_shell_contact_neighbors_[k]);
    }
}
//=================================================================================================//
ContactRelationFromFluidToShell::ContactRelationFromFluidToShell(SPHBody &sph_body, RealBodyVector contact_bodies,
                                                                 const StdVec<bool> &normal_corrections)
    : ContactRelationCrossResolution(sph_body, contact_bodies)
{
    if (contact_bodies.size() != normal_corrections.size())
    {
        throw std::runtime_error("ContactRelationToShell: sizes of normal_corrections and contact_bodies are different!");
    }

    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        get_contact_neighbors_.push_back(
            neighbor_builder_contact_to_shell_ptrs_keeper_.createPtr<NeighborBuilderContactFromFluidToShell>(
                sph_body_, *contact_bodies_[k], normal_corrections[k]));
    }
}
//=================================================================================================//
void ContactRelationFromFluidToShell::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, sph_body_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
SurfaceContactRelation::SurfaceContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies, StdVec<bool> normal_corrections)
    : ContactRelationCrossResolution(sph_body, std::move(contact_bodies)),
      body_surface_layer_(shape_surface_ptr_keeper_.createPtr<BodySurfaceLayer>(sph_body)),
      body_part_particles_(body_surface_layer_->body_part_particles_)
{
    // Check if the source body is a shell body
    // Fix invalid list of particles ids with shell (in case body shape is absent by the time of construction)
    if (dynamic_cast<SurfaceParticles *>(&sph_body.getBaseParticles()))
    {
        body_part_particles_.resize(sph_body.getBaseParticles().TotalRealParticles());
        std::iota(body_part_particles_.begin(), body_part_particles_.end(), 0);
    }
    // check if the normal_corrections size is correct
    if (!normal_corrections.empty() && normal_corrections.size() != contact_bodies_.size())
        throw std::runtime_error("SolidShellSurfaceContactRelation: sizes of normal_corrections and shell_contact_bodies are different!");

    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        if (dynamic_cast<SurfaceParticles *>(&contact_bodies_[k]->getBaseParticles()) == nullptr)
        {
            // solid neighbors
            get_contact_neighbors_.push_back(
                solid_neighbor_builder_contact_ptrs_keeper_
                    .createPtr<NeighborBuilderSurfaceContactFromSolid>(sph_body_, *contact_bodies_[k]));
        }
        else
        {
            // shell neighbors
            const bool normal_correction = normal_corrections.empty() ? false : normal_corrections[k];
            get_contact_neighbors_.push_back(
                shell_neighbor_builder_contact_ptrs_keeper_
                    .createPtr<NeighborBuilderSurfaceContactFromShell>(sph_body_, *contact_bodies_[k], normal_correction));
        }
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
        Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
        target_cell_linked_lists_[k]->searchNeighborsByMesh(
            mesh, *body_surface_layer_, contact_configuration_[k],
            *get_search_depths_[k], *get_contact_neighbors_[k]);
    }
}
//=================================================================================================//
} // namespace SPH
