#include "inner_body_relation.h"
#include "base_particle_dynamics.h"
#include "base_particles.hpp"
#include "cell_linked_list.hpp"

#include "tree_body.h"
namespace SPH
{
//=================================================================================================//
InnerRelation::InnerRelation(RealBody &real_body)
    : BaseInnerRelation(real_body), get_inner_neighbor_(real_body),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
//=================================================================================================//
void InnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh &mesh = cell_linked_list_.getMesh();
    cell_linked_list_.searchNeighborsByMesh(mesh, sph_body_, inner_configuration_,
                                            get_single_search_depth_, get_inner_neighbor_);
}
//=================================================================================================//
AdaptiveInnerRelation::
    AdaptiveInnerRelation(RealBody &real_body)
    : BaseInnerRelation(real_body),
      get_adaptive_inner_neighbor_(real_body),
      multi_level_cell_linked_list_(
          DynamicCast<MultilevelCellLinkedList>(this, real_body.getCellLinkedList()))
{
    Mesh *meshes = multi_level_cell_linked_list_.getMeshes();
    for (size_t l = 0; l != multi_level_cell_linked_list_.ResolutionLevels(); ++l)
    {
        get_multi_level_search_depth_.push_back(
            adaptive_search_depth_ptr_vector_keeper_
                .createPtr<SearchDepthAdaptive>(real_body, meshes[l]));
    }
}
//=================================================================================================//
void AdaptiveInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh *meshes = multi_level_cell_linked_list_.getMeshes();
    for (size_t l = 0; l != multi_level_cell_linked_list_.ResolutionLevels(); ++l)
    {
        multi_level_cell_linked_list_.searchNeighborsByMesh(
            meshes[l], sph_body_, inner_configuration_,
            *get_multi_level_search_depth_[l], get_adaptive_inner_neighbor_);
    }
}
//=================================================================================================//
SelfSurfaceContactRelation::
    SelfSurfaceContactRelation(RealBody &real_body)
    : BaseInnerRelation(real_body),
      body_surface_layer_(real_body),
      body_part_particles_(body_surface_layer_.body_part_particles_),
      get_self_contact_neighbor_(real_body),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
//=================================================================================================//
void SelfSurfaceContactRelation::resetNeighborhoodCurrentSize()
{
    particle_for(execution::ParallelPolicy(), body_part_particles_,
                 [&](size_t index_i)
                 {
                     inner_configuration_[index_i].current_size_ = 0;
                 });
}
//=================================================================================================//
void SelfSurfaceContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh &mesh = cell_linked_list_.getMesh();
    cell_linked_list_.searchNeighborsByMesh(
        mesh, body_surface_layer_, inner_configuration_,
        get_single_search_depth_, get_self_contact_neighbor_);
}
//=================================================================================================//
TreeInnerRelation::TreeInnerRelation(RealBody &real_body)
    : InnerRelation(real_body),
      generative_tree_(DynamicCast<TreeBody>(this, real_body)) {}
//=================================================================================================//
void TreeInnerRelation::updateConfiguration()
{
    generative_tree_.buildParticleConfiguration(inner_configuration_);
}
//=================================================================================================//
ShellInnerRelationWithContactKernel::ShellInnerRelationWithContactKernel(RealBody &real_body, RealBody &contact_body)
    : BaseInnerRelation(real_body),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())),
      get_contact_search_depth_(contact_body, cell_linked_list_.getMesh()),
      get_inner_neighbor_with_contact_kernel_(real_body, contact_body) {}
//=================================================================================================//
void ShellInnerRelationWithContactKernel::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh &mesh = cell_linked_list_.getMesh();
    cell_linked_list_.searchNeighborsByMesh(
        mesh, sph_body_, inner_configuration_,
        get_contact_search_depth_, get_inner_neighbor_with_contact_kernel_);
}
//=================================================================================================//
ShellSelfContactRelation::
    ShellSelfContactRelation(RealBody &real_body)
    : BaseInnerRelation(real_body),
      get_shell_self_contact_neighbor_(real_body),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
//=================================================================================================//
void ShellSelfContactRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh &mesh = cell_linked_list_.getMesh();
    cell_linked_list_.searchNeighborsByMesh(
        mesh, sph_body_, inner_configuration_,
        get_single_search_depth_, get_shell_self_contact_neighbor_);
}
//=================================================================================================//
void AdaptiveSplittingInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh *meshes = multi_level_cell_linked_list_.getMeshes();
    for (size_t l = 0; l != multi_level_cell_linked_list_.ResolutionLevels(); ++l)
    {
        multi_level_cell_linked_list_.searchNeighborsByMesh(
            meshes[l], sph_body_, inner_configuration_,
            *get_multi_level_search_depth_[l], get_adaptive_splitting_inner_neighbor_);
    }
}
//=================================================================================================//
} // namespace SPH
