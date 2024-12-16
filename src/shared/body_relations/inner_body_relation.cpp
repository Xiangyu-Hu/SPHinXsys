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
    cell_linked_list_.searchNeighborsByParticles(
        sph_body_, inner_configuration_,
        get_single_search_depth_, get_inner_neighbor_);
}
//=================================================================================================//
AdaptiveInnerRelation::
    AdaptiveInnerRelation(RealBody &real_body)
    : BaseInnerRelation(real_body), total_levels_(0),
      get_adaptive_inner_neighbor_(real_body)
{
    MultilevelCellLinkedList &multi_level_cell_linked_list =
        DynamicCast<MultilevelCellLinkedList>(this, real_body.getCellLinkedList());
    cell_linked_list_levels_ = multi_level_cell_linked_list.getMeshLevels();
    total_levels_ = cell_linked_list_levels_.size();
    for (size_t l = 0; l != total_levels_; ++l)
    {
        get_multi_level_search_depth_.push_back(
            adaptive_search_depth_ptr_vector_keeper_
                .createPtr<SearchDepthAdaptive>(real_body, cell_linked_list_levels_[l]));
    }
}
//=================================================================================================//
void AdaptiveInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t l = 0; l != total_levels_; ++l)
    {
        cell_linked_list_levels_[l]->searchNeighborsByParticles(
            sph_body_, inner_configuration_,
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
    cell_linked_list_.searchNeighborsByParticles(
        body_surface_layer_, inner_configuration_,
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
      get_contact_search_depth_(contact_body, &cell_linked_list_),
      get_inner_neighbor_with_contact_kernel_(real_body, contact_body) {}
//=================================================================================================//
void ShellInnerRelationWithContactKernel::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    cell_linked_list_.searchNeighborsByParticles(
        sph_body_, inner_configuration_,
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

    cell_linked_list_.searchNeighborsByParticles(
        sph_body_, inner_configuration_,
        get_single_search_depth_, get_shell_self_contact_neighbor_);
}
//=================================================================================================//
void AdaptiveSplittingInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    for (size_t l = 0; l != total_levels_; ++l)
    {
        cell_linked_list_levels_[l]->searchNeighborsByParticles(
            sph_body_, inner_configuration_,
            *get_multi_level_search_depth_[l], get_adaptive_splitting_inner_neighbor_);
    }
}
//=================================================================================================//
void SplittingInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    cell_linked_list_.searchNeighborsByParticles(
        sph_body_, inner_configuration_,
        get_single_search_depth_, get_splitting_inner_neighbor_);
}
} // namespace SPH
