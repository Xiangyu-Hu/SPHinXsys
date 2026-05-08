/**
 * @file em_matched_stencil_body_relations.h
 * @brief Optional inner/contact BodyRelation variants that align cell-stencil depth
 *        between InnerRelation and ContactRelation for EM split/merged consistency studies.
 */
#ifndef EM_MATCHED_STENCIL_BODY_RELATIONS_H
#define EM_MATCHED_STENCIL_BODY_RELATIONS_H

#include "contact_body_relation.h"
#include "inner_body_relation.h"

namespace SPH
{
namespace electromagnetics
{

/**
 * @brief Inner neighbor search using SearchDepthContact (same rule as ContactRelation).
 */
class InnerRelationContactSearchDepth : public InnerRelation
{
    SearchDepthContact depth_;

  public:
    explicit InnerRelationContactSearchDepth(RealBody &real_body)
        : InnerRelation(real_body), depth_(real_body, this->getCellLinkedList().getMesh())
    {
    }

    void updateConfiguration() override
    {
        resetNeighborhoodCurrentSize();
        Mesh &mesh = cell_linked_list_.getMesh();
        cell_linked_list_.searchNeighborsByMesh(mesh, sph_body_, inner_configuration_, depth_, get_inner_neighbor_);
    }
};

/**
 * @brief Contact neighbor search using SearchDepthSingleResolution (same rule as InnerRelation).
 */
class ContactRelationInnerSearchDepth : public ContactRelation
{
  public:
    explicit ContactRelationInnerSearchDepth(SPHBody &sph_body, RealBodyVector contact_bodies)
        : ContactRelation(sph_body, contact_bodies)
    {
    }

    void updateConfiguration() override
    {
        resetNeighborhoodCurrentSize();
        SearchDepthSingleResolution inner_like_depth;
        for (size_t k = 0; k < contact_bodies_.size(); ++k)
        {
            Mesh &mesh = target_cell_linked_lists_[k]->getMesh();
            target_cell_linked_lists_[k]->searchNeighborsByMesh(mesh, sph_body_, contact_configuration_[k],
                                                                inner_like_depth, *get_contact_neighbors_[k]);
        }
    }
};

} // namespace electromagnetics
} // namespace SPH

#endif // EM_MATCHED_STENCIL_BODY_RELATIONS_H
