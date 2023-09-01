/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	contact_body_relation.h
 * @brief 	The topological relations between bodies,
 * 			and the corresponding local topologies (particle configurations)
 * 			are constructed in these classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONTACT_BODY_RELATION_H
#define CONTACT_BODY_RELATION_H

#include "base_body_relation.h"
#include "inner_body_relation.h"

namespace SPH
{
/**
 * @class ContactRelationCrossResolution
 * @brief The relation between a SPH body and its contact SPH bodies
 */
class ContactRelationCrossResolution : public BaseContactRelation
{
  protected:
    UniquePtrsKeeper<SearchDepthContact> search_depth_ptrs_keeper_;

  public:
    template <typename... Args>
    ContactRelationCrossResolution(SPHBody &sph_body, Args &&...args)
        : BaseContactRelation(sph_body, std::forward<Args>(args)...)
    {
        for (size_t k = 0; k != contact_bodies_.size(); ++k)
        {
            CellLinkedList *target_cell_linked_list =
                DynamicCast<CellLinkedList>(this, &contact_bodies_[k]->getCellLinkedList());
            target_cell_linked_lists_.push_back(target_cell_linked_list);
            get_search_depths_.push_back(
                search_depth_ptrs_keeper_.createPtr<SearchDepthContact>(
                    sph_body_, target_cell_linked_list));
        }
        resizeConfiguration();
    };
    virtual ~ContactRelationCrossResolution(){};

  protected:
    StdVec<CellLinkedList *> target_cell_linked_lists_;
    StdVec<SearchDepthContact *> get_search_depths_;
};

/**
 * @class ContactRelation
 * @brief The relation between a SPH body and its contact SPH bodies
 */
class ContactRelation : public ContactRelationCrossResolution
{
  protected:
    UniquePtrsKeeper<NeighborBuilderContact> neighbor_builder_contact_ptrs_keeper_;

  public:
    ContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies);
    virtual ~ContactRelation(){};
    virtual void updateConfiguration() override;

  protected:
    StdVec<NeighborBuilderContact *> get_contact_neighbors_;
};

/**
 * @class SurfaceContactRelation
 * @brief The relation between a solid body and its contact solid bodies
 */
class SurfaceContactRelation : public ContactRelationCrossResolution
{
  protected:
    UniquePtrsKeeper<NeighborBuilderSurfaceContact> neighbor_builder_contact_ptrs_keeper_;
    UniquePtrKeeper<BodySurfaceLayer> shape_surface_ptr_keeper_;

  public:
    BodySurfaceLayer *body_surface_layer_;

    SurfaceContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies);
    SurfaceContactRelation(SelfSurfaceContactRelation &solid_body_relation_self_contact,
                           RealBodyVector contact_bodies)
        : SurfaceContactRelation(*solid_body_relation_self_contact.real_body_, contact_bodies){};
    virtual ~SurfaceContactRelation(){};
    virtual void updateConfiguration() override;

  protected:
    IndexVector &body_part_particles_;
    StdVec<NeighborBuilderSurfaceContact *> get_contact_neighbors_;

    virtual void resetNeighborhoodCurrentSize() override;
};

/**
 * @class ContactRelationToBodyPart
 * @brief The relation between a SPH body and a vector of body parts.
 */
class ContactRelationToBodyPart : public ContactRelationCrossResolution
{
  protected:
    UniquePtrsKeeper<NeighborBuilderContactBodyPart> neighbor_builder_contact_ptrs_keeper_;

  public:
    StdVec<NeighborBuilderContactBodyPart *> get_part_contact_neighbors_;

    ContactRelationToBodyPart(SPHBody &sph_body, BodyPartVector contact_body_parts_);
    virtual ~ContactRelationToBodyPart(){};

    virtual void updateConfiguration() override;
};

/**
 * @class AdaptiveContactRelation
 * @brief The relation between a SPH body and its contact SPH bodies in Multi-level Mesh
 */
class AdaptiveContactRelation : public BaseContactRelation
{
  private:
    UniquePtrsKeeper<SearchDepthAdaptiveContact> adaptive_search_depth_ptr_vector_keeper_;
    UniquePtrsKeeper<NeighborBuilderContactAdaptive> neighbor_builder_contact_adaptive_ptr_vector_keeper_;

  protected:
    StdVec<StdVec<SearchDepthAdaptiveContact *>> get_multi_level_search_range_;
    StdVec<StdVec<CellLinkedList *>> cell_linked_list_levels_;
    StdVec<StdVec<NeighborBuilderContactAdaptive *>> get_contact_neighbors_adaptive_;

  public:
    AdaptiveContactRelation(SPHBody &body, RealBodyVector contact_bodies);
    virtual ~AdaptiveContactRelation(){};

    virtual void updateConfiguration() override;
};
} // namespace SPH
#endif // CONTACT_BODY_RELATION_H