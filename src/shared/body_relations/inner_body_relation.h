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
 * @file 	inner_body_relation.h
 * @brief 	The topological relations within the body.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef INNER_BODY_RELATION_H
#define INNER_BODY_RELATION_H

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class InnerRelation
 * @brief The first concrete relation within a SPH body
 */
class InnerRelation : public BaseInnerRelation
{
  protected:
    SearchDepthSingleResolution get_single_search_depth_;
    NeighborBuilderInner get_inner_neighbor_;
    CellLinkedList &cell_linked_list_;

  public:
    explicit InnerRelation(RealBody &real_body);
    virtual ~InnerRelation(){};

    CellLinkedList &getCellLinkedList() { return cell_linked_list_; };
    virtual void updateConfiguration() override;
};

/**
 * @class AdaptiveInnerRelation
 * @brief The relation within a SPH body with smoothing length adaptation
 */
class AdaptiveInnerRelation : public BaseInnerRelation
{
  private:
    UniquePtrsKeeper<SearchDepthAdaptive> adaptive_search_depth_ptr_vector_keeper_;

  protected:
    StdVec<SearchDepthAdaptive *> get_multi_level_search_depth_;
    NeighborBuilderInnerAdaptive get_adaptive_inner_neighbor_;
    MultilevelCellLinkedList &multi_level_cell_linked_list_;

  public:
    explicit AdaptiveInnerRelation(RealBody &real_body);
    virtual ~AdaptiveInnerRelation(){};

    virtual void updateConfiguration() override;
};

/**
 * @class SelfSurfaceContactRelation
 * @brief The relation for self contact of a solid body
 */
class SelfSurfaceContactRelation : public BaseInnerRelation
{
  public:
    BodySurfaceLayer body_surface_layer_;

    explicit SelfSurfaceContactRelation(RealBody &real_body);
    virtual ~SelfSurfaceContactRelation(){};
    virtual void updateConfiguration() override;

  protected:
    IndexVector &body_part_particles_;
    SearchDepthSingleResolution get_single_search_depth_;
    NeighborBuilderSelfContact get_self_contact_neighbor_;
    CellLinkedList &cell_linked_list_;

    virtual void resetNeighborhoodCurrentSize() override;
};

class TreeBody;
/**
 * @class TreeInnerRelation
 * @brief The relation within a reduced SPH body, viz. network
 */
class TreeInnerRelation : public InnerRelation
{
  protected:
    TreeBody &generative_tree_;

  public:
    explicit TreeInnerRelation(RealBody &real_body);
    virtual ~TreeInnerRelation(){};

    virtual void updateConfiguration() override;
};

/**
 * @class ShellInnerRelationWithContactKernel
 * @brief Shell inner relation with the cut-off radius of the contact body
 *  This class is used in fluid-shell interaction problems to compute shell curvature with the cut-off radius of fluid
 */
class ShellInnerRelationWithContactKernel : public BaseInnerRelation
{
  private:
    CellLinkedList &cell_linked_list_;
    SearchDepthContact get_contact_search_depth_;
    ShellNeighborBuilderInnerWithContactKernel get_inner_neighbor_with_contact_kernel_;

  public:
    explicit ShellInnerRelationWithContactKernel(RealBody &real_body, RealBody &contact_body);
    void updateConfiguration() override;
};

/**
 * @class ShellSelfContactRelation
 * @brief The relation for self contact of a shell
 */
class ShellSelfContactRelation : public BaseInnerRelation
{
  public:
    explicit ShellSelfContactRelation(RealBody &real_body);
    void updateConfiguration() override;

  private:
    SearchDepthSingleResolution get_single_search_depth_;
    NeighborBuilderShellSelfContact get_shell_self_contact_neighbor_;
    CellLinkedList &cell_linked_list_;
};

/**
 * @class AdaptiveSplittingInnerRelation
 * @brief The relation within a SPH body with smoothing length adaptation for splitting algorithm
 *        a particle can only see neighbors with ascending ids or higher levels
 */
class AdaptiveSplittingInnerRelation : public AdaptiveInnerRelation
{
  public:
    explicit AdaptiveSplittingInnerRelation(RealBody &real_body)
        : AdaptiveInnerRelation(real_body),
          get_adaptive_splitting_inner_neighbor_(real_body){};
    void updateConfiguration() override;

  private:
    NeighborBuilderSplitInnerAdaptive get_adaptive_splitting_inner_neighbor_;
};
} // namespace SPH
#endif // INNER_BODY_RELATION_H