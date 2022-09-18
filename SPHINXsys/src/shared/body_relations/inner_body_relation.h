/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	inner_body_relation.h
 * @brief 	The topological relations within the body.
 * @author	Xiangyu Hu
 */

#ifndef INNER_BODY_RELATION_H
#define INNER_BODY_RELATION_H

#include "base_body_relation.h"

namespace SPH
{
	/**
	 * @class BodyRelationInner
	 * @brief The first concrete relation within a SPH body
	 */
	class BodyRelationInner : public BaseBodyRelationInner
	{
	protected:
		SearchDepthSingleResolution get_single_search_depth_;
		NeighborBuilderInner get_inner_neighbor_;
		CellLinkedList *cell_linked_list_;

	public:
		explicit BodyRelationInner(RealBody &real_body);
		virtual ~BodyRelationInner(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BodyRelationInnerVariableSmoothingLength
	 * @brief The relation within a SPH body with smoothing length adaptation
	 */
	class BodyRelationInnerVariableSmoothingLength : public BaseBodyRelationInner
	{
	private:
		UniquePtrKeepers<SearchDepthVariableSmoothingLength> search_variable_smoothinglength_ptr_vector_keeper_;

	protected:
		size_t total_levels_;
		StdVec<SearchDepthVariableSmoothingLength *> get_multi_level_search_depth_;
		NeighborBuilderInnerVariableSmoothingLength get_inner_neighbor_variable_smoothing_length_;
		StdVec<CellLinkedList *> cell_linked_list_levels_;

	public:
		explicit BodyRelationInnerVariableSmoothingLength(RealBody &real_body);
		virtual ~BodyRelationInnerVariableSmoothingLength(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class SolidBodyRelationSelfContact
	 * @brief The relation for self contact of a solid body
	 * TODO: better called BodySurfaceSelfContact
	 */
	class SolidBodyRelationSelfContact : public BaseBodyRelationInner
	{
	public:
		BodySurfaceLayer body_surface_layer_;

		explicit SolidBodyRelationSelfContact(RealBody &real_body);
		virtual ~SolidBodyRelationSelfContact(){};
		BodyPartByParticle &getDynamicsRange() { return body_surface_layer_; };

		virtual void updateConfiguration() override;

	protected:
		IndexVector &body_part_particles_;
		SearchDepthSingleResolution get_single_search_depth_;
		NeighborBuilderSelfContact get_self_contact_neighbor_;
		CellLinkedList *cell_linked_list_;

		virtual void resetNeighborhoodCurrentSize() override;
	};

	/**
	 * @class TreeBodyRelationInner
	 * @brief The relation within a reduced SPH body, viz. network
	 */
	class TreeBodyRelationInner : public BodyRelationInner
	{
	protected:
		TreeBody &generative_tree_;

	public:
		explicit TreeBodyRelationInner(RealBody &real_body)
			: BodyRelationInner(real_body),
			  generative_tree_(DynamicCast<TreeBody>(this, real_body)){};
		virtual ~TreeBodyRelationInner(){};

		virtual void updateConfiguration() override;
	};
}
#endif //INNER_BODY_RELATION_H