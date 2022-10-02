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
 * @file 	base_body_relation.h
 * @brief 	Base classes on body and particle topology relations.
 * @author	Xiangyu Hu
 */

#ifndef BASE_BODY_RELATION_H
#define BASE_BODY_RELATION_H

#include "complex_body.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "neighborhood.h"
#include "base_geometry.h"

namespace SPH
{
	/** a small functor for obtaining search range for the simplest case */
	struct SearchDepthSingleResolution
	{
		int operator()(size_t particle_index) const { return 1; };
	};

	/** @brief a small functor for obtaining search depth across resolution
	 * @details Note that the search depth is defined on the target cell linked list.
	 */
	struct SearchDepthMultiResolution
	{
		int search_depth_;
		SearchDepthMultiResolution(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
			: search_depth_(1)
		{
			Real inv_grid_spacing_ = 1.0 / target_cell_linked_list->GridSpacing();
			Kernel *kernel_ = sph_body.sph_adaptation_->getKernel();
			search_depth_ = 1 + (int)floor(kernel_->CutOffRadius() * inv_grid_spacing_);
		};
		int operator()(size_t particle_index) const { return search_depth_; };
	};

	/** @brief a small functor for obtaining search depth for variable smoothing length
	 * @details Note that the search depth is defined on the target cell linked list.
	 */
	struct AdaptiveSearchDepth
	{
		Real inv_grid_spacing_;
		Kernel *kernel_;
		StdLargeVec<Real> &h_ratio_;
		AdaptiveSearchDepth(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
			: inv_grid_spacing_(1.0 / target_cell_linked_list->GridSpacing()),
			  kernel_(sph_body.sph_adaptation_->getKernel()),
			  h_ratio_(*sph_body.getBaseParticles().getVariableByName<Real>("SmoothingLengthRatio")){};
		int operator()(size_t particle_index) const
		{
			return 1 + (int)floor(kernel_->CutOffRadius(h_ratio_[particle_index]) * inv_grid_spacing_);
		};
	};

	/**
	 * @class SPHRelation
	 * @brief The abstract class for all relations within a SPH body or with its contact SPH bodies
	 */
	class SPHRelation
	{
	public:
		SPHBody &sph_body_;
		BaseParticles &base_particles_;
		SPHBody &getDynamicsRange() { return sph_body_; };

		explicit SPHRelation(SPHBody &sph_body);
		virtual ~SPHRelation(){};

		void subscribeToBody() { sph_body_.body_relations_.push_back(this); };
		virtual void updateConfigurationMemories() = 0;
		virtual void updateConfiguration() = 0;
	};

	/**
	 * @class BaseInnerRelation
	 * @brief The abstract relation within a SPH body
	 */
	class BaseInnerRelation : public SPHRelation
	{
	protected:
		virtual void resetNeighborhoodCurrentSize();

	public:
		RealBody *real_body_;
		ParticleConfiguration inner_configuration_; /**< inner configuration for the neighbor relations. */
		explicit BaseInnerRelation(RealBody &real_body);
		virtual ~BaseInnerRelation(){};

		virtual void updateConfigurationMemories() override;
	};

	/**
	 * @class BaseContactRelation
	 * @brief The base relation between a SPH body and its contact SPH bodies
	 */
	class BaseContactRelation : public SPHRelation
	{
	protected:
		UniquePtrKeepers<SearchDepthMultiResolution> search_depth_multi_resolution_ptr_vector_keeper_;
		UniquePtrKeepers<NeighborBuilderContact> neighbor_relation_contact_ptr_vector_keeper_;

	protected:
		StdVec<CellLinkedList *> target_cell_linked_lists_;
		StdVec<SearchDepthMultiResolution *> get_search_depths_;
		StdVec<NeighborBuilderContact *> get_contact_neighbors_;

		virtual void resetNeighborhoodCurrentSize();

	public:
		RealBodyVector contact_bodies_;
		ContactParticleConfiguration contact_configuration_; /**< Configurations for particle interaction between bodies. */

		BaseContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies);
		BaseContactRelation(SPHBody &sph_body, BodyPartVector contact_body_parts);
		virtual ~BaseContactRelation(){};

		virtual void updateConfigurationMemories() override;
	};
}
#endif // BASE_BODY_RELATION_H