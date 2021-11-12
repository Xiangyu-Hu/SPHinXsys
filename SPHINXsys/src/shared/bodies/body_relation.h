/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	body_relation.h
 * @brief 	The topological relations between bodies are described here.
 * @author	Xiangyu Hu
 * @version 0.3.0
 * 			-- Add reduced body relation for network. Chi ZHANG
 */

#ifndef BODY_RELATION_H
#define BODY_RELATION_H

#include "base_body.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "neighbor_relation.h"
#include "base_geometry.h"

namespace SPH
{
	/** a small functor for obtaining particle index for container index */
	struct SPHBodyParticlesIndex
	{
		size_t operator()(size_t particle_index) const { return particle_index; };
	};

	/** a small functor for obtaining particle index for body part container index */
	struct BodyPartParticlesIndex
	{
		IndexVector &body_part_particles_;
		BodyPartParticlesIndex(IndexVector &body_part_particles) : body_part_particles_(body_part_particles){};
		size_t operator()(size_t particle_entry) const { return body_part_particles_[particle_entry]; };
	};

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
		SearchDepthMultiResolution(SPHBody *body, CellLinkedList *target_cell_linked_list) : search_depth_(1)
		{
			Real inv_grid_spacing_ = 1.0 / target_cell_linked_list->GridSpacing();
			Kernel *kernel_ = body->particle_adaptation_->getKernel();
			search_depth_ = 1 + (int)floor(kernel_->CutOffRadius() * inv_grid_spacing_);
		};
		int operator()(size_t particle_index) const { return search_depth_; };
	};

	/** @brief a small functor for obtaining search depth for variable smoothing length 
	 * @details Note that the search depth is defined on the target cell linked list.
	 */
	struct SearchDepthVariableSmoothingLength
	{
		Real inv_grid_spacing_;
		Kernel *kernel_;
		StdLargeVec<Real> &h_ratio_;
		SearchDepthVariableSmoothingLength(SPHBody *body, CellLinkedList *target_cell_linked_list)
			: inv_grid_spacing_(1.0 / target_cell_linked_list->GridSpacing()),
			  kernel_(body->particle_adaptation_->getKernel()),
			  h_ratio_(*body->base_particles_->getVariableByName<indexScalar, Real>("SmoothingLengthRatio")){};
		int operator()(size_t particle_index) const
		{
			return 1 + (int)floor(kernel_->CutOffRadius(h_ratio_[particle_index]) * inv_grid_spacing_);
		};
	};

	/**
	 * @class SPHBodyRelation
	 * @brief The abstract class for all relations within a SPH body or with its contact SPH bodies
	 */
	class SPHBodyRelation
	{
	public:
		SPHBody *sph_body_;
		BaseParticles *base_particles_;

		SPHBodyRelation(SPHBody *sph_body);
		virtual ~SPHBodyRelation(){};

		void subscribeToBody() { sph_body_->body_relations_.push_back(this); };
		virtual void updateConfigurationMemories() = 0;
		virtual void updateConfiguration() = 0;
	};

	/**
	 * @class BaseBodyRelationInner
	 * @brief The abstract relation within a SPH body
	 */
	class BaseBodyRelationInner : public SPHBodyRelation
	{
	protected:
		virtual void resetNeighborhoodCurrentSize();

	public:
		RealBody *real_body_;
		ParticleConfiguration inner_configuration_; /**< inner configuration for the neighbor relations. */

		BaseBodyRelationInner(RealBody *real_body);
		virtual ~BaseBodyRelationInner(){};

		virtual void updateConfigurationMemories() override;
	};

	/**
	 * @class BodyRelationInner
	 * @brief The first concrete relation within a SPH body
	 */
	class BodyRelationInner : public BaseBodyRelationInner
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;
		SearchDepthSingleResolution get_single_search_depth_;
		NeighborRelationInner get_inner_neighbor_;
		CellLinkedList *cell_linked_list_;

	public:
		BodyRelationInner(RealBody *real_body);
		virtual ~BodyRelationInner(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BodyRelationInnerVariableSmoothingLength
	 * @brief The relation within a SPH body with smoothing length adaptation
	 */
	class BodyRelationInnerVariableSmoothingLength : public BaseBodyRelationInner
	{
	protected:
		size_t total_levels_;
		SPHBodyParticlesIndex get_particle_index_;
		StdVec<SearchDepthVariableSmoothingLength *> get_multi_level_search_depth_;
		NeighborRelationInnerVariableSmoothingLength get_inner_neighbor_variable_smoothing_length_;
		StdVec<CellLinkedList *> cell_linked_list_levels_;

	public:
		BodyRelationInnerVariableSmoothingLength(RealBody *real_body);
		virtual ~BodyRelationInnerVariableSmoothingLength(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class SolidBodyRelationSelfContact
	 * @brief The relation for self contact of a solid body
	 */
	class SolidBodyRelationSelfContact : public BaseBodyRelationInner
	{
	public:
		ShapeSurfaceLayer body_surface_layer_;

		SolidBodyRelationSelfContact(RealBody *real_body);
		virtual ~SolidBodyRelationSelfContact(){};

		virtual void updateConfiguration() override;

	protected:
		IndexVector &body_part_particles_;
		BodyPartParticlesIndex get_body_part_particle_index_;
		SearchDepthSingleResolution get_single_search_depth_;
		NeighborRelationSelfContact get_self_contact_neighbor_;
		CellLinkedList *cell_linked_list_;

		virtual void resetNeighborhoodCurrentSize() override;
	};

	/**
	 * @class BaseBodyRelationContact
	 * @brief The base relation between a SPH body and its contact SPH bodies
	 */
	class BaseBodyRelationContact : public SPHBodyRelation
	{
	protected:
		StdVec<CellLinkedList *> target_cell_linked_lists_;
		StdVec<SearchDepthMultiResolution *> get_search_depths_;
		StdVec<NeighborRelationContact *> get_contact_neighbors_;

		virtual void resetNeighborhoodCurrentSize();

	public:
		RealBodyVector contact_bodies_;
		ContatcParticleConfiguration contact_configuration_; /**< Configurations for particle interaction between bodies. */

		BaseBodyRelationContact(SPHBody *body, RealBodyVector contact_bodies);
		BaseBodyRelationContact(SPHBody *body, BodyPartVector contact_body_parts);
		virtual ~BaseBodyRelationContact(){};

		virtual void updateConfigurationMemories() override;
	};

	/**
	 * @class BodyRelationContact
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class BodyRelationContact : public BaseBodyRelationContact
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;

		void initialization();

	public:
		BodyRelationContact(SPHBody *body, RealBodyVector contact_bodies);
		BodyRelationContact(SPHBody *body, BodyPartVector contact_body_parts);
		virtual ~BodyRelationContact(){};
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SolidBodyRelationContact
	 * @brief The relation between a solid body and its contact solid bodies
	 */
	class SolidBodyRelationContact : public BaseBodyRelationContact
	{
	public:
		ShapeSurfaceLayer body_surface_layer_;

		SolidBodyRelationContact(SPHBody *sph_body, RealBodyVector contact_bodies);
		SolidBodyRelationContact(SolidBodyRelationSelfContact *solid_body_relation_self_contact,
								 RealBodyVector contact_bodies);
		virtual ~SolidBodyRelationContact(){};

		virtual void updateConfiguration() override;
	protected:
		IndexVector &body_part_particles_;
		BodyPartParticlesIndex get_body_part_particle_index_;

		void initialization();
		virtual void resetNeighborhoodCurrentSize() override;
	};
	
	/**
	 * @class GenerativeBodyRelationInner
	 * @brief The relation within a reduced SPH body, viz. network
	 */
	class GenerativeBodyRelationInner : public BodyRelationInner
	{
	protected:
		GenerativeStructure *generative_structure_;

	public:
		GenerativeBodyRelationInner(RealBody *real_body)
			: BodyRelationInner(real_body),
			  generative_structure_(real_body->generative_structure_){};
		virtual ~GenerativeBodyRelationInner(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BodyPartRelationContact
	 * @brief The relation between a Body part with a SPH body. 
	 */
	class BodyPartRelationContact : public BodyRelationContact
	{

	public:
		BodyPart *body_part_;
		IndexVector &body_part_particles_;
		BodyPartParticlesIndex get_body_part_particle_index_;

		BodyPartRelationContact(BodyPart *body_part, RealBodyVector contact_bodies);
		virtual ~BodyPartRelationContact(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BodyRelationContactToBodyPart
	 * @brief The relation between a SPH body and a vector of body parts. 
	 */
	class BodyRelationContactToBodyPart : public BodyRelationContact
	{

	public:
		BodyPartVector contact_body_parts_;
		StdVec<NeighborRelationContactBodyPart *> get_part_contact_neighbors_;

		BodyRelationContactToBodyPart(RealBody *real_body, BodyPartVector contact_body_parts);
		virtual ~BodyRelationContactToBodyPart(){};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class ComplexBodyRelation
	 * @brief The relation combined an inner and a contactbody relation.
	 * The interaction is in a inner-boundary-condition fashion. Here inner interaction is
	 * different from contact interaction.
	 */
	class ComplexBodyRelation : public SPHBodyRelation
	{
	public:
		BaseBodyRelationInner *inner_relation_;
		BaseBodyRelationContact *contact_relation_;
		RealBodyVector contact_bodies_;
		ParticleConfiguration &inner_configuration_;
		ContatcParticleConfiguration &contact_configuration_;

		ComplexBodyRelation(BaseBodyRelationInner *inner_relation, BaseBodyRelationContact *contact_relation);
		ComplexBodyRelation(RealBody *real_body, RealBodyVector contact_bodies);
		ComplexBodyRelation(BaseBodyRelationInner *inner_relation, RealBodyVector contact_bodies);
		ComplexBodyRelation(RealBody *real_body, BodyPartVector contact_body_parts);
		virtual ~ComplexBodyRelation()
		{
			delete inner_relation_;
			delete contact_relation_;
		};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};
}
#endif //BODY_RELATION_H