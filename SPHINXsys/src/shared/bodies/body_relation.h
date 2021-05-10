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

#pragma once

#include "base_body.h"
#include "mesh_cell_linked_list.h"
#include "neighbor_relation.h"
#include "base_geometry.h"

namespace SPH
{
	/** a small functor for obtaining particle index for container index */
	struct SPHBodyParticlesIndex
	{
		size_t operator () (size_t particle_index) const { return particle_index; };
	};

	/** a small functor for obtaining particle index for body part container index */
	struct BodyPartParticlesIndex
	{
		IndexVector& body_part_particles_;
		BodyPartParticlesIndex(IndexVector& body_part_particles) : body_part_particles_(body_part_particles) {};
		size_t operator () (size_t particle_entry) const {return body_part_particles_[particle_entry]; };
	};
	/** a small functor for checking particle is body part */
	struct BodyPartParticlesCheck
	{
		IndexVector& body_part_particles_;
		BodyPartParticlesCheck(IndexVector& body_part_particles) : body_part_particles_(body_part_particles) {};
		bool operator () (size_t particle_entry) const 
		{
			if (std::count(body_part_particles_.begin(), body_part_particles_.end(), particle_entry))
			{
				return true;
			}
			else {
				return false;
			}
		};
	};

	/** a small functor for obtaining search range for the simplest case */
		struct SearchRangeSingleResolution
	{
		int operator () (size_t particle_index) const { return 1; };
	};

	/** a small functor for obtaining search range across resolution */
	struct SearchRangeMultiResolution
	{
		int search_range_;
		SearchRangeMultiResolution(SPHBody* body, SPHBody* contact_body) : search_range_(0)
		{
			int body_refinement_level = body->particle_adaptation_->GlobalRefinementLevel();
			int contact_body_refinement_level = contact_body->particle_adaptation_->GlobalRefinementLevel();
			search_range_ = body_refinement_level >= contact_body_refinement_level ? 
				1 : powerN(2, contact_body_refinement_level - body_refinement_level);
		};
		int operator () (size_t particle_index) const { return search_range_; };
	};

	/** a small functor for obtaining search for variable smoothing length */
	struct SearchRangeVariableSmoothingLength
	{
		Real inv_grid_spacing_;
		Kernel* kernel_;
		StdLargeVec<Real>& h_ratio_;
		SearchRangeVariableSmoothingLength(SPHBody* body, MeshCellLinkedList* target_mesh_cell_linked_list) : 
			inv_grid_spacing_(1.0 / target_mesh_cell_linked_list->GridSpacing()), 
			kernel_(body->particle_adaptation_->getKernel()), 
			h_ratio_(*body->base_particles_->getVariableByName<indexScalar, Real>("SmoothingLengthRatio")) {};
		int operator () (size_t particle_index) const 
		{ 
			return 1 + (int)floor(kernel_->CutOffRadius(h_ratio_[particle_index]) * inv_grid_spacing_); 
		};
	};
	
	/**
	 * @class SPHBodyRelation
	 * @brief The relation within a SPH body or with its contact SPH bodies
	 */
	class SPHBodyRelation
	{
	public:
		SPHBody* sph_body_;
		BaseParticles* base_particles_;

		SPHBodyRelation(SPHBody* sph_body);
		virtual ~SPHBodyRelation() {};

		void subscribe_to_body() { sph_body_->body_relations_.push_back(this); };
		virtual void updateConfigurationMemories() = 0;
		virtual void updateConfiguration() = 0;
	};

	/**
	 * @class BaseInnerBodyRelation
	 * @brief The base relation within a SPH body
	 */
	class BaseInnerBodyRelation : public SPHBodyRelation
	{
	protected:
		virtual void resetNeighborhoodCurrentSize();
	public:
		RealBody* real_body_;
		ParticleConfiguration inner_configuration_; /**< inner configuration for the neighbor relations. */

		BaseInnerBodyRelation(RealBody* real_body);
		virtual ~BaseInnerBodyRelation() {};

		virtual void updateConfigurationMemories() override;
	};

	/**
	 * @class InnerBodyRelation
	 * @brief The relation within a SPH body
	 */
	class InnerBodyRelation : public BaseInnerBodyRelation
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;
		SearchRangeSingleResolution get_single_search_range_;
		NeighborRelationInner get_inner_neighbor_;
		MeshCellLinkedList* mesh_cell_linked_list_;

	public:
		InnerBodyRelation(RealBody* real_body);
		virtual ~InnerBodyRelation() {};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class InnerBodyRelationVariableSmoothingLength
	 * @brief The relation within a SPH body with smoothing length adaptation
	 */
	class InnerBodyRelationVariableSmoothingLength : public BaseInnerBodyRelation
	{
	protected:
		size_t total_levels_;
		SPHBodyParticlesIndex get_particle_index_;
		StdVec<SearchRangeVariableSmoothingLength*> get_multi_level_search_range_;
		NeighborRelationInnerVariableSmoothingLength get_inner_neighbor_variable_smoothing_length_;
		StdVec<MeshCellLinkedList*> mesh_cell_linked_list_levels_;
	public:
		InnerBodyRelationVariableSmoothingLength(RealBody* real_body);
		virtual ~InnerBodyRelationVariableSmoothingLength() {};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BaseContactBodyRelation
	 * @brief The base relation between a SPH body and its contact SPH bodies
	 */
	class BaseContactBodyRelation : public SPHBodyRelation
	{
	protected:
		virtual void resetNeighborhoodCurrentSize();
	public:
		RealBodyVector contact_bodies_;
		ContatcParticleConfiguration contact_configuration_; /**< Configurations for particle interaction between bodies. */

		BaseContactBodyRelation(SPHBody* body, RealBodyVector contact_bodies);
		BaseContactBodyRelation(SPHBody* body, BodyPartVector contact_bodyparts);
		virtual ~BaseContactBodyRelation() {};

		virtual void updateConfigurationMemories() override;
	};
	
	/**
	 * @class ContactBodyRelation
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class ContactBodyRelation : public BaseContactBodyRelation
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;
		StdVec<MeshCellLinkedList*> target_mesh_cell_linked_lists_;
		StdVec<SearchRangeMultiResolution*> get_search_ranges_;
		StdVec<NeighborRelationContact*> get_contact_neighbors_;
	public:
		ContactBodyRelation(SPHBody* body, RealBodyVector contact_bodies);
		ContactBodyRelation(SPHBody* body, BodyPartVector contact_bodyparts);
		virtual ~ContactBodyRelation() {};
		void initializaiton();
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SolidContactBodyRelation
	 * @brief The relation between a solid body and its contact solid bodies
	 */
	class SolidContactBodyRelation : public ContactBodyRelation
	{
	protected:
		IndexVector& body_part_particles_;
		BodyPartParticlesIndex get_body_part_particle_index_;

		virtual void resetNeighborhoodCurrentSize() override;
	public:
		ShapeSurfaceLayer body_surface_layer_;

		SolidContactBodyRelation(SPHBody* body, RealBodyVector contact_bodies);
		virtual ~SolidContactBodyRelation() {};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class ReducedInnerBodyRelation
	 * @brief The relation within a reduced SPH body, viz. network
	 */
	class ReducedInnerBodyRelation : public InnerBodyRelation
	{
	protected:
		Tree* tree_;
	public:
		ReducedInnerBodyRelation(RealBody* real_body):InnerBodyRelation(real_body),tree_(real_body->tree_){};
		virtual ~ReducedInnerBodyRelation() {};

		virtual void updateConfiguration() override;
	};
	
	/**
	 * @class PartContactBodyRelation
	 * @brief The relation between a Body part with a SPH body. 
	 */
	class PartContactBodyRelation : public ContactBodyRelation
	{

	public:
		BodyPart* body_part_;
		BodyPartParticlesIndex get_body_part_particle_index_;
		IndexVector& body_part_particles_;

		PartContactBodyRelation(BodyPart* body_part, RealBodyVector contact_bodies);
		virtual ~PartContactBodyRelation() {};

		virtual void updateConfiguration() override;
	};

	/**
	 * @class BodyContactPartRelation
	 * @brief The relation between a Body part with a SPH body. 
	 */
	class BodyContactPartRelation : public ContactBodyRelation
	{

	public:
		BodyPartVector contact_bodyParts_;

		BodyContactPartRelation(RealBody* real_body, BodyPartVector contact_bodyparts);
		virtual ~BodyContactPartRelation() {};

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
		BaseInnerBodyRelation* inner_relation_;
		BaseContactBodyRelation* contact_relation_;
		RealBodyVector contact_bodies_;
		ParticleConfiguration& inner_configuration_;
		ContatcParticleConfiguration& contact_configuration_;

		ComplexBodyRelation(BaseInnerBodyRelation* inner_relation, BaseContactBodyRelation* contact_relation);
		ComplexBodyRelation(RealBody* real_body, RealBodyVector contact_bodies);
		ComplexBodyRelation(BaseInnerBodyRelation* inner_relation, RealBodyVector contact_bodies);
		ComplexBodyRelation(RealBody* real_body, BodyPartVector contact_bodyparts);
		virtual ~ComplexBodyRelation() {
			delete inner_relation_;
			delete contact_relation_;
		};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration()  override;
	};
}
