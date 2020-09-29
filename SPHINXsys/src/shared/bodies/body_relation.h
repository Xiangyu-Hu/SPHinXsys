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
 */

#pragma once
#include "neighbor_relation.h"
#include "base_body.h"

namespace SPH
{
	/**
	 * @class SPHBodyBaseRelation
	 * @brief The relation within a SPH body or with its contact SPH bodies
	 */
	class SPHBodyBaseRelation
	{
	public:
		SPHBody* sph_body_;
		SplitCellLists& split_cell_lists_;
		BaseParticles* base_particles_;
		BaseMeshCellLinkedList* mesh_cell_linked_list_;

		SPHBodyBaseRelation(SPHBody* sph_body);
		virtual ~SPHBodyBaseRelation() {};

		void subscribe_to_body() { sph_body_->body_relations_.push_back(this); };
		virtual void updateConfigurationMemories() = 0;
		virtual void updateConfiguration() = 0;
	protected:
		virtual void createNeighborRelation(Neighborhood& neighborhood,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
		virtual void initializeNeighborRelation(Neighborhood& neighborhood, size_t current_count_of_neighbors, 
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);

	};

	/**
	 * @class SPHBodyInnerRelation
	 * @brief The relation within a SPH body
	 */
	class SPHBodyInnerRelation : public SPHBodyBaseRelation
	{
	public:
		/** inner configuration for the neighbor relations. */
		ParticleConfiguration inner_configuration_;

		SPHBodyInnerRelation(SPHBody* sph_body);
		virtual ~SPHBodyInnerRelation() {};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SPHBodyContactRelation
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class SPHBodyContactRelation : public SPHBodyBaseRelation
	{
	protected:
		StdVec<BaseMeshCellLinkedList*> target_mesh_cell_linked_lists_;
	public:
		SPHBodyVector contact_sph_bodies_;

		/** Configurations for particle interaction between bodies. */
		ContatcParticleConfiguration contact_configuration_;

		SPHBodyContactRelation(SPHBody* body, SPHBodyVector relation_bodies);
		virtual ~SPHBodyContactRelation() {};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SPHBodyComplexRelation
	 * @brief The relation within a SPH body and with its contact SPH bodies.
	 * The interaction is in a inner-boundary-condition fashion. Here inner interaction is
	 * different from conact interaction
	 */
	class SPHBodyComplexRelation : public SPHBodyBaseRelation
	{
	protected:
		SPHBodyInnerRelation* inner_relation_;
		SPHBodyContactRelation* contact_relation_;
	public:
		SPHBodyVector contact_sph_bodies_;

		/** inner configuration for the neighbor relations. */
		ParticleConfiguration& inner_configuration_;
		/** Configurations for updated Lagrangian formulation. **/
		ContatcParticleConfiguration& contact_configuration_;

		SPHBodyComplexRelation(SPHBody* body, SPHBodyVector contact_sph_bodies);
		SPHBodyComplexRelation(SPHBodyInnerRelation* body_inner_relation, SPHBodyVector contact_sph_bodies);
		virtual ~SPHBodyComplexRelation() {
			delete inner_relation_;
			delete contact_relation_;
		};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration()  override;
	};
}
