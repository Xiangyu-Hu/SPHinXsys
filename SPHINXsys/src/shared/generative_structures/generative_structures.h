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
* @file generative_structures.h
* @brief In generative structures, the particles and, if necessary, 
* and their neighnor relations or particle configuration
* can be generated directly form the information of the structure.   
* @author	Xiangyu Hu
*/

#ifndef GENERATIVE_STRUCTURES_H
#define GENERATIVE_STRUCTURES_H

#include "base_geometry.h"
#include "neighbor_relation.h"

namespace SPH
{

	class SPHBody;
	class BaseParticles;

	/**
	 * @class GenerativeStructure
	 * @brief Abstract class as interface for all generative structures.
	 * It is linked with particles by particle position and volume
	 */
	class GenerativeStructure
	{
	public:
		explicit GenerativeStructure(SPHBody *sph_body);
		virtual ~GenerativeStructure(){};

		virtual void buildParticleConfiguration(BaseParticles &base_particles,
												ParticleConfiguration &particle_configuration) = 0;

	protected:
		SPHBody *sph_body_;
		Real spacing_ref_;
		BaseParticles *base_particles_;
		NeighborRelationInner neighbor_relation_inner_;

	public:
		StdLargeVec<Vecd> &pos_n_; /**< current position */
		StdLargeVec<Real> &Vol_;   /**< particle volume */
	};

	/**
	 * @class GenerativeTree
	 * @brief The tree is composed of a root (the first branch) 
	 * and other branch generated sequentially.
	 */
	class GenerativeTree : public GenerativeStructure
	{
	public:
		class Branch;

	private:
		UniquePtrVectorKeeper<Branch> branches_ptr_keeper_;

	public:
		StdVec<Branch *> branches_;	   /**< Contanier of all branches */
		IndexVector branch_locations_; /**< in which branch are the particles located */
		size_t last_branch_id_;
		Branch *root_;

		explicit GenerativeTree(SPHBody *sph_body);
		virtual ~GenerativeTree(){};

		Branch *createANewBranch(size_t parent_id)
		{
			return branches_ptr_keeper_.createPtr<Branch>(parent_id, this);
		};
		void growAParticleOnBranch(Branch *branch, const Vecd &new_point, const Vecd &end_direction);
		size_t BranchLocation(size_t particle_idx);
		Branch *LastBranch() { return branches_[last_branch_id_]; };

		virtual void buildParticleConfiguration(BaseParticles &base_particles,
												ParticleConfiguration &particle_configuration) override;
		size_t ContainerSize() { return branches_.size(); };
	};

	/**
	 * @class GenerativeTree::Branch
	 * @brief Each branch (excapt the root) has a parent and several children, and geometric information.
	 * It is a realized edge and has multi inner particles.
	 * The first is the last particle from the parent or root,
	 * and the last is the first particle of all its child branches.
	 * Many connected branches compose a tree. */
	class GenerativeTree::Branch : public Edge<size_t, IndexVector>
	{
	public:
		/** construct the root branch  */
		explicit Branch(GenerativeTree *tree);
		/** construct an branch connecting with its parent */
		Branch(size_t parent_id, GenerativeTree *tree);
		virtual ~Branch(){};

		Vecd end_direction_; /**< the direction pointing to the last particle */
		/** The indexes of particle within this branch.
		 * The first is the last particle from the parent or root,
		 * and the last is the first of all its child branches. */
		IndexVector inner_particles_;
		bool is_terminated_; /**< whether is an terminate branch or not */
	};

}
#endif //GENERATIVE_STRUCTURES_H