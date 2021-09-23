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
 * @file 	neighbor_relation.h
 * @brief 	There are the classes for neighboring particles. 
 * It saves the information for carring out pair
 * interaction, and also considered as the topology of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef NEIGHBOR_RELATION_H
#define NEIGHBOR_RELATION_H


#include "base_data_package.h"
#include "all_kernels.h"

namespace SPH {

	class SPHBody;
	class BodyPart;

	/**
	 * @class Neighborhood
	 * @brief A neighborhood around particle i.
	 */
	class Neighborhood
	{
	public:
		size_t current_size_; 		/**< the current number of neighors */
		size_t allocated_size_; 	/**< the limit of neighors does not require memory allocation  */

		StdLargeVec<size_t> j_;		/**< index of the neighbor particle. */
		StdLargeVec<Real> W_ij_;	/**< kernel value or particle volume contribution */
		StdLargeVec<Real> dW_ij_;	/**< derivative of kernel function or inter-particle surface contribution */
		StdLargeVec<Real> r_ij_;	/**< distance between j and i. */
		StdLargeVec<Vecd> e_ij_;	/**< unit vector pointing from j to i or inter-particle surface direction */

		Neighborhood() : current_size_(0), allocated_size_(0) {};
		~Neighborhood() {};

		void removeANeighbor(size_t neighbor_n);
	};

	/** Inner neighborhoods for all particles in a body for a inner body relation. */
	using ParticleConfiguration = StdLargeVec<Neighborhood>;
	/** All contact neighborhoods for all particles in a body for a contact body relation. */
	using ContatcParticleConfiguration = StdVec<ParticleConfiguration>;

	/**
	 * @class NeighborRelation
	 * @brief Base neighbor relation between particles i and j.
	 */
	class NeighborRelation
	{
	protected:
		Kernel* kernel_;
		//----------------------------------------------------------------------
		//	Below are for constant smoothing length.
		//----------------------------------------------------------------------
		void createRelation(Neighborhood& neighborhood,	Real& distance, 
			Vecd& displacement, size_t j_index) const;
		void initializeRelation(Neighborhood& neighborhood, Real& distance, 
			Vecd& displacement, size_t j_index) const;
		//----------------------------------------------------------------------
		//	Below are for variable smoothing length.
		//----------------------------------------------------------------------
		void createRelation(Neighborhood& neighborhood, Real& distance, 
			Vecd& displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const;
		void initializeRelation(Neighborhood& neighborhood, Real& distance, 
			Vecd& displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const;
	public:
		NeighborRelation() : kernel_(nullptr) {};
		virtual ~NeighborRelation() {};
	};

	/**
	 * @class NeighborRelationInner
	 * @brief A inner neighbor relation functor between particles i and j.
	 */
	class NeighborRelationInner : public NeighborRelation
	{
	public:
		NeighborRelationInner(SPHBody* body);
		void operator () (Neighborhood& neighborhood, 
			Vecd& displacement, size_t i_index, size_t j_index) const;
	};

	/**
	 * @class NeighborRelationInnerVariableSmoothingLength
	 * @brief A inner neighbor relation functor between particles i and j.
	 */
	class NeighborRelationInnerVariableSmoothingLength : public NeighborRelation
	{
	public:
		NeighborRelationInnerVariableSmoothingLength(SPHBody* body);
		void operator () (Neighborhood& neighborhood, 
			Vecd& displacement, size_t i_index, size_t j_index) const;
	protected:
		StdLargeVec<Real>& h_ratio_;
	};
	
	/**
	 * @class NeighborRelationSelfContact
	 * @brief A solid contact neighbor relation functor between particles i and j.
	 */
	class NeighborRelationSelfContact : public NeighborRelation
	{
	public:
		explicit NeighborRelationSelfContact(SPHBody* body);
		virtual ~NeighborRelationSelfContact() {};
		void operator () (Neighborhood& neighborhood, 
			Vecd& displacement, size_t i_index, size_t j_index) const;
	protected:
		StdLargeVec<Vecd>& pos_0_;
	};

	/**
	 * @class NeighborRelationContact
	 * @brief A contact neighbor relation functor between particles i and j.
	 */
	class NeighborRelationContact : public NeighborRelation
	{
	public:
		NeighborRelationContact(SPHBody* body, SPHBody* contact_body);
		virtual ~NeighborRelationContact() {};
		void operator () (Neighborhood& neighborhood,
			Vecd& displacement, size_t i_index, size_t j_index) const;
	};

	/**
	 * @class NeighborRelationSolidContact
	 * @brief A solid contact neighbor relation functor between particles i and j.
	 */
	class NeighborRelationSolidContact : public NeighborRelationContact
	{
	public:
		NeighborRelationSolidContact(SPHBody* body, SPHBody* contact_body);
		virtual ~NeighborRelationSolidContact() {};
	};

	/**
	 * @class NeighborRelationContactBodyPart
	 * @brief A contact neighbor relation functor between particles i and j.
	 */
	class NeighborRelationContactBodyPart : public NeighborRelation
	{
	public:
		NeighborRelationContactBodyPart(SPHBody* body, BodyPart* contact_body_part);
		virtual ~NeighborRelationContactBodyPart() {};
		void operator () (Neighborhood& neighborhood,
			Vecd& displacement, size_t i_index, size_t j_index) const;
	protected:
		StdLargeVec<int>& part_indicator_;	/**< indicator of the body part */
	};	
}
#endif //NEIGHBOR_RELATION_H