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
 * @file 	neighborhood.h
 * @brief 	There are the classes for particle neighborhood.
 * It saves the information for carrying out inter-particle (or pair) interaction,
 * and also considered as the local configuration of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "all_kernels.h"

namespace SPH
{

	class SPHBody;
	class BodyPart;

	/**
	 * @class Neighborhood
	 * @brief A neighborhood around particle i.
	 */
	class Neighborhood
	{
	public:
		size_t current_size_;	/**< the current number of neighbors */
		size_t allocated_size_; /**< the limit of neighbors does not require memory allocation  */

		StdLargeVec<size_t> j_;	  /**< index of the neighbor particle. */
		StdLargeVec<Real> W_ij_;  /**< kernel value or particle volume contribution */
		StdLargeVec<Real> dW_ijV_j_; /**< derivative of kernel function or inter-particle surface contribution */
		StdLargeVec<Real> r_ij_;  /**< distance between j and i. */
		StdLargeVec<Vecd> e_ij_;  /**< unit vector pointing from j to i or inter-particle surface direction */

		Neighborhood() : current_size_(0), allocated_size_(0){};
		~Neighborhood(){};

		void removeANeighbor(size_t neighbor_n);
	};

	/** Neighborhoods for all particles in a body for a inner body relation. */
	using ParticleConfiguration = StdLargeVec<Neighborhood>;
	/** Neighborhoods for all particles in a body for a contact body relation. */
	using ContactParticleConfiguration = StdVec<ParticleConfiguration>;

	/**
	 * @class NeighborBuilder
	 * @brief Base class for building neighbor between particles i and j.
	 */
	class NeighborBuilder
	{
	protected:
		Kernel *kernel_;
		//----------------------------------------------------------------------
		//	Below are for constant smoothing length.
		//----------------------------------------------------------------------
		void createNeighbor(Neighborhood &neighborhood, const Real &distance,
							const Vecd &displacement, size_t j_index, const Real Vol_j);
		void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
								const Vecd &displacement, size_t j_index, const Real Vol_j);
		//----------------------------------------------------------------------
		//	Below are for variable smoothing length.
		//----------------------------------------------------------------------
		void createNeighbor(Neighborhood &neighborhood, const Real &distance,
							const Vecd &displacement, size_t j_index, const Real Vol_j, Real i_h_ratio, Real h_ratio_min);
		void initializeNeighbor(Neighborhood &neighborhood, const Real &distance,
								const Vecd &displacement, size_t j_index, const Real Vol_j, Real i_h_ratio, Real h_ratio_min);

	public:
		NeighborBuilder() : kernel_(nullptr){};
		virtual ~NeighborBuilder(){};
	};

	/**
	 * @class NeighborBuilderInner
	 * @brief A inner neighbor relation functor between particles i and j.
	 */
	class NeighborBuilderInner : public NeighborBuilder
	{
	public:
		explicit NeighborBuilderInner(SPHBody &body);
		void operator()(Neighborhood &neighborhood,
						const Vecd &pos_i, size_t index_i, const ListData &list_data_j);
	};

	/**
	 * @class NeighborBuilderInnerVariableSmoothingLength
	 * @brief A inner neighbor relation functor between particles i and j
	 * when the particles have different smoothing lengths.
	 */
	class NeighborBuilderInnerVariableSmoothingLength : public NeighborBuilder
	{
	public:
		explicit NeighborBuilderInnerVariableSmoothingLength(SPHBody &body);
		void operator()(Neighborhood &neighborhood,
						const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

	protected:
		StdLargeVec<Real> &h_ratio_;
	};

	/**
	 * @class NeighborBuilderSelfContact
	 * @brief A self-contact neighbor builder between particles i and j.
	 */
	class NeighborBuilderSelfContact : public NeighborBuilder
	{
	public:
		explicit NeighborBuilderSelfContact(SPHBody &body);
		virtual ~NeighborBuilderSelfContact(){};
		void operator()(Neighborhood &neighborhood,
						const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

	protected:
		StdLargeVec<Vecd> &pos0_;
	};

	/**
	 * @class NeighborBuilderContact
	 * @brief A contact neighbor relation functor between particles i and j.
	 */
	class NeighborBuilderContact : public NeighborBuilder
	{
	public:
		NeighborBuilderContact(SPHBody &body, SPHBody &contact_body);
		virtual ~NeighborBuilderContact(){};
		void operator()(Neighborhood &neighborhood,
						const Vecd &pos_i, size_t index_i, const ListData &list_data_j);
	};

	/**
	 * @class NeighborBuilderSolidContact
	 * @brief A solid contact neighbor relation functor between particles i and j.
	 */
	class NeighborBuilderSolidContact : public NeighborBuilderContact
	{
	private:
		UniquePtrKeeper<Kernel> kernel_keeper_;

	public:
		NeighborBuilderSolidContact(SPHBody &body, SPHBody &contact_body);
		virtual ~NeighborBuilderSolidContact(){};
	};

	/**
	 * @class NeighborBuilderContactBodyPart
	 * @brief A contact neighbor relation functor between particles i and j.
	 */
	class NeighborBuilderContactBodyPart : public NeighborBuilder
	{
	public:
		NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part);
		virtual ~NeighborBuilderContactBodyPart(){};
		void operator()(Neighborhood &neighborhood,
						const Vecd &pos_i, size_t index_i, const ListData &list_data_j);

	protected:
		StdLargeVec<int> part_indicator_; /**< indicator of the body part */
	};
}
#endif // NEIGHBORHOOD_H