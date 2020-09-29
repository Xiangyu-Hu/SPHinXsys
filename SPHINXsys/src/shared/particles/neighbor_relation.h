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
 * @brief 	There are the classes for a neighboring particle. 
 * It saves the information for carring out pair
 * interaction, and also considered as the topology of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once
#include "base_data_package.h"
#include "base_kernel.h"

using namespace std;

namespace SPH {
	/**
	 * @class Neighborhood
	 * @brief A neighborhood around particle i.
	 */
	class Neighborhood
	{
	public:
		/** the currnet number of neighors */
		size_t current_size_;
		/** the limit of neighors does not require memory allocation  */
		size_t memory_size_;

		StdLargeVec<size_t> j_;		/**< index of the neighbor particle. */
		StdLargeVec<Real> W_ij_;	/**< kernel value */
		StdLargeVec<Real> dW_ij_;	/**< derivative of kernel function */
		StdLargeVec<Real> r_ij_;	/**< distance between j and i. */
		StdLargeVec<Vecd> e_ij_;	/**< unit vector pointing from j to i */

		/** default constructor */
		Neighborhood()
			: current_size_(0), memory_size_(0) {};
		~Neighborhood() {};

		void addANeighbor(Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
	};

	/** A neighborhoods for all particles in a body. */
	using ParticleConfiguration = StdLargeVec<Neighborhood>;
	/** All contact neighborhoods for all particles in a body. */
	using ContatcParticleConfiguration = StdVec<ParticleConfiguration>;
}
