/**
 * @file 	neighboring_particle.h
 * @brief 	There are the classes for a neighboring partcle. 
 * It saves the information for carring out pair
 * interaction, and also considered as the topology of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once
#include "base_data_package.h"

namespace SPH {
	/**
	 * @brief preclaimed class.
	 */
	class Kernel;
	/**
	 * @class NeighboringParticle
	 * @brief A neighboring particle j of particle i.
	 */
	class NeighboringParticle
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** Derivative of kernel function and kernel fucntion values. */
		Real dW_ij_, W_ij_;
		/** Unit vector pointing from j to i. */
		Vecd e_ij_;
		/** Distance between i and j. */
		Real r_ij_;	

		/** default constrcutor*/
		NeighboringParticle() : j_(0), dW_ij_(0.0), W_ij_(0.0), 
			e_ij_(FisrtAxisVector(Vecd(0))), r_ij_(1.0) {};
		/**
		* @brief Constructor.
		* @param[in] kernel Specific kernel.
		* @param[in] r_ij Postion vector of interacing particles.
		* @param[in] j_index Index of particle j
		*/
		NeighboringParticle(Kernel &kernel,	Vecd &r_ij, size_t j_index);
		~NeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t j_index);

		/** compute gradient of the kernel function. */
		Vecd getNablaWij() { return dW_ij_ * e_ij_; };
	};
}