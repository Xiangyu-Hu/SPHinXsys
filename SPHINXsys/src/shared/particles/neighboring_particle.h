/**
 * @file 	neighboring_particle.h
 * @brief 	This is the classe neighboring partcle. It save the information for pair
 *			interaction, and also considered as the topology of the particles.
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
	 * @brief All the interacing particle j of particle i.
	 */
	class NeighboringParticle
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** Norm of particle postion,e_ij_ pointing from i to j. */
		Vecd e_ij_;
		/** Particle distance, Weighted fucntion and gardient of weighted fucntion. */
		Real r_ij_, W_ij_, dW_ij_;	

		/** default constrcutor*/
		NeighboringParticle() {};
		/**
		* @brief Constructor.
		* @param[in] kernel Specific kernel.
		* @param[in] r_ij Postion vector of interacing particles.
		* @param[in] i_index Index of particle i
		* @param[in] j_index Index of particle j
		*/
		NeighboringParticle(Kernel &kernel,	Vecd &r_ij, size_t j_index);
		~NeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t j_index);
	};

	/**
	  * @class ReferenceNeighboringParticle
	  * @brief All the interacing particle j of particle i from reference configuration.
	  */
	class ReferenceNeighboringParticle
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** Norm of particle postion,e_ij_ pointing from i to j. */
		Vecd e_ij_;	
		/** Particle distrane, Weighted fucntion and gardient of weighted fucntion. */
		Real r_ij_, W_ij_, dW_ij_;	

		/** default constrcutor*/
		ReferenceNeighboringParticle() {};
		/**
		 * @brief Constructor.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		ReferenceNeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index);
		~ReferenceNeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t j_index);
	};
}