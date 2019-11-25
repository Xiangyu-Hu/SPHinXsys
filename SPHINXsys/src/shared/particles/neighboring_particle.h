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
		/** derivative of kernel function
		 *  and kernel fucntion values. */
		Real dW_ij_, W_ij_;
		/** Unit vector pointing from j to i. */
		Vecd e_ij_;
		/** Particle distance. */
		Real r_ij_;	

		/** default constrcutor*/
		NeighboringParticle() {};
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
	};

	/**
	  * @class ReferenceNeighboringParticle
	  * @brief A neigboring particle j of particle i from reference configuration.
	  */
	class ReferenceNeighboringParticle
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** derivative of kernel function
		 *  and kernel fucntion values. */
		Real dW_ij_, W_ij_;
		/** Unit vector pointing from j to i. */
		Vecd e_ij_;
		/** Particle distance. */
		Real r_ij_;

		/** default constrcutor*/
		ReferenceNeighboringParticle() {};
		/**
		 * @brief Constructor.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] j_index Index of particle j
		 */
		ReferenceNeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index);
		~ReferenceNeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t j_index);
	};
	/**
	  * @class ReferenceNeighboringParticleDiffusion
	  * @brief Interparticle averaged diffusion tensor
	  */
	class ReferenceNeighboringParticleDiffusion
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** Diffusion tensor */
		Matd diffusion_ij_;

		/** default constrcutor*/
		ReferenceNeighboringParticleDiffusion() {};
		/**
		 * @brief Constructor.
		 * @param[in] j_index Index of particle j
		 * @param[in] diff_ij Inter average diffusion tensor
		 */
		ReferenceNeighboringParticleDiffusion(size_t j_index, Matd &diff_ij);
		~ReferenceNeighboringParticleDiffusion() {};
		/**
		 * @brief Reset the neighboring particles diffusion tensor.
		 * @param[in] j_index Index of particle j
		 * @param[in] diff_ij Inter average diffusion tensor
		 */
		void Reset(size_t j_index, Matd &diff_ij);
	};
}