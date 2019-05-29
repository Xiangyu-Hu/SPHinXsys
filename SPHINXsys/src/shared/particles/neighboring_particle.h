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
	 * @brief Friend class.
	 */
	class Kernel;
	/**
	 * @class NeighboringParticle
	 * @brief All the interacing particle j of particle i.
	 */
	class NeighboringParticle
	{
	public:
		size_t i_, j_;		/**< Index of the origina and neighbor particles. */
		Vecd r_ij_, e_ij_;	/**< Vecor and norm of particle postion,r_ij_ and e_ij_ pointing from i to j. */
		Real W_ij_, dW_ij_;	/**< Weighted fucntion and gardient of weighted fucntion. */
		Vecd gradW_ij_;		/**< Kernel gradient in weak form. */
		Vecd str_gradW_ij_;	/**< kernel gradient in strong form. */

		/**
		 * @brief Default construcutor. 
		 */
		NeighboringParticle() {};
		/**
		 * @brief Constructor.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		NeighboringParticle(Kernel &kernel,	Vecd &r_ij, size_t i_index, size_t j_index);
		/**
		 * @brief Destructor.
		 */
		~NeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t i_index, size_t j_index);
	};

	/**
	 * @class NeighboringParticleCofiguration 
	 * @brief interacing configuration between particle i and particle j.
	 */
	class NeighboringParticleCofiguration
	{
	public:
		Vecd r_ij_, e_ij_;	/**< Vecor and norm of particle postion,r_ij_ and e_ij_ pointing from i to j. */
		Real W_ij_, dW_ij_;	/**< Weighted fucntion and gardient of weighted fucntion. */
		Vecd gradW_ij_;		/**< Kernel gradient in weak form. */

		/**
		 * @brief Default construcutor.
		 */
		NeighboringParticleCofiguration() {};
		/**
		 * @brief Constructor.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 */
		NeighboringParticleCofiguration(Kernel &kernel, Vecd &r_ij);
		/**
		 * @brief Destructor.
		 */
		~NeighboringParticleCofiguration() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 */
		void Reset(Kernel &kernel, Vecd &r_ij);
	};

	/**
	  * @class ReferenceNeighboringParticle
	  * @brief All the interacing particle j of particle i from reference configuration.
	  */
	class ReferenceNeighboringParticle
	{
	public:
		size_t j_;		/**< Index of the origina and neighbor particles. */
		Vecd r_ij_;	/**< Vecor and norm of particle postion,r_ij_ and e_ij_ pointing from i to j. */
		Real W_ij_, dW_ij_;	/**< Weighted fucntion and gardient of weighted fucntion. */
		Vecd gradW_ij_;		/**< Kernel gradient in weak form. */

		/**
		 * @brief Default construcutor.
		 */
		ReferenceNeighboringParticle() {};
		/**
		 * @brief Constructor.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		ReferenceNeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t i_index, size_t j_index);
		/**
		 * @brief Destructor.
		 */
		~ReferenceNeighboringParticle() {};
		/**
		 * @brief Reset the neighboring particles.
		 * @param[in] kernel Specific kernel.
		 * @param[in] r_ij Postion vector of interacing particles.
		 * @param[in] i_index Index of particle i
		 * @param[in] j_index Index of particle j
		 */
		void Reset(Kernel &kernel, Vecd &r_ij, size_t i_index, size_t j_index);
	};
}