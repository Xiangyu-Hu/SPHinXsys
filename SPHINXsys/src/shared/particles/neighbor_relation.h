/**
 * @file 	neighbor_relation.h
 * @brief 	There are the classes for a neighboring partcle. 
 * It saves the information for carring out pair
 * interaction, and also considered as the topology of the particles.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once
#include "base_data_package.h"

using namespace std;

namespace SPH {
	/**
	 * @brief preclaimed class.
	 */
	class Kernel;
	class BaseParticleData;

	/**
	 * @class NeighborRelation
	 * @brief A neighboring realtion for a particle j around particle i.
	 */
	class BaseNeighborRelation
	{
	public:
		/** Index of the neighbor particle. */
		size_t j_;
		/** kernel function value. */
		Real W_ij_;
		/** Derivative of kernel function. */
		Real dW_ij_;
		/** Unit vector pointing from j to i. */
		Vecd e_ij_;
		/** Distance between i and j. */
		Real r_ij_;

		/** default constructor*/
		BaseNeighborRelation() {};
		/**
		* @brief Constructor.
		* @param[in] base_particles Particles with geometric informaiton.
		* @param[in] kernel Specific kernel.
		* @param[in] r_ij Postion vector pointing from j to i.
		* @param[in] i_index Index of particle i
		* @param[in] j_index Index of particle j
		*/
		BaseNeighborRelation(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
		~BaseNeighborRelation() {};
		/**
		 * @brief Reset the neighboring particles.
		* @param[in] base_particles Particles with geometric informaiton.
		* @param[in] kernel Specific kernel.
		* @param[in] r_ij Postion vector pointing from j to i.
		* @param[in] i_index Index of particle i
		* @param[in] j_index Index of particle j
		*/
		virtual void resetRelation(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index) = 0;
		/** reset from a symmetric realtion. */
		void resetSymmetricRelation(BaseNeighborRelation* symmetic_neigbor_relation, size_t i_index);

		/** compute gradient of the kernel function. */
		virtual Vecd getNablaWij() { return dW_ij_ * e_ij_; };
	};

	/**
	 * @class NeighborRelation
	 * @brief A neighboring realtion for a particle j around particle i.
	 */
	class NeighborRelation : public BaseNeighborRelation
	{
	public:
		/** Constructor. */
		NeighborRelation(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
		/** Constructor with symmetic condition.*/
		NeighborRelation(BaseNeighborRelation* symmetic_neigbor_relation, size_t i_index);
		~NeighborRelation() {};

		/** Reset the neighboring particles. */
		virtual void resetRelation(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index) override;
	};

	/**
	 * @class NeighborRelationWithVariableSmoothingLength
	 * @brief A neighboring particle j of particle i.
	 */
	class NeighborRelationWithVariableSmoothingLength : public BaseNeighborRelation
	{
	public:
		/** Constructor. */
		NeighborRelationWithVariableSmoothingLength(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
		/** Constructor with symmetic condition.*/
		NeighborRelationWithVariableSmoothingLength(BaseNeighborRelation* symmetic_neigbor_relation, size_t i_index);
		~NeighborRelationWithVariableSmoothingLength() {};

		/** Reset the neighboring particles. */
		virtual void resetRelation(StdLargeVec<BaseParticleData>& base_particle_data,
			Kernel& kernel, Vecd& r_ij, size_t i_index, size_t j_index) override;
	};
}