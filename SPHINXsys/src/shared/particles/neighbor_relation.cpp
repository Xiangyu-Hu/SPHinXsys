/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#include "neighbor_relation.h"
#include "base_particles.h"
#include "base_kernel.h"
//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	BaseNeighborRelation::BaseNeighborRelation(StdLargeVec<BaseParticleData>& base_particle_data,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
		: j_(j_index), e_ij_(normalize(vec_r_ij)), r_ij_(vec_r_ij.norm()), W_ij_(0), dW_ij_(0) {}
	//=================================================================================================//
	void BaseNeighborRelation::resetSymmetricRelation(BaseNeighborRelation*
		symmetic_neigbor_relation, size_t i_index)
	{
		j_ = i_index;
		e_ij_ = -symmetic_neigbor_relation->e_ij_;
		r_ij_ = symmetic_neigbor_relation->r_ij_;
		W_ij_ = symmetic_neigbor_relation->W_ij_;
		dW_ij_ = symmetic_neigbor_relation->dW_ij_;
	}
	//=================================================================================================//
	NeighborRelation::NeighborRelation(StdLargeVec<BaseParticleData>& base_particle_data,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
		: BaseNeighborRelation(base_particle_data, kernel, vec_r_ij, i_index, j_index)
	{
		W_ij_ = kernel.W(vec_r_ij);
		dW_ij_ = kernel.dW(vec_r_ij);
	}
	//=================================================================================================//
	NeighborRelation::NeighborRelation(BaseNeighborRelation* symmetic_neigbor_relation, size_t i_index)
		: BaseNeighborRelation()
	{
		resetSymmetricRelation(symmetic_neigbor_relation, i_index);
	}
	//=================================================================================================//
	void NeighborRelation::resetRelation(StdLargeVec<BaseParticleData>& base_particle_data,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		j_ 		= j_index;
		e_ij_ 	= normalize(vec_r_ij);
		r_ij_ 	= vec_r_ij.norm();
		W_ij_ 	= kernel.W(vec_r_ij);
		dW_ij_ 	= kernel.dW(vec_r_ij);
	}
	//=================================================================================================//
	NeighborRelationWithVariableSmoothingLength
		::NeighborRelationWithVariableSmoothingLength(StdLargeVec<BaseParticleData>& base_particle_data,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
		: BaseNeighborRelation(base_particle_data, kernel, vec_r_ij, i_index, j_index)
	{
		Real inv_smoothing_length
			= 1.0 / SMAX(base_particle_data[i_index].smoothing_length_, 
				base_particle_data[j_index].smoothing_length_);
		W_ij_ = kernel.W(inv_smoothing_length, vec_r_ij);
		dW_ij_ = kernel.dW(inv_smoothing_length, vec_r_ij);
	}
	//=================================================================================================//
	void NeighborRelationWithVariableSmoothingLength
		::resetRelation(StdLargeVec<BaseParticleData>& base_particle_data,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		j_ = j_index;
		e_ij_ = normalize(vec_r_ij);
		r_ij_ = vec_r_ij.norm();
		Real inv_smoothing_length
			= 1.0 / SMAX(base_particle_data[i_index].smoothing_length_, 
				base_particle_data[j_index].smoothing_length_);
		W_ij_ = kernel.W(inv_smoothing_length, vec_r_ij);
		dW_ij_ = kernel.dW(inv_smoothing_length, vec_r_ij);
	}
	//=================================================================================================//
	NeighborRelationWithVariableSmoothingLength
		::NeighborRelationWithVariableSmoothingLength(BaseNeighborRelation* symmetic_neigbor_relation, size_t i_index)
		: BaseNeighborRelation()
	{
		resetSymmetricRelation(symmetic_neigbor_relation, i_index);
	}
	//=================================================================================================//
}
//=================================================================================================//
