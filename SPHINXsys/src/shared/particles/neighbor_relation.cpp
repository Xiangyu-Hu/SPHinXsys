/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#include "neighbor_relation.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	CommonRelation::CommonRelation()
		: e_ij_(0), dW_ij_(0), r_ij_(0), j_(0)
	{

	}
	//=================================================================================================//
	CommonRelation::CommonRelation(Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
		: CommonRelation()
	{
		j_ = j_index;
		r_ij_ = vec_r_ij.norm();
		e_ij_ = vec_r_ij / (r_ij_ + TinyReal);
		dW_ij_ = kernel.dW(vec_r_ij);
	}
	//=================================================================================================//
}
//=================================================================================================//
