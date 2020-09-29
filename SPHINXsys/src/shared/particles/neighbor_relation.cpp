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
	void Neighborhood::addANeighbor(Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		j_.push_back(j_index);
		W_ij_.push_back(kernel.W(vec_r_ij));
		dW_ij_.push_back(kernel.dW(vec_r_ij));
		Real r_ij = vec_r_ij.norm();
		r_ij_.push_back(r_ij);
		e_ij_.push_back(vec_r_ij / (r_ij + TinyReal));
	}
	//=================================================================================================//
}
//=================================================================================================//
