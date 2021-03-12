/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "neighbor_relation.h"

namespace SPH
{
	//=================================================================================================//
	NeighborRelation::NeighborRelation() :
		kernel_(NULL), cutoff_radius_(0.0) {}
	//=================================================================================================//
	void NeighborRelation::createRelation(Neighborhood& neighborhood,
		Kernel* kernel, Real& distance, Vecd& displacement, size_t j_index) const
	{
		neighborhood.j_.push_back(j_index);
		neighborhood.W_ij_.push_back(kernel->W(distance, displacement));
		neighborhood.dW_ij_.push_back(kernel->dW(distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborRelation::initializeRelation(Neighborhood& neighborhood, 
		Kernel* kernel, Real& distance, Vecd& displacement, size_t j_index) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = kernel->W(distance, displacement);
		neighborhood.dW_ij_[current_size] = kernel->dW(distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
	//=================================================================================================//
	NeighborRelationVariableSmoothingLength::NeighborRelationVariableSmoothingLength() :
		kernel_(NULL) {}
	//=================================================================================================//
	void NeighborRelationVariableSmoothingLength::
		createRelation(Neighborhood& neighborhood, Kernel* kernel, Real& distance, 
			Vecd& displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		neighborhood.j_.push_back(j_index);
		Real weight = distance < kernel->CutOffRadius(i_h_ratio) ?
			kernel->W(i_h_ratio, distance, displacement) : 0.0;
		neighborhood.W_ij_.push_back(weight);
		neighborhood.dW_ij_.push_back(kernel->dW(h_ratio_min, distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborRelationVariableSmoothingLength::
		initializeRelation(Neighborhood& neighborhood, Kernel* kernel, Real& distance, 
			Vecd& displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = distance < kernel->CutOffRadius(i_h_ratio) ?
			kernel->W(i_h_ratio, distance, displacement) : 0.0;
		neighborhood.dW_ij_[current_size] = kernel->dW(h_ratio_min, distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
	//=================================================================================================//
	NeighborRelationInner::NeighborRelationInner(SPHBody* body) :
		NeighborRelation() 
	{
		kernel_ = body->particle_adaptation_->getKernel();
		cutoff_radius_ =kernel_->CutOffRadius();
	}
	//=================================================================================================//
	void NeighborRelationInner::operator () (Neighborhood& neighborhood,
		Vecd& displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < cutoff_radius_ && i_index != j_index)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_ ?
				createRelation(neighborhood, kernel_, distance, displacement, j_index)
				: initializeRelation(neighborhood, kernel_, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborRelationInnerVariableSmoothingLength::
		NeighborRelationInnerVariableSmoothingLength(SPHBody* body) :
		NeighborRelationVariableSmoothingLength(),
		h_ratio_(body->base_particles_->h_ratio_)
	{
		kernel_ = body->particle_adaptation_->getKernel();
	}
	//=================================================================================================//
	void NeighborRelationInnerVariableSmoothingLength::operator () (Neighborhood& neighborhood, 
		Vecd& displacement, size_t i_index, size_t j_index) const
	{
		Real i_h_ratio = h_ratio_[i_index];
		Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[j_index]);
		Real cutoff_radius = kernel_->CutOffRadius(h_ratio_min);
		Real distance = displacement.norm();
		if (distance < cutoff_radius && i_index != j_index)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_ ?
				createRelation(neighborhood, kernel_, distance, displacement, j_index, i_h_ratio, h_ratio_min)
				: initializeRelation(neighborhood, kernel_, distance, displacement, j_index, i_h_ratio, h_ratio_min);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborRelationContact::NeighborRelationContact(SPHBody* body, SPHBody* contact_body) :
		NeighborRelation()
	{
		Kernel* source_kernel = body->particle_adaptation_->getKernel();
		Kernel* target_kernel = contact_body->particle_adaptation_->getKernel();
		kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
		cutoff_radius_ = kernel_->CutOffRadius();
	}
	//=================================================================================================//
	void NeighborRelationContact::operator () (Neighborhood& neighborhood,
		Vecd& displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < cutoff_radius_)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_ ?
				createRelation(neighborhood, kernel_, distance, displacement, j_index)
				: initializeRelation(neighborhood, kernel_, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
}
//=================================================================================================//
