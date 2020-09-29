/**
 * @file 	body_relation.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	hi ZHang and Xiangyu Hu
 * @version	0.1
 * 			0.2.0
 * 			Cell splitting algorithm are added.
 * 			Chi Zhang
 */
#include "base_kernel.h"
#include "body_relation.h"
#include "base_particles.h"

namespace SPH
{
	//=================================================================================================//
	SPHBodyBaseRelation::SPHBodyBaseRelation(SPHBody* sph_body)
		: sph_body_(sph_body), split_cell_lists_(sph_body->split_cell_lists_), base_particles_(sph_body->base_particles_),
		mesh_cell_linked_list_(sph_body->mesh_cell_linked_list_)
	{
	}
	//=================================================================================================//
	void SPHBodyBaseRelation::createNeighborRelation(Neighborhood& neighborhood,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		neighborhood.addANeighbor(kernel, vec_r_ij, i_index, j_index);
		neighborhood.memory_size_++;
	}
	//=================================================================================================//
	void SPHBodyBaseRelation::initializeNeighborRelation(Neighborhood& neighborhood, 
		size_t current_count_of_neighbors, Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		neighborhood.j_[current_count_of_neighbors] = j_index;
		neighborhood.W_ij_[current_count_of_neighbors] = kernel.W(vec_r_ij);
		neighborhood.dW_ij_[current_count_of_neighbors] = kernel.dW(vec_r_ij);
		Real r_ij = vec_r_ij.norm();
		neighborhood.r_ij_[current_count_of_neighbors] = r_ij;
		neighborhood.e_ij_[current_count_of_neighbors] = vec_r_ij / (r_ij + TinyReal);
	}
	//=================================================================================================//
	SPHBodyInnerRelation::SPHBodyInnerRelation(SPHBody* sph_body)
		: SPHBodyBaseRelation(sph_body)
	{
		subscribe_to_body();
		updateConfigurationMemories();
	};
	//=================================================================================================//
	void SPHBodyInnerRelation::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	SPHBodyContactRelation::SPHBodyContactRelation(SPHBody* sph_body, SPHBodyVector contact_sph_bodies)
		: SPHBodyBaseRelation(sph_body), contact_sph_bodies_(contact_sph_bodies) {
		for (size_t k = 0; k != contact_sph_bodies_.size(); ++k) {
			target_mesh_cell_linked_lists_.push_back(contact_sph_bodies_[k]->mesh_cell_linked_list_);
		}
		subscribe_to_body();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void SPHBodyContactRelation::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		contact_configuration_.resize(contact_sph_bodies_.size());
		for (size_t k = 0; k != contact_sph_bodies_.size(); ++k) {
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	SPHBodyComplexRelation::SPHBodyComplexRelation(SPHBody* body, SPHBodyVector contact_sph_bodies)
		: SPHBodyBaseRelation(body),
		inner_relation_(new SPHBodyInnerRelation(body)),
		contact_relation_(new SPHBodyContactRelation(body, contact_sph_bodies)),
		contact_sph_bodies_(contact_sph_bodies),
		inner_configuration_(inner_relation_->inner_configuration_),
		contact_configuration_(contact_relation_->contact_configuration_) 
	{
		updateConfigurationMemories();
	};
	//=================================================================================================//
	SPHBodyComplexRelation::SPHBodyComplexRelation(SPHBodyInnerRelation* body_inner_relation, 
		SPHBodyVector contact_sph_bodies) : SPHBodyBaseRelation(body_inner_relation->sph_body_),
		inner_relation_(body_inner_relation),
		contact_relation_(new SPHBodyContactRelation(sph_body_, contact_sph_bodies)),
		contact_sph_bodies_(contact_sph_bodies),
		inner_configuration_(inner_relation_->inner_configuration_),
		contact_configuration_(contact_relation_->contact_configuration_)
	{
		updateConfigurationMemories();
	};	
	//=================================================================================================//
	void SPHBodyComplexRelation::updateConfigurationMemories()
	{
		inner_relation_->updateConfigurationMemories();
		contact_relation_->updateConfigurationMemories();
	}
	//=================================================================================================//
	void SPHBodyComplexRelation::updateConfiguration()
	{
		inner_relation_->updateConfiguration();
		contact_relation_->updateConfiguration();
	}
	//=================================================================================================//
}