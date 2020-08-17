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
	SPHBodyBaseRelation::SPHBodyBaseRelation(SPHBody* body)
		: body_(body), split_cell_lists_(body->split_cell_lists_), base_particles_(body->base_particles_),
		base_mesh_cell_linked_list_(body->base_mesh_cell_linked_list_) 
	{
	}
	//=================================================================================================//
	void SPHBodyBaseRelation::createNeighborRelation(Neighborhood& neighborhood,
		Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		neighborhood.common_relation_list_.push_back(CommonRelation(kernel, vec_r_ij, i_index, j_index));
		neighborhood.kernel_value_list_.push_back(kernel.W(vec_r_ij));
		neighborhood.memory_size_++;
	}
	//=================================================================================================//
	void SPHBodyBaseRelation::initializeNeighborRelation(Neighborhood& neighborhood, 
		size_t current_count_of_neighbors, Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index)
	{
		CommonRelation& common_relation = neighborhood.common_relation_list_[current_count_of_neighbors];
		common_relation.j_ = j_index;
		common_relation.r_ij_ = vec_r_ij.norm();
		common_relation.e_ij_ = vec_r_ij / (common_relation.r_ij_ + TinyReal);
		common_relation.dW_ij_ = kernel.dW(vec_r_ij);
		neighborhood.kernel_value_list_[current_count_of_neighbors] = kernel.W(vec_r_ij);
	}
	//=================================================================================================//
	SPHBodyInnerRelation::SPHBodyInnerRelation(SPHBody* body)
		: SPHBodyBaseRelation(body)
	{
		subscribe_to_body();
		updateConfigurationMemories();
	};
	//=================================================================================================//
	void SPHBodyInnerRelation::updateConfigurationMemories()
	{
		size_t updated_size = body_->base_particles_->real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	SPHBodyContactRelation::SPHBodyContactRelation(SPHBody* body, SPHBodyVector relation_bodies)
		: SPHBodyBaseRelation(body), relation_bodies_(relation_bodies) {
		for (size_t k = 0; k != relation_bodies_.size(); ++k) {
			target_mesh_cell_linked_lists_.push_back(relation_bodies_[k]->base_mesh_cell_linked_list_);
		}
		subscribe_to_body();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void SPHBodyContactRelation::updateConfigurationMemories()
	{
		size_t updated_size = body_->base_particles_->real_particles_bound_;
		contact_configuration_.resize(relation_bodies_.size());
		for (size_t k = 0; k != relation_bodies_.size(); ++k) {
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	bool SPHBodyContactRelation::checkNeighbor(Real particle_distance_sqr, Real cutoff_radius_sqr,
		BaseParticleData& base_particle_data_i, BaseParticleData& base_particle_data_j) 
	{
		return particle_distance_sqr < cutoff_radius_sqr ? true : false;
	}
	//=================================================================================================//
	bool SPHBodyCollisionRelation::checkNeighbor(Real particle_distance, Real cutoff_radius,
		BaseParticleData& base_particle_data_i, BaseParticleData& base_particle_data_j) 
	{
		return particle_distance < cutoff_radius
			&& (base_particle_data_i.pos_0_ - base_particle_data_j.pos_0_).norm() > cutoff_radius
			? true : false;
	}
	//=================================================================================================//
	SPHBodyComplexRelation::SPHBodyComplexRelation(SPHBody* body, SPHBodyVector contact_bodies)
		: SPHBodyBaseRelation(body),
		inner_relation_(new SPHBodyInnerRelation(body)),
		inner_configuration_(inner_relation_->inner_configuration_),
		contact_relation_(new SPHBodyContactRelation(body, contact_bodies)),
		relation_bodies_(contact_bodies),
		contact_configuration_(contact_relation_->contact_configuration_) 
	{
		updateConfigurationMemories();
	};
	//=================================================================================================//
	SPHBodyComplexRelation::SPHBodyComplexRelation(SPHBodyInnerRelation* body_inner_relation, 
		SPHBodyVector contact_bodies) : SPHBodyBaseRelation(body_inner_relation->body_),
		inner_relation_(body_inner_relation),
		inner_configuration_(inner_relation_->inner_configuration_),
		contact_relation_(new SPHBodyContactRelation(body_, contact_bodies)),
		relation_bodies_(contact_bodies),
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