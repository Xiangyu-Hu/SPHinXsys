/**
 * @file 	body_relation.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "body_relation.h"

#include "base_kernel.h"
#include "base_particles.h"
#include "mesh_cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	SPHBodyRelation::SPHBodyRelation(SPHBody* sph_body)
		: sph_body_(sph_body), base_particles_(sph_body->base_particles_) {}
	//=================================================================================================//
	BaseInnerBodyRelation::BaseInnerBodyRelation(RealBody* real_body)
		: SPHBodyRelation(real_body), real_body_(real_body)
	{
		subscribe_to_body();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseInnerBodyRelation::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	void  BaseInnerBodyRelation::resetNeighborhoodCurrentSize()
	{
		parallel_for(blocked_range<size_t>(0, base_particles_->total_real_particles_),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {
					inner_configuration_[num].current_size_ = 0;
				}
			}, ap);
	}
	//=================================================================================================//
	InnerBodyRelation::InnerBodyRelation(RealBody* real_body)
		: BaseInnerBodyRelation(real_body), get_inner_neighbor_(real_body),
		mesh_cell_linked_list_(dynamic_cast<MeshCellLinkedList*>(real_body->mesh_cell_linked_list_)) {}
	//=================================================================================================//
	void InnerBodyRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		mesh_cell_linked_list_->searchNeighborsByParticles(base_particles_->total_real_particles_, *base_particles_,
			inner_configuration_, get_particle_index_, get_single_search_range_, get_inner_neighbor_);
	}
	//=================================================================================================//
	InnerBodyRelationVariableSmoothingLength::
		InnerBodyRelationVariableSmoothingLength(RealBody* real_body)
		: BaseInnerBodyRelation(real_body), total_levels_(0),
		get_inner_neighbor_variable_smoothing_length_(real_body)
	{
		MultilevelMeshCellLinkedList* multi_level_mesh_cell_linked_list =
			dynamic_cast<MultilevelMeshCellLinkedList*>(real_body->mesh_cell_linked_list_);
		mesh_cell_linked_list_levels_ = multi_level_mesh_cell_linked_list->getMeshLevels();
		total_levels_ = mesh_cell_linked_list_levels_.size();
		for (size_t l = 0; l != total_levels_; ++l) {
			get_multi_level_search_range_.push_back(
				new SearchRangeVariableSmoothingLength(real_body, mesh_cell_linked_list_levels_[l]));
		}
	}
	//=================================================================================================//
	void InnerBodyRelationVariableSmoothingLength::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		for (size_t l = 0; l != total_levels_; ++l) {
			mesh_cell_linked_list_levels_[l]->searchNeighborsByParticles(base_particles_->total_real_particles_,
				*base_particles_, inner_configuration_, get_particle_index_,
				*get_multi_level_search_range_[l], get_inner_neighbor_variable_smoothing_length_);
		}
	}	
	//=================================================================================================//
	BaseContactBodyRelation::BaseContactBodyRelation(SPHBody* sph_body, RealBodyVector contact_sph_bodies)
		: SPHBodyRelation(sph_body), contact_bodies_(contact_sph_bodies)
	{
		subscribe_to_body();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseContactBodyRelation::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		contact_configuration_.resize(contact_bodies_.size());
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	void  BaseContactBodyRelation::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			parallel_for(blocked_range<size_t>(0, base_particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
					for (size_t num = r.begin(); num != r.end(); ++num) {
						contact_configuration_[k][num].current_size_ = 0;
					}
				}, ap);
		}
	}
	//=================================================================================================//
	ContactBodyRelation::ContactBodyRelation(SPHBody* sph_body, RealBodyVector contact_sph_bodies)
		: BaseContactBodyRelation(sph_body, contact_sph_bodies)
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			target_mesh_cell_linked_lists_.push_back(dynamic_cast<MeshCellLinkedList*>(contact_bodies_[k]->mesh_cell_linked_list_));
			get_search_ranges_.push_back(new SearchRangeMultiResolution(sph_body, contact_sph_bodies[k]));
			get_contact_neighbors_.push_back(new NeighborRelationContact(sph_body, contact_sph_bodies[k]));
		}
	}
	//=================================================================================================//
	void ContactBodyRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			target_mesh_cell_linked_lists_[k]->searchNeighborsByParticles(total_real_particles, 
				*base_particles_, contact_configuration_[k],
				get_particle_index_, *get_search_ranges_[k], *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	SolidContactBodyRelation::SolidContactBodyRelation(SPHBody* sph_body, RealBodyVector contact_sph_bodies)
		: ContactBodyRelation(sph_body, contact_sph_bodies), 
		body_part_particles_(body_surface_layer_.body_part_particles_),
		get_body_part_particle_index_(body_part_particles_),
		body_surface_layer_(ShapeSurfaceLayer(sph_body)) {}
	//=================================================================================================//
	void  SolidContactBodyRelation::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			parallel_for(blocked_range<size_t>(0, body_part_particles_.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t num = r.begin(); num != r.end(); ++num) {
						size_t index_i = get_body_part_particle_index_(num);
						contact_configuration_[k][index_i].current_size_ = 0;
					}
				}, ap);
		}
	}
	//=================================================================================================//
	void SolidContactBodyRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = body_part_particles_.size();
		for (size_t k = 0; k != contact_bodies_.size(); ++k) {
			target_mesh_cell_linked_lists_[k]->searchNeighborsByParticles(total_real_particles, 
				*base_particles_, contact_configuration_[k],
				get_body_part_particle_index_,*get_search_ranges_[k], *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(BaseInnerBodyRelation* inner_relation, BaseContactBodyRelation* contact_relation) : 
		SPHBodyRelation(inner_relation->sph_body_),
		inner_relation_(inner_relation), contact_relation_(contact_relation),
		contact_bodies_(contact_relation->contact_bodies_),
		inner_configuration_(inner_relation->inner_configuration_),
		contact_configuration_(contact_relation->contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(RealBody* real_body, RealBodyVector contact_bodies) :
		ComplexBodyRelation(new InnerBodyRelation(real_body), new ContactBodyRelation(real_body, contact_bodies)) {}
	//=================================================================================================//
	ComplexBodyRelation::
		ComplexBodyRelation(BaseInnerBodyRelation* inner_relation, RealBodyVector contact_bodies) :
		ComplexBodyRelation(inner_relation, new ContactBodyRelation(inner_relation->sph_body_, contact_bodies)) {}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfigurationMemories()
	{
		inner_relation_->updateConfigurationMemories();
		contact_relation_->updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfiguration()
	{
		inner_relation_->updateConfiguration();
		contact_relation_->updateConfiguration();
	}
	//=================================================================================================//
}
