/**
 * @file 	body_relation.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "contact_body_relation.h"

#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	BodyRelationContact::BodyRelationContact(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: BaseBodyRelationContact(sph_body, contact_sph_bodies)
	{
		initialization();
	}
	//=================================================================================================//
	BodyRelationContact::BodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts)
		: BaseBodyRelationContact(sph_body, contact_body_parts)
	{
		initialization();
	}
	//=================================================================================================//
	void BodyRelationContact::initialization()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			CellLinkedList *target_cell_linked_list =
				DynamicCast<CellLinkedList>(this, contact_bodies_[k]->cell_linked_list_);
			target_cell_linked_lists_.push_back(target_cell_linked_list);
			get_search_depths_.push_back(
				search_depth_multi_resolution_ptr_vector_keeper_.createPtr<SearchDepthMultiResolution>(*sph_body_, target_cell_linked_list));
			get_contact_neighbors_.push_back(
				neighbor_relation_contact_ptr_vector_keeper_.createPtr<NeighborRelationContact>(sph_body_, contact_bodies_[k]));
		}
	}
	//=================================================================================================//
	void BodyRelationContact::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]
				->searchNeighborsByParticles(total_real_particles,
											 *base_particles_, contact_configuration_[k],
											 get_particle_index_, *get_search_depths_[k],
											 *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	SolidBodyRelationContact::SolidBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies)
		: BaseBodyRelationContact(sph_body, contact_bodies),
		  body_surface_layer_(shape_surface_ptr_keeper_.createPtr<BodySurfaceLayer>(sph_body)),
		  body_part_particles_(body_surface_layer_->body_part_particles_),
		  get_body_part_particle_index_(body_part_particles_)
	{
		initialization();
	}
	//=================================================================================================//
	SolidBodyRelationContact::
		SolidBodyRelationContact(SolidBodyRelationSelfContact &solid_body_relation_self_contact,
								 RealBodyVector contact_bodies)
		: BaseBodyRelationContact(*solid_body_relation_self_contact.real_body_, contact_bodies),
		  body_surface_layer_(&solid_body_relation_self_contact.body_surface_layer_),
		  body_part_particles_(body_surface_layer_->body_part_particles_),
		  get_body_part_particle_index_(body_part_particles_)
	{
		initialization();
	}
	//=================================================================================================//
	void SolidBodyRelationContact::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			parallel_for(
				blocked_range<size_t>(0, body_part_particles_.size()),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t num = r.begin(); num != r.end(); ++num)
					{
						size_t index_i = get_body_part_particle_index_(num);
						contact_configuration_[k][index_i].current_size_ = 0;
					}
				},
				ap);
		}
	}
	//=================================================================================================//
	void SolidBodyRelationContact::initialization()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			CellLinkedList *target_cell_linked_list =
				DynamicCast<CellLinkedList>(this, contact_bodies_[k]->cell_linked_list_);
			target_cell_linked_lists_.push_back(target_cell_linked_list);
			get_search_depths_.push_back(
				search_depth_multi_resolution_ptr_vector_keeper_.createPtr<SearchDepthMultiResolution>(*sph_body_, target_cell_linked_list));
			get_contact_neighbors_.push_back(
				neighbor_relation_contact_ptr_vector_keeper_.createPtr<NeighborRelationSolidContact>(sph_body_, contact_bodies_[k]));
		}
	}
	//=================================================================================================//
	void SolidBodyRelationContact::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = body_part_particles_.size();
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]
				->searchNeighborsByParticles(total_real_particles,
											 *base_particles_, contact_configuration_[k],
											 get_body_part_particle_index_, *get_search_depths_[k],
											 *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	BodyPartRelationContact::BodyPartRelationContact(BodyPart &body_part, RealBodyVector contact_bodies)
		: BodyRelationContact(*body_part.getSPHBody(), contact_bodies), body_part_(&body_part),
		  body_part_particles_(DynamicCast<BodyPartByParticle>(this, body_part).body_part_particles_),
		  get_body_part_particle_index_(DynamicCast<BodyPartByParticle>(this, body_part).body_part_particles_)
	{
	}
	//=================================================================================================//
	void BodyPartRelationContact::updateConfiguration()
	{
		size_t number_of_particles = body_part_particles_.size();
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]
				->searchNeighborsByParticles(number_of_particles,
											 *base_particles_, contact_configuration_[k],
											 get_body_part_particle_index_, *get_search_depths_[k],
											 *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	BodyRelationContactToBodyPart::BodyRelationContactToBodyPart(RealBody &real_body, BodyPartVector contact_body_parts)
		: BodyRelationContact(real_body, contact_body_parts), contact_body_parts_(contact_body_parts)
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			get_part_contact_neighbors_.push_back(
				neighbor_relation_contact_body_part_ptr_vector_keeper_
					.createPtr<NeighborRelationContactBodyPart>(sph_body_, contact_body_parts[k]));
		}
	}
	//=================================================================================================//
	void BodyRelationContactToBodyPart::updateConfiguration()
	{
		size_t number_of_particles = base_particles_->total_real_particles_;
		for (size_t k = 0; k != contact_body_parts_.size(); ++k)
		{
			target_cell_linked_lists_[k]
				->searchNeighborsByParticles(number_of_particles,
											 *base_particles_, contact_configuration_[k],
											 get_particle_index_, *get_search_depths_[k],
											 *get_part_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
}
