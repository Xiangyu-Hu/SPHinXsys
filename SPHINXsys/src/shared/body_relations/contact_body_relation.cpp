/**
 * @file 	body_relation.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "contact_body_relation.h"

#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	ContactRelation::ContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: BaseContactRelation(sph_body, contact_sph_bodies)
	{
		initialization();
	}
	//=================================================================================================//
	ContactRelation::ContactRelation(SPHBody &sph_body, BodyPartVector contact_body_parts)
		: BaseContactRelation(sph_body, contact_body_parts)
	{
		initialization();
	}
	//=================================================================================================//
	void ContactRelation::initialization()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			CellLinkedList *target_cell_linked_list =
				DynamicCast<CellLinkedList>(this, contact_bodies_[k]->cell_linked_list_);
			target_cell_linked_lists_.push_back(target_cell_linked_list);
			get_search_depths_.push_back(
				search_depth_multi_resolution_ptr_vector_keeper_.createPtr<SearchDepthMultiResolution>(sph_body_, target_cell_linked_list));
			get_contact_neighbors_.push_back(
				neighbor_relation_contact_ptr_vector_keeper_.createPtr<NeighborBuilderContact>(sph_body_, *contact_bodies_[k]));
		}
	}
	//=================================================================================================//
	void ContactRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]->searchNeighborsByParticles(
				sph_body_, contact_configuration_[k],
				*get_search_depths_[k], *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================////=================================================================================================//
	AdaptiveContactRelation::AdaptiveContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: BaseContactRelation(sph_body, contact_sph_bodies)
	{
		cell_linked_list_levels_.resize(contact_bodies_.size());
		get_multi_level_search_range_.resize(contact_bodies_.size());

		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			MultilevelCellLinkedList *target_multi_level_mesh_cell_linked_lists_ =
				dynamic_cast<MultilevelCellLinkedList *>(contact_bodies_[k]->cell_linked_list_);
			cell_linked_list_levels_[k] = target_multi_level_mesh_cell_linked_lists_->getMeshLevels();

			for (size_t l = 0; l != cell_linked_list_levels_[k].size(); ++l)
			{
				get_multi_level_search_range_[k].push_back(
					adaptive_search_depth_ptr_vector_keeper_.createPtr<AdaptiveSearchDepth>(
						*contact_sph_bodies[k], cell_linked_list_levels_[k][l]));
			}

			get_contact_neighbors_.push_back(
				neighbor_relation_contact_ptr_vector_keeper_.createPtr<NeighborBuilderContact>(
					sph_body, *contact_sph_bodies[k]));
		}
	}
	//=================================================================================================//
	void AdaptiveContactRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			for (size_t l = 0; l != cell_linked_list_levels_[k].size(); ++l)
			{
				cell_linked_list_levels_[k][l]->searchNeighborsByParticles(
					sph_body_, contact_configuration_[k],
					*get_multi_level_search_range_[k][l], *get_contact_neighbors_[k]);
			}
		}
	}
	//=================================================================================================//
	SurfaceContactRelation::SurfaceContactRelation(SPHBody &sph_body, RealBodyVector contact_bodies)
		: BaseContactRelation(sph_body, contact_bodies),
		  body_surface_layer_(shape_surface_ptr_keeper_.createPtr<BodySurfaceLayer>(sph_body)),
		  body_part_particles_(body_surface_layer_->body_part_particles_)
	{
		initialization();
	}
	//=================================================================================================//
	SurfaceContactRelation::
		SurfaceContactRelation(SelfSurfaceContactRelation &solid_body_relation_self_contact,
							   RealBodyVector contact_bodies)
		: BaseContactRelation(*solid_body_relation_self_contact.real_body_, contact_bodies),
		  body_surface_layer_(&solid_body_relation_self_contact.body_surface_layer_),
		  body_part_particles_(body_surface_layer_->body_part_particles_)
	{
		initialization();
	}
	//=================================================================================================//
	void SurfaceContactRelation::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			parallel_for(
				blocked_range<size_t>(0, body_part_particles_.size()),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t num = r.begin(); num != r.end(); ++num)
					{
						size_t index_i = body_surface_layer_->getParticleIndex(num);
						contact_configuration_[k][index_i].current_size_ = 0;
					}
				},
				ap);
		}
	}
	//=================================================================================================//
	void SurfaceContactRelation::initialization()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			CellLinkedList *target_cell_linked_list =
				DynamicCast<CellLinkedList>(this, contact_bodies_[k]->cell_linked_list_);
			target_cell_linked_lists_.push_back(target_cell_linked_list);
			get_search_depths_.push_back(
				search_depth_multi_resolution_ptr_vector_keeper_.createPtr<SearchDepthMultiResolution>(sph_body_, target_cell_linked_list));
			get_contact_neighbors_.push_back(
				neighbor_relation_contact_ptr_vector_keeper_.createPtr<NeighborBuilderSolidContact>(sph_body_, *contact_bodies_[k]));
		}
	}
	//=================================================================================================//
	void SurfaceContactRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = body_part_particles_.size();
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]->searchNeighborsByParticles(
				*body_surface_layer_, contact_configuration_[k],
				*get_search_depths_[k], *get_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
	ContactRelationToBodyPart::ContactRelationToBodyPart(RealBody &real_body, BodyPartVector contact_body_parts)
		: ContactRelation(real_body, contact_body_parts), contact_body_parts_(contact_body_parts)
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			get_part_contact_neighbors_.push_back(
				neighbor_builder_contact_body_part_ptr_vector_keeper_
					.createPtr<NeighborBuilderContactBodyPart>(sph_body_, *contact_body_parts[k]));
		}
	}
	//=================================================================================================//
	void ContactRelationToBodyPart::updateConfiguration()
	{
		size_t number_of_particles = base_particles_.total_real_particles_;
		for (size_t k = 0; k != contact_body_parts_.size(); ++k)
		{
			target_cell_linked_lists_[k]->searchNeighborsByParticles(
				sph_body_, contact_configuration_[k],
				*get_search_depths_[k], *get_part_contact_neighbors_[k]);
		}
	}
	//=================================================================================================//
}
