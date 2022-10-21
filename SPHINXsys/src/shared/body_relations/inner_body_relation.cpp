#include "inner_body_relation.h"
#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	BodyRelationInner::BodyRelationInner(RealBody &real_body)
		: BaseBodyRelationInner(real_body), get_inner_neighbor_(real_body),
		  cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.cell_linked_list_)) {}
	//=================================================================================================//
	void BodyRelationInner::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		cell_linked_list_
			->searchNeighborsByParticles(base_particles_->total_real_particles_,
										 *base_particles_, inner_configuration_,
										 get_particle_index_, get_single_search_depth_,
										 get_inner_neighbor_);
	}
	//=================================================================================================//
	BodyRelationInnerVariableSmoothingLength::
		BodyRelationInnerVariableSmoothingLength(RealBody &real_body)
		: BaseBodyRelationInner(real_body), total_levels_(0),
		  get_inner_neighbor_variable_smoothing_length_(real_body)
	{
		MultilevelCellLinkedList *multi_level_cell_linked_list =
			DynamicCast<MultilevelCellLinkedList>(this, real_body.cell_linked_list_);
		cell_linked_list_levels_ = multi_level_cell_linked_list->getMeshLevels();
		total_levels_ = cell_linked_list_levels_.size();
		for (size_t l = 0; l != total_levels_; ++l)
		{
			get_multi_level_search_depth_.push_back(
				search_variable_smoothinglength_ptr_vector_keeper_
					.createPtr<SearchDepthVariableSmoothingLength>(real_body, cell_linked_list_levels_[l]));
		}
	}
	//=================================================================================================//
	void BodyRelationInnerVariableSmoothingLength::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		for (size_t l = 0; l != total_levels_; ++l)
		{
			cell_linked_list_levels_[l]
				->searchNeighborsByParticles(base_particles_->total_real_particles_,
											 *base_particles_, inner_configuration_, get_particle_index_,
											 *get_multi_level_search_depth_[l],
											 get_inner_neighbor_variable_smoothing_length_);
		}
	}
	//=================================================================================================//
	SolidBodyRelationSelfContact::
		SolidBodyRelationSelfContact(RealBody &real_body)
		: BaseBodyRelationInner(real_body),
		  body_surface_layer_(real_body),
		  body_part_particles_(body_surface_layer_.body_part_particles_),
		  get_body_part_particle_index_(body_part_particles_),
		  get_self_contact_neighbor_(real_body),
		  cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.cell_linked_list_)) {}
	//=================================================================================================//
	void SolidBodyRelationSelfContact::resetNeighborhoodCurrentSize()
	{
		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					size_t index_i = get_body_part_particle_index_(num);
					inner_configuration_[index_i].current_size_ = 0;
				}
			},
			ap);
	}
	//=================================================================================================//
	void SolidBodyRelationSelfContact::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = body_part_particles_.size();
		cell_linked_list_
			->searchNeighborsByParticles(total_real_particles,
										 *base_particles_, inner_configuration_,
										 get_body_part_particle_index_, get_single_search_depth_,
										 get_self_contact_neighbor_);
	}
	//=================================================================================================//
	void TreeBodyRelationInner::updateConfiguration()
	{
		generative_tree_.buildParticleConfiguration(inner_configuration_);
	}
	//=================================================================================================//
}
