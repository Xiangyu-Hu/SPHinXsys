/**
 * @file 	body_relation.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "body_relation.h"

#include "base_kernel.h"
#include "base_particles.h"
#include "cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	SPHBodyRelation::SPHBodyRelation(SPHBody &sph_body)
		: sph_body_(&sph_body), base_particles_(sph_body.base_particles_) {}
	//=================================================================================================//
	BaseBodyRelationInner::BaseBodyRelationInner(RealBody &real_body)
		: SPHBodyRelation(real_body), real_body_(&real_body)
	{
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseBodyRelationInner::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	void BaseBodyRelationInner::resetNeighborhoodCurrentSize()
	{
		parallel_for(
			blocked_range<size_t>(0, base_particles_->total_real_particles_),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					inner_configuration_[num].current_size_ = 0;
				}
			},
			ap);
	}
	//=================================================================================================//
	BodyRelationInner::BodyRelationInner(RealBody &real_body)
		: BaseBodyRelationInner(real_body), get_inner_neighbor_(&real_body),
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
		  get_inner_neighbor_variable_smoothing_length_(&real_body)
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
		  get_self_contact_neighbor_(&real_body),
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
	BaseBodyRelationContact::BaseBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: SPHBodyRelation(sph_body), contact_bodies_(contact_sph_bodies)
	{
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	BaseBodyRelationContact::BaseBodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts)
		: SPHBodyRelation(sph_body)
	{
		for (size_t k = 0; k != contact_body_parts.size(); ++k)
		{
			contact_bodies_.push_back(DynamicCast<RealBody>(this, contact_body_parts[k]->getSPHBody()));
		}
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseBodyRelationContact::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_->base_particles_->real_particles_bound_;
		contact_configuration_.resize(contact_bodies_.size());
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	void BaseBodyRelationContact::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			parallel_for(
				blocked_range<size_t>(0, base_particles_->total_real_particles_),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t num = r.begin(); num != r.end(); ++num)
					{
						contact_configuration_[k][num].current_size_ = 0;
					}
				},
				ap);
		}
	}
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
	void GenerativeBodyRelationInner::updateConfiguration()
	{
		generative_structure_->buildParticleConfiguration(*base_particles_, inner_configuration_);
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
	ComplexBodyRelation::
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation)
		: SPHBodyRelation(*inner_relation.sph_body_),
		  inner_relation_(inner_relation),
		  contact_relation_(contact_relation),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(RealBody &real_body, RealBodyVector contact_bodies)
		: SPHBodyRelation(real_body),
		  inner_relation_(base_body_relation_inner_ptr_keeper_.createRef<BodyRelationInner>(real_body)),
		  contact_relation_(base_body_relation_contact_ptr_keeper_
								.createRef<BodyRelationContact>(real_body, contact_bodies)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, RealBodyVector contact_bodies)
		: SPHBodyRelation(*inner_relation.sph_body_),
		  inner_relation_(inner_relation),
		  contact_relation_(base_body_relation_contact_ptr_keeper_.createRef<BodyRelationContact>(
			  DynamicCast<RealBody>(this, *sph_body_), contact_bodies)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	ComplexBodyRelation::ComplexBodyRelation(RealBody &real_body, BodyPartVector contact_body_parts)
		: SPHBodyRelation(real_body),
		  inner_relation_(base_body_relation_inner_ptr_keeper_.createRef<BodyRelationInner>(real_body)),
		  contact_relation_(base_body_relation_contact_ptr_keeper_
								.createRef<BodyRelationContactToBodyPart>(real_body, contact_body_parts)),
		  contact_bodies_(contact_relation_.contact_bodies_),
		  inner_configuration_(inner_relation_.inner_configuration_),
		  contact_configuration_(contact_relation_.contact_configuration_)
	{
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfigurationMemories()
	{
		inner_relation_.updateConfigurationMemories();
		contact_relation_.updateConfigurationMemories();
	}
	//=================================================================================================//
	void ComplexBodyRelation::updateConfiguration()
	{
		inner_relation_.updateConfiguration();
		contact_relation_.updateConfiguration();
	}
	//=================================================================================================//
}
