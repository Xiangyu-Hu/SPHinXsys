/**
 * @file 	general_bounding.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "general_bounding.h"

namespace SPH
{
	//=================================================================================================//
	BoundingAlongAxis::
		BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis)
		: BaseDynamics<void>(), LocalDynamics(real_body),
		  DataDelegateSimple<SPHBody, BaseParticles>(real_body),
		  axis_(axis), bounding_bounds_(bounding_bounds), pos_(particles_->pos_),
		  cell_linked_list_(real_body.cell_linked_list_),
		  cut_off_radius_max_(real_body.sph_adaptation_->getKernel()->CutOffRadius()) {}
	//=================================================================================================//
	Vecd BasePeriodicCondition::
		setPeriodicTranslation(BoundingBox &bounding_bounds, int axis)
	{
		Vecd periodic_translation(0);
		periodic_translation[axis] =
			bounding_bounds.second[axis] - bounding_bounds.first[axis];
		return periodic_translation;
	}
	//=================================================================================================//
	BasePeriodicCondition::
		BasePeriodicCondition(RealBody &real_body, BoundingBox bounding_bounds, int axis)
		: periodic_translation_(setPeriodicTranslation(bounding_bounds, axis))
	{
		bound_cells_.resize(2);
		BaseCellLinkedList *cell_linked_list = real_body.cell_linked_list_;
		cell_linked_list->tagBoundingCells(bound_cells_, bounding_bounds, axis);
		if (periodic_translation_.norm() < real_body.sph_adaptation_->ReferenceSpacing())
		{
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			std::cout << "\n Periodic bounding failure: bounds not defined!" << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::checkLowerBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] < bounding_bounds_.first[axis_])
			pos_[index_i][axis_] += periodic_translation_[axis_];
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] > bounding_bounds_.second[axis_])
			pos_[index_i][axis_] -= periodic_translation_[axis_];
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::exec(Real dt)
	{
		setupDynamics(dt);

		particle_for(
			bound_cells_[0],
			[&](size_t i, Real delta)
			{ checkLowerBound(i, delta); },
			dt);

		particle_for(
			bound_cells_[1],
			[&](size_t i, Real delta)
			{ checkUpperBound(i, delta); },
			dt);
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		particle_parallel_for(
			bound_cells_[0],
			[&](size_t i, Real delta)
			{ checkLowerBound(i, delta); },
			dt);

		particle_parallel_for(
			bound_cells_[1],
			[&](size_t i, Real delta)
			{ checkUpperBound(i, delta); },
			dt);
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::exec(Real dt)
	{
		setupDynamics(dt);

		cell_list_for(
			bound_cells_[0],
			[&](CellList *cell_ist, Real delta)
			{ checkLowerBound(cell_ist, delta); },
			dt);

		cell_list_for(
			bound_cells_[1],
			[&](CellList *cell_ist, Real delta)
			{ checkUpperBound(cell_ist, delta); },
			dt);
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::
		PeriodicCellLinkedList::checkUpperBound(CellList *cell_list, Real dt)
	{
		ListDataVector &cell_list_data = cell_list->cell_list_data_;
		for (size_t num = 0; num < cell_list_data.size(); ++num)
		{
			Vecd particle_position = cell_list_data[num].second;
			if (particle_position[axis_] < bounding_bounds_.second[axis_] &&
				particle_position[axis_] > (bounding_bounds_.second[axis_] - cut_off_radius_max_))
			{
				Vecd translated_position = particle_position - periodic_translation_;
				/** insert ghost particle to cell linked list */
				cell_linked_list_->InsertACellLinkedListDataEntry(cell_list_data[num].first, translated_position);
			}
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::
		PeriodicCellLinkedList::checkLowerBound(CellList *cell_list, Real dt)
	{
		ListDataVector &cell_list_data = cell_list->cell_list_data_;
		for (size_t num = 0; num < cell_list_data.size(); ++num)
		{
			Vecd particle_position = cell_list_data[num].second;
			if (particle_position[axis_] > bounding_bounds_.first[axis_] &&
				particle_position[axis_] < (bounding_bounds_.first[axis_] + cut_off_radius_max_))
			{
				Vecd translated_position = particle_position + periodic_translation_;
				/** insert ghost particle to cell linked list */
				cell_linked_list_->InsertACellLinkedListDataEntry(cell_list_data[num].first, translated_position);
			}
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		CreatPeriodicGhostParticles::setupDynamics(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_.size(); ++i)
			ghost_particles_[i].clear();
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		CreatPeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] > bounding_bounds_.first[axis_] &&
			particle_position[axis_] < (bounding_bounds_.first[axis_] + cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[0].push_back(expected_particle_index);
			Vecd translated_position = particle_position + periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		CreatPeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] < bounding_bounds_.second[axis_] &&
			particle_position[axis_] > (bounding_bounds_.second[axis_] - cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[1].push_back(expected_particle_index);
			Vecd translated_position = particle_position - periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		UpdatePeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_[index_i] += periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		UpdatePeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_[index_i] -= periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		UpdatePeriodicGhostParticles::exec(Real dt)
	{
		setupDynamics(dt);

		particle_for(
			ghost_particles_[0],
			[&](size_t i, Real delta)
			{ checkLowerBound(i, delta); },
			dt);

		particle_for(
			ghost_particles_[1],
			[&](size_t i, Real delta)
			{ checkUpperBound(i, delta); },
			dt);
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::
		UpdatePeriodicGhostParticles::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		particle_parallel_for(
			ghost_particles_[0],
			[&](size_t i, Real delta)
			{ checkLowerBound(i, delta); },
			dt);

		particle_parallel_for(
			ghost_particles_[1],
			[&](size_t i, Real delta)
			{ checkUpperBound(i, delta); },
			dt);
	}
	//=================================================================================================//
	OpenBoundaryConditionAlongAxis::
		OpenBoundaryConditionAlongAxis(RealBody &real_body, BoundingBox bounding_bounds,
									   int axis, bool positive)
		: particle_type_transfer_(this->bound_cells_, real_body, bounding_bounds, axis, positive)
	{
		bound_cells_.resize(2);
		BaseCellLinkedList *cell_linked_list = real_body.cell_linked_list_;
		cell_linked_list->tagBoundingCells(bound_cells_, bounding_bounds, axis);
	}
	//=================================================================================================//
	void OpenBoundaryConditionAlongAxis ::
		ParticleTypeTransfer::checkLowerBound(size_t index_i, Real dt)
	{
		while (index_i < particles_->total_real_particles_ && pos_[index_i][axis_] < bounding_bounds_.first[axis_])
		{
			particles_->switchToBufferParticle(index_i);
		}
	}
	//=================================================================================================//
	void OpenBoundaryConditionAlongAxis ::
		ParticleTypeTransfer::checkUpperBound(size_t index_i, Real dt)
	{
		while (index_i < particles_->total_real_particles_ && pos_[index_i][axis_] > bounding_bounds_.second[axis_])
		{
			particles_->switchToBufferParticle(index_i);
		}
	}
	//=================================================================================================//
	void OpenBoundaryConditionAlongAxis::ParticleTypeTransfer::exec(Real dt)
	{
		setupDynamics(dt);

		// check lower bound
		CellLists &lower_bound_cells = bound_cells_[0];
		for (size_t i = 0; i != lower_bound_cells.size(); ++i)
		{
			IndexVector &particle_indexes = lower_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num)
				checking_bound_(particle_indexes[num], dt);
		}

		// check upper bound
		CellLists &upper_bound_cells = bound_cells_[1];
		for (size_t i = 0; i != upper_bound_cells.size(); ++i)
		{
			IndexVector &particle_indexes = upper_bound_cells[i]->real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num)
				checking_bound_(particle_indexes[num], dt);
		}
	}
	//=================================================================================================//
	MirrorConditionAlongAxis::MirrorBounding::
		MirrorBounding(CellLists &bound_cells, RealBody &real_body,
					   BoundingBox bounding_bounds, int axis, bool positive)
		: BoundingAlongAxis(real_body, bounding_bounds, axis),
		  bound_cells_(bound_cells), vel_(particles_->vel_)
	{
		checking_bound_ =
			positive ? std::bind(&MirrorConditionAlongAxis::MirrorBounding::checkUpperBound, this, _1, _2)
					 : std::bind(&MirrorConditionAlongAxis::MirrorBounding::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	MirrorConditionAlongAxis::CreatingMirrorGhostParticles::
		CreatingMirrorGhostParticles(IndexVector &ghost_particles, CellLists &bound_cells, RealBody &real_body,
									 BoundingBox bounding_bounds, int axis, bool positive)
		: MirrorBounding(bound_cells, real_body, bounding_bounds, axis, positive),
		  ghost_particles_(ghost_particles) {}
	//=================================================================================================//
	MirrorConditionAlongAxis::UpdatingMirrorGhostStates::
		UpdatingMirrorGhostStates(IndexVector &ghost_particles, CellLists &bound_cells, RealBody &real_body,
								  BoundingBox bounding_bounds, int axis, bool positive)
		: MirrorBounding(bound_cells, real_body, bounding_bounds, axis, positive), ghost_particles_(ghost_particles)
	{
		checking_bound_update_ =
			positive ? std::bind(&MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkUpperBound, this, _1, _2)
					 : std::bind(&MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis ::MirrorBounding::checkLowerBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] < bounding_bounds_.first[axis_])
		{
			mirrorAlongAxis(index_i, bounding_bounds_.first, axis_);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding ::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] > bounding_bounds_.second[axis_])
		{
			mirrorAlongAxis(index_i, bounding_bounds_.second, axis_);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding::
		mirrorAlongAxis(size_t particle_index_i, Vecd body_bound, int axis)
	{
		pos_[particle_index_i][axis] = 2.0 * body_bound[axis] - pos_[particle_index_i][axis];
		vel_[particle_index_i][axis] *= -1.0;
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding::exec(Real dt)
	{
		setupDynamics(dt);
		for (size_t i = 0; i != bound_cells_.size(); ++i)
		{
			ListDataVector &list_data = bound_cells_[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num)
				checking_bound_(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding ::parallel_exec(Real dt)
	{
		setupDynamics(dt);
		parallel_for(
			blocked_range<size_t>(0, bound_cells_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ListDataVector &list_data = bound_cells_[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num)
						checking_bound_(list_data[num].first, dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::CreatingMirrorGhostParticles ::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] > bounding_bounds_.first[axis_] &&
			particle_position[axis_] < (bounding_bounds_.first[axis_] + cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorAlongAxis(expected_particle_index, bounding_bounds_.first, axis_);
			Vecd translated_position = particles_->pos_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::CreatingMirrorGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] < bounding_bounds_.second[axis_] &&
			particle_position[axis_] > (bounding_bounds_.second[axis_] - cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorAlongAxis(expected_particle_index, bounding_bounds_.second, axis_);
			Vecd translated_position = particles_->pos_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_->InsertACellLinkedListDataEntry(expected_particle_index, translated_position);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorAlongAxis(index_i, bounding_bounds_.first, axis_);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates ::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorAlongAxis(index_i, bounding_bounds_.second, axis_);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates ::exec(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_.size(); ++i)
		{
			checking_bound_update_(ghost_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates ::parallel_exec(Real dt)
	{
		parallel_for(
			blocked_range<size_t>(0, ghost_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					checking_bound_update_(ghost_particles_[i], dt);
				}
			},
			ap);
	}
	//=================================================================================================//
}
//=================================================================================================//