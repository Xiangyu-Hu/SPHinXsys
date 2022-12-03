#include "general_bounding.h"

namespace SPH
{
	//=================================================================================================//
	BoundingAlongAxis::BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis)
		: BaseDynamics<void>(), LocalDynamics(real_body),
		  GeneralDataDelegateSimple(real_body),
		  axis_(axis), bounding_bounds_(bounding_bounds),
		  pos_(particles_->pos_),
		  cell_linked_list_(real_body.getCellLinkedList()),
		  cut_off_radius_max_(real_body.sph_adaptation_->getKernel()->CutOffRadius()) {}
	//=================================================================================================//
	Vecd BasePeriodicCondition::setPeriodicTranslation(BoundingBox &bounding_bounds, int axis)
	{
		Vecd periodic_translation = Vecd::Zero();
		periodic_translation[axis] =
			bounding_bounds.second_[axis] - bounding_bounds.first_[axis];
		return periodic_translation;
	}
	//=================================================================================================//
	BasePeriodicCondition::BasePeriodicCondition(RealBody &real_body, BoundingBox bounding_bounds, int axis)
		: periodic_translation_(setPeriodicTranslation(bounding_bounds, axis))
	{
		bound_cells_data_.resize(2);
		BaseCellLinkedList &cell_linked_list = real_body.getCellLinkedList();
		cell_linked_list.tagBoundingCells(bound_cells_data_, bounding_bounds, axis);
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
		if (pos_[index_i][axis_] < bounding_bounds_.first_[axis_])
			pos_[index_i][axis_] += periodic_translation_[axis_];
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] > bounding_bounds_.second_[axis_])
			pos_[index_i][axis_] -= periodic_translation_[axis_];
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::exec(Real dt)
	{
		setupDynamics(dt);

		particle_for(bound_cells_data_[0].first,
					 [&](size_t i)
					 { checkLowerBound(i, dt); });

		particle_for(bound_cells_data_[1].first,
					 [&](size_t i)
					 { checkUpperBound(i, dt); });
	}
	//=================================================================================================//
	void BasePeriodicCondition::PeriodicBounding::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		particle_parallel_for(bound_cells_data_[0].first,
							  [&](size_t i)
							  { checkLowerBound(i, dt); });

		particle_parallel_for(bound_cells_data_[1].first,
							  [&](size_t i)
							  { checkUpperBound(i, dt); });
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::
		PeriodicCellLinkedList::checkUpperBound(ListDataVector &cell_list_data, Real dt)
	{
		for (size_t num = 0; num < cell_list_data.size(); ++num)
		{
			Vecd particle_position = std::get<1>(cell_list_data[num]);
			if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
				particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
			{
				Vecd translated_position = particle_position - periodic_translation_;
				/** insert ghost particle to cell linked list */
				mutex_cell_list_entry_.lock();
				cell_linked_list_.InsertListDataEntry(std::get<0>(cell_list_data[num]),
													   translated_position, std::get<2>(cell_list_data[num]));
				mutex_cell_list_entry_.unlock();
			}
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::
		PeriodicCellLinkedList::checkLowerBound(ListDataVector &cell_list_data, Real dt)
	{
		for (size_t num = 0; num < cell_list_data.size(); ++num)
		{
			Vecd particle_position = std::get<1>(cell_list_data[num]);
			if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
				particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
			{
				Vecd translated_position = particle_position + periodic_translation_;
				/** insert ghost particle to cell linked list */
				mutex_cell_list_entry_.lock();
				cell_linked_list_.InsertListDataEntry(std::get<0>(cell_list_data[num]),
													   translated_position, std::get<2>(cell_list_data[num]));
				mutex_cell_list_entry_.unlock();
			}
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::exec(Real dt)
	{
		setupDynamics(dt);

		cell_list_for(bound_cells_data_[0].second,
					  [&](ListDataVector *cell_ist)
					  { checkLowerBound(*cell_ist, dt); });

		cell_list_for(bound_cells_data_[1].second,
					  [&](ListDataVector *cell_ist)
					  { checkUpperBound(*cell_ist, dt); });
	}
	//=================================================================================================//
	void PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		cell_list_parallel_for(bound_cells_data_[0].second,
							   [&](ListDataVector *cell_ist)
							   { checkLowerBound(*cell_ist, dt); });

		cell_list_parallel_for(bound_cells_data_[1].second,
							   [&](ListDataVector *cell_ist)
							   { checkUpperBound(*cell_ist, dt); });
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::setupDynamics(Real dt)
	{
		for (size_t i = 0; i != ghost_particles_.size(); ++i)
			ghost_particles_[i].clear();
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		Real particle_volumetric = Vol_[index_i];
		if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
			particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
		{
			mutex_create_ghost_particle_.lock();
			size_t ghost_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[0].push_back(ghost_particle_index);
			pos_[ghost_particle_index] = particle_position + periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_.InsertListDataEntry(ghost_particle_index,
												   pos_[ghost_particle_index], particle_volumetric);
			mutex_create_ghost_particle_.unlock();
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		Real particle_volumetric = Vol_[index_i];
		if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
			particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
		{
			mutex_create_ghost_particle_.lock();
			size_t ghost_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_[1].push_back(ghost_particle_index);
			pos_[ghost_particle_index] = particle_position - periodic_translation_;
			/** insert ghost particle to cell linked list */
			cell_linked_list_.InsertListDataEntry(ghost_particle_index,
												   pos_[ghost_particle_index], particle_volumetric);
			mutex_create_ghost_particle_.unlock();
		}
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_[index_i] += periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		pos_[index_i] -= periodic_translation_;
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::exec(Real dt)
	{
		setupDynamics(dt);

		particle_for(ghost_particles_[0],
					 [&](size_t i)
					 { checkLowerBound(i, dt); });

		particle_for(ghost_particles_[1],
					 [&](size_t i)
					 { checkUpperBound(i, dt); });
	}
	//=================================================================================================//
	void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::parallel_exec(Real dt)
	{
		setupDynamics(dt);

		particle_parallel_for(ghost_particles_[0],
							  [&](size_t i)
							  { checkLowerBound(i, dt); });

		particle_parallel_for(ghost_particles_[1],
							  [&](size_t i)
							  { checkUpperBound(i, dt); });
	}
	//=================================================================================================//
	MirrorConditionAlongAxis::MirrorBounding::
		MirrorBounding(CellLists bound_cells_data, RealBody &real_body,
					   BoundingBox bounding_bounds, int axis, bool positive)
		: BoundingAlongAxis(real_body, bounding_bounds, axis),
		  bound_cells_data_(bound_cells_data), vel_(particles_->vel_)
	{
		checking_bound_ =
			positive ? std::bind(&MirrorConditionAlongAxis::MirrorBounding::checkUpperBound, this, _1, _2)
					 : std::bind(&MirrorConditionAlongAxis::MirrorBounding::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	MirrorConditionAlongAxis::CreatingMirrorGhostParticles::
		CreatingMirrorGhostParticles(IndexVector &ghost_particles,
									 CellLists bound_cells_data,
									 RealBody &real_body,
									 BoundingBox bounding_bounds, int axis, bool positive)
		: MirrorBounding(bound_cells_data, real_body, bounding_bounds, axis, positive),
		  ghost_particles_(ghost_particles) {}
	//=================================================================================================//
	MirrorConditionAlongAxis::UpdatingMirrorGhostStates::
		UpdatingMirrorGhostStates(IndexVector &ghost_particles,
								  CellLists bound_cells_data, RealBody &real_body,
								  BoundingBox bounding_bounds, int axis, bool positive)
		: MirrorBounding(bound_cells_data, real_body, bounding_bounds, axis, positive),
		  ghost_particles_(ghost_particles)
	{
		checking_bound_update_ =
			positive ? std::bind(&MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkUpperBound, this, _1, _2)
					 : std::bind(&MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkLowerBound, this, _1, _2);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis ::MirrorBounding::checkLowerBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] < bounding_bounds_.first_[axis_])
		{
			mirrorAlongAxis(index_i, bounding_bounds_.first_, axis_);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding ::checkUpperBound(size_t index_i, Real dt)
	{
		if (pos_[index_i][axis_] > bounding_bounds_.second_[axis_])
		{
			mirrorAlongAxis(index_i, bounding_bounds_.second_, axis_);
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
		DataListsInCells &cell_list_data = bound_cells_data_.second;
		for (size_t i = 0; i != cell_list_data.size(); ++i)
		{
			ListDataVector &list_data = *cell_list_data[i];
			for (size_t num = 0; num < list_data.size(); ++num)
				checking_bound_(std::get<0>(list_data[num]), dt);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::MirrorBounding ::parallel_exec(Real dt)
	{
		setupDynamics(dt);
		DataListsInCells &cell_list_data = bound_cells_data_.second;
		parallel_for(
			blocked_range<size_t>(0, cell_list_data.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ListDataVector &list_data = *cell_list_data[i];
					for (size_t num = 0; num < list_data.size(); ++num)
						checking_bound_(std::get<0>(list_data[num]), dt);
				}
			},
			ap);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::CreatingMirrorGhostParticles ::checkLowerBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
			particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorAlongAxis(expected_particle_index, bounding_bounds_.first_, axis_);
			Vecd translated_position = particles_->pos_[expected_particle_index];
			Real particle_volumetric = particles_->Vol_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_.InsertListDataEntry(
				expected_particle_index, translated_position, particle_volumetric);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::CreatingMirrorGhostParticles::checkUpperBound(size_t index_i, Real dt)
	{
		Vecd particle_position = pos_[index_i];
		if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
			particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
		{
			size_t expected_particle_index = particles_->insertAGhostParticle(index_i);
			ghost_particles_.push_back(expected_particle_index);
			/** mirror boundary condition */
			mirrorAlongAxis(expected_particle_index, bounding_bounds_.second_, axis_);
			Vecd translated_position = particles_->pos_[expected_particle_index];
			Real particle_volumetric = particles_->Vol_[expected_particle_index];
			/** insert ghost particle to cell linked list */
			cell_linked_list_.InsertListDataEntry(
				expected_particle_index, translated_position, particle_volumetric);
		}
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates::checkLowerBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorAlongAxis(index_i, bounding_bounds_.first_, axis_);
	}
	//=================================================================================================//
	void MirrorConditionAlongAxis::UpdatingMirrorGhostStates ::checkUpperBound(size_t index_i, Real dt)
	{
		particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
		mirrorAlongAxis(index_i, bounding_bounds_.second_, axis_);
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