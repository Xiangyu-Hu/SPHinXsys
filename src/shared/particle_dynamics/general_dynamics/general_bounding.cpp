#include "general_bounding.h"

namespace SPH
{
//=================================================================================================//
BoundingAlongAxis::BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis)
    : BaseDynamics<void>(real_body), LocalDynamics(real_body),
      GeneralDataDelegateSimple(real_body),
      axis_(axis), bounding_bounds_(bounding_bounds),
      pos_(particles_->pos_),
      cell_linked_list_(real_body.getCellLinkedList()),
      cut_off_radius_max_(real_body.sph_adaptation_->getKernel()->CutOffRadius()) {}
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

    particle_for(execution::ParallelPolicy(), bound_cells_data_[0].second,
                 [&](ListDataVector *cell_ist)
                 { checkLowerBound(*cell_ist, dt); });

    particle_for(execution::ParallelPolicy(), bound_cells_data_[1].second,
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

    particle_for(execution::ParallelPolicy(), ghost_particles_[0],
                 [&](size_t i)
                 { checkLowerBound(i, dt); });

    particle_for(execution::ParallelPolicy(), ghost_particles_[1],
                 [&](size_t i)
                 { checkUpperBound(i, dt); });
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//