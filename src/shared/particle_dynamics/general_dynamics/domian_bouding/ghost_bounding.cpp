#include "ghost_bounding.h"

namespace SPH
{
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::setupDynamics(Real dt)
{
    particles_->total_ghost_particles_ = 0;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::exec(Real dt)
{
    setupDynamics(dt);

    ghost_ranges_.[0].first =
        particle_for(execution::ParallelPolicy(), bound_cells_data_[0].first,
                     [&](size_t i)
                     { checkLowerBound(i, dt); });

    particle_for(execution::ParallelPolicy(), bound_cells_data_[1].first,
                 [&](size_t i)
                 { checkUpperBound(i, dt); });
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