#include "ghost_bounding.h"

namespace SPH
{
//=================================================================================================//
Ghost<PeriodicAlongAxis>::Ghost(BoundingBoxd bounding_bounds, int axis)
    : Ghost<Base>(), PeriodicAlongAxis(bounding_bounds, axis) {}
//=================================================================================================//
void Ghost<PeriodicAlongAxis>::reserveGhostParticles(BaseParticles &base_particles, Real particle_spacing)
{
    ghost_size_ = calculateGhostSize(particle_spacing);

    lower_ghost_bound_.first = base_particles.allocateGhostParticles(ghost_size_);
    upper_ghost_bound_.first = base_particles.allocateGhostParticles(ghost_size_);

    is_particles_reserved_ = true;
}
//=================================================================================================//
size_t Ghost<PeriodicAlongAxis>::calculateGhostSize(Real particle_spacing)
{
    int next_axis = NextAxis(axis_);
    Real bound_size = bounding_bounds_.upper_[next_axis] - bounding_bounds_.lower_[next_axis];
    Real ghost_width = 4.0;
    return std::ceil(2.0 * ghost_width * ABS(bound_size) / particle_spacing);
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::
    PeriodicConditionUsingGhostParticles(RealBody &real_body, Ghost<PeriodicAlongAxis> &ghost_boundary)
    : BasePeriodicCondition<execution::ParallelPolicy>(real_body, ghost_boundary),
      bounding_(bound_cells_data_, real_body, ghost_boundary),
      ghost_creation_(bound_cells_data_, real_body, ghost_boundary),
      ghost_update_(bound_cells_data_, real_body, ghost_boundary)
{
    ghost_boundary.checkParticlesReserved();
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::
    CreatPeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                Ghost<PeriodicAlongAxis> &ghost_boundary)
    : PeriodicBounding(bound_cells_data, real_body, ghost_boundary),
      ghost_boundary_(ghost_boundary),
      lower_ghost_bound_(ghost_boundary.LowerGhostBound()),
      upper_ghost_bound_(ghost_boundary.UpperGhostBound()),
      cell_linked_list_(real_body.getCellLinkedList()),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::setupDynamics(Real dt)
{
    PeriodicBounding::setupDynamics(dt);
    lower_ghost_bound_.second = lower_ghost_bound_.first;
    upper_ghost_bound_.second = upper_ghost_bound_.first;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::exec(Real dt)
{
    setupDynamics(dt);

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
    if (particle_position[axis_] > bounding_bounds_.lower_[axis_] &&
        particle_position[axis_] < (bounding_bounds_.lower_[axis_] + cut_off_radius_max_))
    {
        mutex_create_ghost_particle_.lock();
        particles_->updateGhostParticle(lower_ghost_bound_.second, index_i);
        pos_[lower_ghost_bound_.second] = particle_position + periodic_translation_;
        /** insert ghost particle to cell linked list */
        cell_linked_list_.InsertListDataEntry(lower_ghost_bound_.second, pos_[lower_ghost_bound_.second]);
        lower_ghost_bound_.second++;
        ghost_boundary_.checkWithinGhostSize(lower_ghost_bound_);
        mutex_create_ghost_particle_.unlock();
    }
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
{
    Vecd particle_position = pos_[index_i];
    if (particle_position[axis_] < bounding_bounds_.upper_[axis_] &&
        particle_position[axis_] > (bounding_bounds_.upper_[axis_] - cut_off_radius_max_))
    {
        mutex_create_ghost_particle_.lock();
        particles_->updateGhostParticle(upper_ghost_bound_.second, index_i);
        pos_[upper_ghost_bound_.second] = particle_position - periodic_translation_;
        /** insert ghost particle to cell linked list */
        cell_linked_list_.InsertListDataEntry(upper_ghost_bound_.second, pos_[upper_ghost_bound_.second]);
        upper_ghost_bound_.second++;
        ghost_boundary_.checkWithinGhostSize(upper_ghost_bound_);
        mutex_create_ghost_particle_.unlock();
    }
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::
    UpdatePeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                 Ghost<PeriodicAlongAxis> &ghost_boundary)
    : PeriodicBounding(bound_cells_data, real_body, ghost_boundary),
      ghost_boundary_(ghost_boundary),
      lower_ghost_bound_(ghost_boundary.LowerGhostBound()),
      upper_ghost_bound_(ghost_boundary.UpperGhostBound()),
      sorted_id_(particles_->ParticleSortedIds()) {}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
{
    particles_->copyFromAnotherParticle(index_i, sorted_id_[index_i]);
    pos_[index_i] += periodic_translation_;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
{
    particles_->copyFromAnotherParticle(index_i, sorted_id_[index_i]);
    pos_[index_i] -= periodic_translation_;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::exec(Real dt)
{
    setupDynamics(dt);

    particle_for(execution::ParallelPolicy(), ghost_boundary_.getGhostParticleRange(lower_ghost_bound_),
                 [&](size_t i)
                 { checkLowerBound(i, dt); });

    particle_for(execution::ParallelPolicy(), ghost_boundary_.getGhostParticleRange(upper_ghost_bound_),
                 [&](size_t i)
                 { checkUpperBound(i, dt); });
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//