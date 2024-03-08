#include "bidirectional_buffer.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
template <>
void BidirectionalBuffer::Injection<execution::ParallelPolicy>::update(size_t index_i, Real dt)
{
    size_t sorted_index_i = buffer_particle_list_[index_i];
    if (aligned_box_.checkUpperBound(axis_, pos_n_[sorted_index_i]) && buffer_particle_indicator_[sorted_index_i] == 1)
    {
        mutex_switch_to_real_.lock();
        if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
        {
            std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
                      << "Not enough body buffer particles! Exit the code."
                      << "\n";
            exit(0);
        }
        buffer_particle_indicator_[sorted_index_i] = 0;
        particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
        // {
        //     // update sorted id
        //     size_t last_unsorted_index = particles_->unsorted_id_[particles_->total_real_particles_];
        //     particles_->sorted_id_[last_unsorted_index] = particles_->total_real_particles_;
        // }
        particles_->total_real_particles_ += 1;
        mutex_switch_to_real_.unlock();
        pos_n_[sorted_index_i] = aligned_box_.getUpperPeriodic(axis_, pos_n_[sorted_index_i]);
        Real sound_speed = fluid_.getSoundSpeed(rho_n_[sorted_index_i]);
        if (prescribe_pressure_)
        {
            p_[sorted_index_i] = buffer_.getTargetPressure(dt);
            rho_n_[sorted_index_i] = p_[sorted_index_i] / pow(sound_speed, 2.0) + fluid_.ReferenceDensity();
        }
        previous_surface_indicator_[sorted_index_i] = 1;
    }
}
} // namespace fluid_dynamics
  //=====================================================================================================//
} // namespace SPH
  //=========================================================================================================//