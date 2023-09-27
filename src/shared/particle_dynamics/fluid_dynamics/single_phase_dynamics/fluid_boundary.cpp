#include "fluid_boundary.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BaseFlowBoundaryCondition::BaseFlowBoundaryCondition(BodyPartByCell &body_part)
    : BaseLocalDynamics<BodyPartByCell>(body_part), FluidDataSimple(sph_body_),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
      pos_(particles_->pos_), vel_(particles_->vel_){};
//=================================================================================================//
FlowVelocityBuffer::FlowVelocityBuffer(BodyPartByCell &body_part, Real relaxation_rate)
    : BaseFlowBoundaryCondition(body_part), relaxation_rate_(relaxation_rate){};
//=================================================================================================//
void FlowVelocityBuffer::update(size_t index_i, Real dt)
{
    vel_[index_i] += relaxation_rate_ * (getTargetVelocity(pos_[index_i], vel_[index_i]) - vel_[index_i]);
}
//=================================================================================================//
DampingBoundaryCondition::DampingBoundaryCondition(BodyRegionByCell &body_part)
    : BaseFlowBoundaryCondition(body_part), strength_(5.0),
      damping_zone_bounds_(body_part.body_part_shape_.getBounds()){};
//=================================================================================================//
void DampingBoundaryCondition::update(size_t index_i, Real dt)
{
    Real damping_factor = (pos_[index_i][0] - damping_zone_bounds_.first_[0]) /
                          (damping_zone_bounds_.second_[0] - damping_zone_bounds_.first_[0]);
    vel_[index_i] *= (1.0 - dt * strength_ * damping_factor * damping_factor);
}
//=================================================================================================//
EmitterInflowCondition::
    EmitterInflowCondition(BodyAlignedBoxByParticle &aligned_box_part)
    : BaseLocalDynamics<BodyPartByParticle>(aligned_box_part), FluidDataSimple(sph_body_),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      pos_(particles_->pos_), vel_(particles_->vel_), acc_(particles_->acc_),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
      drho_dt_(*particles_->getVariableByName<Real>("DensityChangeRate")),
      inflow_pressure_(0), rho0_(fluid_.ReferenceDensity()),
      aligned_box_(aligned_box_part.aligned_box_),
      updated_transform_(aligned_box_.getTransform()),
      old_transform_(updated_transform_) {}
//=================================================================================================//
void EmitterInflowCondition ::update(size_t unsorted_index_i, Real dt)
{
    size_t sorted_index_i = sorted_id_[unsorted_index_i];
    Vecd frame_position = old_transform_.shiftBaseStationToFrame(pos_[sorted_index_i]);
    Vecd frame_velocity = old_transform_.xformBaseVecToFrame(vel_[sorted_index_i]);
    pos_[sorted_index_i] = updated_transform_.shiftFrameStationToBase(frame_position);
    vel_[sorted_index_i] = updated_transform_.xformFrameVecToBase(getTargetVelocity(frame_position, frame_velocity));
    rho_[sorted_index_i] = rho0_;
    p_[sorted_index_i] = fluid_.getPressure(rho0_);
}
//=================================================================================================//
EmitterInflowInjection::EmitterInflowInjection(BodyAlignedBoxByParticle &aligned_box_part,
                                               size_t body_buffer_width, int axis)
    : BaseLocalDynamics<BodyPartByParticle>(aligned_box_part), FluidDataSimple(sph_body_),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      pos_(particles_->pos_), rho_(particles_->rho_),
      p_(*particles_->getVariableByName<Real>("Pressure")),
      axis_(axis), aligned_box_(aligned_box_part.aligned_box_)
{
    size_t total_body_buffer_particles = aligned_box_part.body_part_particles_.size() * body_buffer_width;
    particles_->addBufferParticles(total_body_buffer_particles);
    sph_body_.allocateConfigurationMemoriesForBufferParticles();
}
//=================================================================================================//
void EmitterInflowInjection::update(size_t unsorted_index_i, Real dt)
{
    size_t sorted_index_i = sorted_id_[unsorted_index_i];
    if (aligned_box_.checkUpperBound(axis_, pos_[sorted_index_i]))
    {
        mutex_switch_to_real_.lock();
        if (particles_->total_real_particles_ >= particles_->real_particles_bound_)
        {
            std::cout << "EmitterInflowBoundaryCondition::ConstraintAParticle: \n"
                      << "Not enough body buffer particles! Exit the code."
                      << "\n";
            exit(0);
        }
        /** Buffer Particle state copied from real particle. */
        particles_->copyFromAnotherParticle(particles_->total_real_particles_, sorted_index_i);
        /** Realize the buffer particle by increasing the number of real particle in the body.  */
        particles_->total_real_particles_ += 1;
        mutex_switch_to_real_.unlock();
        /** Periodic bounding. */
        pos_[sorted_index_i] = aligned_box_.getUpperPeriodic(axis_, pos_[sorted_index_i]);
        rho_[sorted_index_i] = fluid_.ReferenceDensity();
        p_[sorted_index_i] = fluid_.getPressure(rho_[sorted_index_i]);
    }
}
//=================================================================================================//
DisposerOutflowDeletion::
    DisposerOutflowDeletion(BodyAlignedBoxByCell &aligned_box_part, int axis)
    : BaseLocalDynamics<BodyPartByCell>(aligned_box_part), FluidDataSimple(sph_body_),
      pos_(particles_->pos_), axis_(axis), aligned_box_(aligned_box_part.aligned_box_) {}
//=================================================================================================//
void DisposerOutflowDeletion::update(size_t index_i, Real dt)
{
    mutex_switch_to_buffer_.lock();
    while (aligned_box_.checkUpperBound(axis_, pos_[index_i]) && index_i < particles_->total_real_particles_)
    {
        particles_->switchToBufferParticle(index_i);
    }
    mutex_switch_to_buffer_.unlock();
}
} // namespace fluid_dynamics
} // namespace SPH
