#include "pressure_boundary.h"

namespace SPH
{
//=================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
FlowPressureBuffer::FlowPressureBuffer(BodyPartByCell &body_part, Vecd normal_vector)
    : BaseFlowBoundaryCondition(body_part),
      kernel_sum_(*particles_->getVariableByName<Vecd>("KernelSummation")), direction_(normal_vector){};
//=================================================================================================//
void FlowPressureBuffer::update(size_t index_i, Real dt)
{
    vel_[index_i] += 2.0 * kernel_sum_[index_i] * getTargetPressure(index_i,dt) / rho_[index_i] * dt;

    vel_[index_i] = vel_[index_i].dot(direction_) * direction_;
}
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//