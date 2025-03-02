#ifndef INITILIZATION_DYNAMICS_CK_HPP
#define INITILIZATION_DYNAMICS_CK_HPP

#include "initilization_dynamics_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
ContinuumInitialConditionCK::ContinuumInitialConditionCK(SPHBody &sph_body) 
: LocalDynamics(sph_body),
dv_pos_(this->particles_->template registerStateVariableOnly<Vecd>("Position")),
dv_vel_(this->particles_->template registerStateVariableOnly<Vecd>("Velocity")),
dv_stress_tensor_3D_(this->particles_->template registerStateVariableOnly<Mat3d>("StressTensor3D"))
{}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
ContinuumInitialConditionCK::UpdateKernel::
UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser):
pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy))
{};

} // namespace continuum_dynamics
} // namespace SPH
#endif //INITILIZATION_DYNAMICS_CK_HPP