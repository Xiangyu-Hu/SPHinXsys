#include "relaxation_stepping_ck.h"

namespace SPH
{
//=================================================================================================//
RelaxationScalingCK::RelaxationScalingCK(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      dv_residual_(particles_->getVariableByName<Vecd>("ZeroGradientResidual")),
      h_ref_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
RelaxationScalingCK::FinishDynamics::FinishDynamics(RelaxationScalingCK &encloser)
    : h_ref_(encloser.h_ref_) {}
//=================================================================================================//
Real RelaxationScalingCK::FinishDynamics::Result(Real reduced_value)
{
    return 0.0625 * h_ref_ / (reduced_value + TinyReal);
}
//=================================================================================================//
PositionRelaxationCK::PositionRelaxationCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableByName<Vecd>("Position")),
      residual_(particles_->getVariableByName<Vecd>("ZeroGradientResidual")) {}
//=================================================================================================//
} // namespace SPH
