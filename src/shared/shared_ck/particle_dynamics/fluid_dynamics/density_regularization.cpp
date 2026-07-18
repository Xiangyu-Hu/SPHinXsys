#include "density_regularization.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AverageCompression::AverageCompression(SPHBody &sph_body)
    : BaseDynamicsType(sph_body),
      dv_compression_sum_(particles_->getVariableByName<Real>("CompressionSummation")),
      sv_compression_inv_ref_(particles_->getSingleVariableByName<Real>("InverseReferenceCompression")) {}
//=================================================================================================//
AverageCompression::FinishDynamics::FinishDynamics(AverageCompression &encloser)
    : sv_compression_inv_ref_(encloser.sv_compression_inv_ref_) {}
//=================================================================================================//
Real AverageCompression::FinishDynamics::Result(const ReduceReturnType &reduced_value)
{
    Real average_compression = reduced_value.first / reduced_value.second;
    sv_compression_inv_ref_->setValue(1.0 / average_compression);
    return average_compression;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
