#include "common_compressible_FVM_classes.h"

namespace SPH
{
//=================================================================================================//
CompressibleAcousticTimeStepSizeInFVM::
    CompressibleAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real min_distance_between_nodes, Real acousticCFL)
    : AcousticTimeStep(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      min_distance_between_nodes_(min_distance_between_nodes),
      compressible_fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      acousticCFL_(acousticCFL){};
//=================================================================================================//
Real CompressibleAcousticTimeStepSizeInFVM::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real CompressibleAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
{
    return acousticCFL_ / Real(Dimensions) * min_distance_between_nodes_ / (reduced_value + TinyReal);
}
//=================================================================================================//
} // namespace SPH
