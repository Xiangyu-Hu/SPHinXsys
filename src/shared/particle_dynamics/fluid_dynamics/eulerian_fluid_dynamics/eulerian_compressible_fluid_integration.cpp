#include "eulerian_compressible_fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BaseIntegrationInCompressible::BaseIntegrationInCompressible(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation),
      compressible_fluid_(CompressibleFluid(1.0, 1.4)),
      Vol_(*particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      E_(*particles_->registerSharedVariable<Real>("TotalEnergy")),
      dE_dt_(*particles_->registerSharedVariable<Real>("TotalEnergyChangeRate")),
      dmass_dt_(*particles_->registerSharedVariable<Real>("MassChangeRate")),
      mom_(*particles_->registerSharedVariable<Vecd>("Momentum")),
      force_(*particles_->registerSharedVariable<Vecd>("Force")),
      force_prior_(*particles_->registerSharedVariable<Vecd>("ForcePrior")){};
//=================================================================================================//
CompressibleFluidInitialCondition::CompressibleFluidInitialCondition(SPHBody &sph_body)
    : FluidInitialCondition(sph_body),
      mom_(*particles_->getVariableDataByName<Vecd>("Momentum")),
      rho_(*particles_->getVariableDataByName<Real>("Density")),
      Vol_(*particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableDataByName<Real>("Mass")),
      p_(*particles_->getVariableDataByName<Real>("Pressure")),
      E_(*particles_->getVariableDataByName<Real>("TotalEnergy")) {}
//=================================================================================================//
EulerianCompressibleAcousticTimeStepSize::
    EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body)
    : AcousticTimeStepSize(sph_body),
      rho_(*particles_->getVariableDataByName<Real>("Density")),
      p_(*particles_->getVariableDataByName<Real>("Pressure")),
      vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
      compressible_fluid_(CompressibleFluid(1.0, 1.4)){};
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH