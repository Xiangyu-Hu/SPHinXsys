#include "general_reduce_ck.h"

namespace SPH
{
//=================================================================================================//
TotalKineticEnergyCK::TotalKineticEnergyCK(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity"))
{
    quantity_name_ = "TotalKineticEnergy";
}
//=================================================================================================//
TotalMechanicalEnergyCK::TotalMechanicalEnergyCK(SPHBody &sph_body, Gravity &gravity)
    : TotalKineticEnergyCK(sph_body), gravity_(gravity),
      dv_pos_(particles_->getVariableByName<Vecd>("Position"))
{
    quantity_name_ = "TotalMechanicalEnergy";
}
//=================================================================================================//
} // namespace SPH
