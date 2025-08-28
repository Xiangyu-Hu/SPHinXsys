#include "viscous_dynamics.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
VorticityInner::VorticityInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      vorticity_(particles_->registerStateVariable<AngularVecd>("VorticityInner"))
{
    particles_->addVariableToWrite<AngularVecd>("VorticityInner");
}
//=================================================================================================//
void VorticityInner::interaction(size_t index_i, Real dt)
{
    AngularVecd vorticity = ZeroData<AngularVecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Vecd vel_diff = vel_[index_i] - vel_[index_j];
        vorticity += getCrossProduct(vel_diff, inner_neighborhood.e_ij_[n]) * inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }

    vorticity_[index_i] = vorticity;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
