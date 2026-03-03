
#include "common_weakly_compressible_FVM_classes.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
WCAcousticTimeStepSizeInFVM::WCAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real min_distance_between_nodes, Real acousticCFL)
    : AcousticTimeStep(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      min_distance_between_nodes_(min_distance_between_nodes), acousticCFL_(acousticCFL) {};
//=================================================================================================//
Real WCAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
{
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ / Dimensions * min_distance_between_nodes_ / (reduced_value + TinyReal);
}
//=================================================================================================//
BaseForceFromFluidInFVM::BaseForceFromFluidInFVM(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      force_from_fluid_(nullptr) {};
//=================================================================================================//
ViscousForceFromFluidInFVM::ViscousForceFromFluidInFVM(BaseInnerRelation &inner_relation, StdVec<StdVec<size_t>> each_boundary_type_contact_real_index)
    : BaseForceFromFluidInFVM(inner_relation),
      viscosity_(DynamicCast<Viscosity>(this, particles_->getBaseMaterial())),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      mu_(viscosity_.ReferenceViscosity()),
      each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index)
{
    force_from_fluid_ = particles_->registerStateVariableData<Vecd>("ViscousForceOnSolid");
};
//=================================================================================================//
void ViscousForceFromFluidInFVM::interaction(size_t index_i, Real dt)
{
    for (size_t real_particle_num = 0; real_particle_num != each_boundary_type_contact_real_index_[3].size(); ++real_particle_num)
    {
        Vecd force = Vecd::Zero();
        const Vecd &vel_i = vel_[index_i];
        if (index_i == each_boundary_type_contact_real_index_[3][real_particle_num])
        {
            Real Vol_i = Vol_[index_i];
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];

            Vecd vel_j = -vel_i;
            size_t index_j = inner_neighborhood.j_[2];
            Vecd vel_derivative = (vel_j - vel_i) / (inner_neighborhood.r_ij_[2] + TinyReal);
            force += 2.0 * mu_ * vel_derivative * Vol_i * inner_neighborhood.dW_ij_[2] * Vol_[index_j];
            force_from_fluid_[index_i] = force;
        }
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
