
#include "common_weakly_compressible_FVM_classes.h"

namespace SPH
{
//=================================================================================================//
WCAcousticTimeStepSizeInFVM::WCAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real max_distance_between_nodes, Real acousticCFL)
    : AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
      vel_(particles_->vel_), fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      max_distance_between_nodes_(max_distance_between_nodes), acousticCFL_(acousticCFL){};
//=================================================================================================//
Real WCAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
{
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ / Dimensions * max_distance_between_nodes_ / (reduced_value + TinyReal);
}
//=================================================================================================//
BaseForceFromFluidInFVM::BaseForceFromFluidInFVM(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), fluid_dynamics::FluidDataInner(inner_relation), Vol_(particles_->Vol_){};
//=================================================================================================//
ViscousForceFromFluidInFVM::ViscousForceFromFluidInFVM(BaseInnerRelation &inner_relation, vector<vector<size_t>> each_boundary_type_contact_real_index)
    : BaseForceFromFluidInFVM(inner_relation), fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      vel_(particles_->vel_), mu_(fluid_.ReferenceViscosity()), each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index)
{
    particles_->registerVariable(force_from_fluid_, "ViscousForceFromFluid");
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
            Vecd vel_derivative = (vel_j - vel_i) / (inner_neighborhood.r_ij_[2] + TinyReal);
            force += 2.0 * mu_ * vel_derivative * Vol_i * inner_neighborhood.dW_ijV_j_[2];
            force_from_fluid_[index_i] = force;
        }
    }
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//