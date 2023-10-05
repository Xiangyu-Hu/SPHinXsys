#include "multi_phase_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
MultiPhaseColorFunctionGradient::
    MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), MultiPhaseContactData(contact_relation),
      rho0_(sph_body_.base_material_->ReferenceDensity()), Vol_(particles_->Vol_),
      pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
      indicator_(*particles_->getVariableByName<int>("Indicator"))
{
    particles_->registerVariable(color_grad_, "ColorGradient");
    particles_->registerVariable(surface_norm_, "SurfaceNormal");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_rho0_.push_back(rho0_k);
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
void MultiPhaseColorFunctionGradient::interaction(size_t index_i, Real dt)
{
    Real Vol_i = Vol_[index_i];
    Vecd gradient = Vecd::Zero();
    if (indicator_[index_i])
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real rho0_k = contact_rho0_[k];
            StdLargeVec<Real> &contact_Vol_k = *(contact_Vol_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                /** Norm of interface.*/
                Real rho_ij = rho0_ / (rho0_ + rho0_k);
                Real area_ij = (Vol_i * Vol_i + contact_Vol_k[index_j] * contact_Vol_k[index_j]) *
                               contact_neighborhood.dW_ijV_j_[n] / contact_Vol_k[index_j];
                gradient += rho_ij * area_ij * contact_neighborhood.e_ij_[n] / Vol_i;
            }
        }
    }
    color_grad_[index_i] = gradient;
    surface_norm_[index_i] = gradient / (gradient.norm() + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
