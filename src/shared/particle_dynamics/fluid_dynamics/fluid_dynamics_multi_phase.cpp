#include "fluid_dynamics_multi_phase.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
ViscousAccelerationMultiPhase::ViscousAccelerationMultiPhase(BaseContactRelation &contact_relation)
    : BaseViscousAcceleration<MultiPhaseContactData>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real mu_k = DynamicCast<Fluid>(this, &contact_particles_[k]->getBaseMaterial())->ReferenceViscosity();
        contact_mu_.push_back(Real(2) * (mu_ * mu_k) / (mu_ + mu_k));
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
    }
}
//=================================================================================================//
void ViscousAccelerationMultiPhase::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real contact_mu_k = this->contact_mu_[k];
        StdLargeVec<Vecd> &vel_k = *(this->contact_vel_[k]);
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd vel_derivative = (this->vel_[index_i] - vel_k[index_j]) /
                                  (contact_neighborhood.r_ij_[n] + 0.01 * this->smoothing_length_);
            acceleration += 2.0 * contact_mu_k * vel_derivative * contact_neighborhood.dW_ijV_j_[n];
        }
    }
    acc_prior_[index_i] += acceleration / this->rho_[index_i];
}
//=================================================================================================//
MultiPhaseColorFunctionGradient::
    MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), MultiPhaseData(contact_relation),
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
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//