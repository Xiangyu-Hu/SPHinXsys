#include "fluid_structure_interaction.hpp"

#include "viscosity.h"

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
BaseForceFromFluid::BaseForceFromFluid(BaseContactRelation &contact_relation, const std::string &force_name)
    : ForcePrior(contact_relation.getSPHBody(), force_name), DataDelegateContact(contact_relation),
      solid_(DynamicCast<Solid>(this, sph_body_->getBaseMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      force_from_fluid_(particles_->getVariableDataByName<Vecd>(force_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_fluids_.push_back(DynamicCast<Fluid>(this, &contact_particles_[k]->getBaseMaterial()));
    }
}
//=================================================================================================//
ViscousForceFromFluid::ViscousForceFromFluid(BaseContactRelation &contact_relation)
    : BaseForceFromFluid(contact_relation, "ViscousForceFromFluid"),
      vel_ave_(solid_.AverageVelocity(particles_))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Velocity"));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        Viscosity &viscosity_k = DynamicCast<Viscosity>(this, contact_particles_[k]->getBaseMaterial());
        mu_.push_back(viscosity_k.ReferenceViscosity());
        smoothing_length_.push_back(contact_bodies_[k]->getSPHAdaptation().ReferenceSmoothingLength());
    }
}
//=================================================================================================//
void ViscousForceFromFluid::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    /** Contact interaction. */
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real mu_k = mu_[k];
        Real smoothing_length_k = smoothing_length_[k];
        Vecd *vel_n_k = contact_vel_[k];
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];

            Vecd vel_derivative = 2.0 * (vel_ave_[index_i] - vel_n_k[index_j]) /
                                  (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);
            force += 2.0 * mu_k * vel_derivative * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }

    force_from_fluid_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
InitializeDisplacement::
    InitializeDisplacement(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos_temp_(particles_->registerStateVariableData<Vecd>("TemporaryPosition")) {}
//=================================================================================================//
void InitializeDisplacement::update(size_t index_i, Real dt)
{
    pos_temp_[index_i] = pos_[index_i];
}
//=================================================================================================//
UpdateAverageVelocityAndAcceleration::
    UpdateAverageVelocityAndAcceleration(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos_temp_(particles_->getVariableDataByName<Vecd>("TemporaryPosition")),
      vel_ave_(particles_->getVariableDataByName<Vecd>("AverageVelocity")),
      acc_ave_(particles_->getVariableDataByName<Vecd>("AverageAcceleration")) {}
//=================================================================================================//
void UpdateAverageVelocityAndAcceleration::update(size_t index_i, Real dt)
{
    Vecd updated_vel_ave = (pos_[index_i] - pos_temp_[index_i]) / (dt + Eps);
    acc_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
    vel_ave_[index_i] = updated_vel_ave;
}
//=================================================================================================//
AverageVelocityAndAcceleration::
    AverageVelocityAndAcceleration(SPHBody &sph_body)
    : initialize_displacement_(sph_body),
      update_averages_(sph_body) {}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
