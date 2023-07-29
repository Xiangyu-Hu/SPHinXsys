#include "fluid_structure_interaction.h"

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
BaseForceFromFluid::BaseForceFromFluid(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), FSIContactData(contact_relation),
      Vol_(particles_->Vol_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_fluids_.push_back(DynamicCast<Fluid>(this, &contact_particles_[k]->getBaseMaterial()));
    }
}
//=================================================================================================//
ViscousForceFromFluid::ViscousForceFromFluid(BaseContactRelation &contact_relation)
    : BaseForceFromFluid(contact_relation), vel_ave_(*particles_->AverageVelocity())
{
    particles_->registerVariable(force_from_fluid_, "ViscousForceFromFluid");
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
        mu_.push_back(contact_fluids_[k]->ReferenceViscosity());
        smoothing_length_.push_back(contact_bodies_[k]->sph_adaptation_->ReferenceSmoothingLength());
    }
}
//=================================================================================================//
void TotalForceFromFluid::setupDynamics(Real dt)
{
    if (!force_from_fluid_dynamics_.checkNewlyUpdated())
    {
        force_from_fluid_dynamics_.exec();
    }
    force_from_fluid_dynamics_.setNotNewlyUpdated();
}
//=================================================================================================//
Vecd TotalForceFromFluid::reduce(size_t index_i, Real dt)
{
    return force_from_fluid_[index_i];
}
//=================================================================================================//
InitializeDisplacement::
    InitializeDisplacement(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp)
    : LocalDynamics(sph_body), ElasticSolidDataSimple(sph_body),
      pos_temp_(pos_temp), pos_(particles_->pos_) {}
//=================================================================================================//
void InitializeDisplacement::update(size_t index_i, Real dt)
{
    pos_temp_[index_i] = pos_[index_i];
}
//=================================================================================================//
UpdateAverageVelocityAndAcceleration::
    UpdateAverageVelocityAndAcceleration(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp)
    : LocalDynamics(sph_body), ElasticSolidDataSimple(sph_body),
      pos_temp_(pos_temp), pos_(particles_->pos_),
      vel_ave_(particles_->vel_ave_),
      acc_ave_(particles_->acc_ave_) {}
//=================================================================================================//
void UpdateAverageVelocityAndAcceleration::update(size_t index_i, Real dt)
{
    Vecd updated_vel_ave = (pos_[index_i] - pos_temp_[index_i]) / (dt + Eps);
    acc_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
    vel_ave_[index_i] = updated_vel_ave;
}
//=================================================================================================//
AverageVelocityAndAcceleration::
    AverageVelocityAndAcceleration(SolidBody &solid_body)
    : initialize_displacement_(solid_body, pos_temp_),
      update_averages_(solid_body, pos_temp_)
{
    solid_body.getBaseParticles().registerVariable(pos_temp_, "TemporaryPosition");
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
