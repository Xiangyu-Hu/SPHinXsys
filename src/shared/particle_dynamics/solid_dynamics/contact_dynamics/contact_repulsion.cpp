#include "contact_repulsion.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
SelfContactForce::
    SelfContactForce(SelfSurfaceContactRelation &self_contact_relation)
    : LocalDynamics(self_contact_relation.getSPHBody()),
      SolidDataInner(self_contact_relation),
      solid_(particles_->solid_), mass_(particles_->mass_),
      self_repulsion_density_(*particles_->getVariableByName<Real>("SelfRepulsionDensity")),
      Vol_(particles_->Vol_), force_prior_(particles_->force_prior_),
      vel_(particles_->vel_),
      contact_impedance_(solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness())) {}
//=================================================================================================//
ContactForce::ContactForce(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      repulsion_density_(*particles_->getVariableByName<Real>("RepulsionDensity")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      force_prior_(particles_->force_prior_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("RepulsionDensity"));
    }
}
//=================================================================================================//
ContactForceFromWall::ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactWithWallData(solid_body_contact_relation), solid_(particles_->solid_),
      repulsion_density_(*particles_->getVariableByName<Real>("RepulsionDensity")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      force_prior_(particles_->force_prior_) {}
//=================================================================================================//
ContactForceToWall::ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      force_prior_(particles_->force_prior_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("RepulsionDensity"));
    }
}
//=================================================================================================//
DynamicContactForceWithWall::
    DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_), force_prior_(particles_->force_prior_),
      penalty_strength_(penalty_strength)
{
    impedance_ = solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness());
    reference_pressure_ = solid_.ReferenceDensity() * solid_.ContactStiffness();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
        contact_n_.push_back(&(contact_particles_[k]->n_));
    }
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
