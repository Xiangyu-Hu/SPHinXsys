#include "contact_dynamics.h"

#ifdef max
#undef max
#endif

namespace SPH
{
//=========================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
SelfContactDensitySummation::
    SelfContactDensitySummation(SelfSurfaceContactRelation &self_contact_relation)
    : ContactDensityAccessor(self_contact_relation.base_particles_, "SelfContactDensity"),
      LocalDynamics(self_contact_relation.getSPHBody()),
      SolidDataInner(self_contact_relation),
      mass_(particles_->mass_)
{
    Real dp_1 = self_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    offset_W_ij_ = self_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(dp_1, ZeroVecd);
}
//=================================================================================================//
ContactDensitySummation::
    ContactDensitySummation(SurfaceContactRelation &solid_body_contact_relation)
    : ContactDensityAccessor(solid_body_contact_relation.base_particles_, "ContactDensity"),
      LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation), mass_(particles_->mass_),
      offset_W_ij_(StdVec<Real>(contact_configuration_.size(), 0.0))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
    }

    // we modify the default formulation by an offset, so that exactly touching bodies produce 0 initial force
    // subtract summation of the kernel function of 2 particles at 1 particle distance, and if the result is negative, we take 0
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    // dp1, dp2 half reference spacing
    Real dp_1 = solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing();
    // different resolution: distance = 0.5 * dp1 + 0.5 * dp2
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real dp_2 = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real distance = 0.5 * dp_1 + 0.5 * dp_2;
        offset_W_ij_[k] = solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()->W(distance, ZeroVecd);
    }
}
//=================================================================================================//
ShellContactDensity::ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation)
    : ContactDensityAccessor(solid_body_contact_relation.base_particles_, "ContactDensity"),
      LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation), solid_(particles_->solid_),
      kernel_(solid_body_contact_relation.getSPHBody().sph_adaptation_->getKernel()),
      particle_spacing_(solid_body_contact_relation.getSPHBody().sph_adaptation_->ReferenceSpacing())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real dp_k = solid_body_contact_relation.contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real average_spacing_k = 0.5 * particle_spacing_ + 0.5 * dp_k;
        Real h_ratio_k = particle_spacing_ / average_spacing_k;
        offset_W_ij_.push_back(kernel_->W(h_ratio_k, average_spacing_k, ZeroVecd));

        Real contact_max(0.0);
        for (int l = 0; l != 3; ++l)
        {
            Real temp = three_gaussian_points_[l] * average_spacing_k * 0.5 + average_spacing_k * 0.5;
            Real contact_temp = 2.0 * (kernel_->W(h_ratio_k, temp, ZeroVecd) - offset_W_ij_[k]) *
                                average_spacing_k * 0.5 * three_gaussian_weights_[l];
            contact_max += Dimensions == 2 ? contact_temp : contact_temp * Pi * temp;
        }
        /** a calibration factor to avoid particle penetration into shell structure */
        calibration_factor_.push_back(solid_.ReferenceDensity() / (contact_max + Eps));

        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
SelfContactForce::
    SelfContactForce(SelfSurfaceContactRelation &self_contact_relation)
    : LocalDynamics(self_contact_relation.getSPHBody()),
      SolidDataInner(self_contact_relation),
      solid_(particles_->solid_), mass_(particles_->mass_),
      self_contact_density_(*particles_->getVariableByName<Real>("SelfContactDensity")),
      Vol_(particles_->Vol_), acc_prior_(particles_->acc_prior_),
      vel_(particles_->vel_),
      contact_impedance_(solid_.ReferenceDensity() * sqrt(solid_.ContactStiffness())) {}
//=================================================================================================//
ContactForce::ContactForce(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      contact_density_(*particles_->getVariableByName<Real>("ContactDensity")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("ContactDensity"));
    }
}
//=================================================================================================//
ContactForceFromWall::ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactWithWallData(solid_body_contact_relation), solid_(particles_->solid_),
      contact_density_(*particles_->getVariableByName<Real>("ContactDensity")),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_) {}
//=================================================================================================//
ContactForceToWall::ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      acc_prior_(particles_->acc_prior_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&contact_particles_[k]->solid_);
        contact_contact_density_.push_back(contact_particles_[k]->getVariableByName<Real>("ContactDensity"));
    }
}
//=================================================================================================//
PairwiseFrictionFromWall::
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta)
    : LocalDynamics(contact_relation.getSPHBody()), ContactWithWallData(contact_relation),
      eta_(eta), Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        wall_vel_n_.push_back(&contact_particles_[k]->vel_);
        wall_n_.push_back(&contact_particles_[k]->n_);
    }
}
//=================================================================================================//
DynamicContactForceWithWall::
    DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength)
    : LocalDynamics(solid_body_contact_relation.getSPHBody()),
      ContactDynamicsData(solid_body_contact_relation),
      solid_(particles_->solid_),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      vel_(particles_->vel_), acc_prior_(particles_->acc_prior_),
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
