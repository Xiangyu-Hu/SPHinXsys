#include "contact_repulsion.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
RepulsionForce<Contact<Inner<>>>::
    RepulsionForce(BaseInnerRelation &self_contact_relation)
    : RepulsionForce<Base, DataDelegateInner>(self_contact_relation, "SelfRepulsionForce"),
      ForcePrior(particles_, "SelfRepulsionForce"),
      solid_(DynamicCast<Solid>(this, sph_body_.getBaseMaterial())),
      self_repulsion_factor_(*particles_->getVariableDataByName<Real>("SelfRepulsionFactor")),
      vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
      contact_impedance_(sqrt(solid_.ReferenceDensity() * solid_.ContactStiffness())) {}
//=================================================================================================//
void RepulsionForce<Contact<Inner<>>>::interaction(size_t index_i, Real dt)
{
    Real p_i = self_repulsion_factor_[index_i] * solid_.ContactStiffness();
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real p_star = 0.5 * (p_i + self_repulsion_factor_[index_j] * solid_.ContactStiffness());
        Real impedance_p = 0.5 * contact_impedance_ * (vel_[index_i] - vel_[index_j]).dot(-e_ij);
        // force to mimic pressure
        force -= 2.0 * (p_star + impedance_p) * e_ij * inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }
    repulsion_force_[index_i] = force * particles_->ParticleVolume(index_i);
}
//=================================================================================================//
RepulsionForce<Contact<>>::RepulsionForce(BaseContactRelation &solid_body_contact_relation)
    : RepulsionForce<Base, DataDelegateContact>(solid_body_contact_relation, "RepulsionForce"),
      ForcePrior(particles_, "RepulsionForce"),
      solid_(DynamicCast<Solid>(this, sph_body_.getBaseMaterial())),
      repulsion_factor_(*particles_->getVariableDataByName<Real>("RepulsionFactor"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&DynamicCast<Solid>(this, contact_bodies_[k]->getBaseMaterial()));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_repulsion_factor_.push_back(contact_particles_[k]->getVariableDataByName<Real>("RepulsionFactor"));
    }
}
//=================================================================================================//
void RepulsionForce<Contact<>>::interaction(size_t index_i, Real dt)
{
    Real sigma_i = repulsion_factor_[index_i];
    Vecd force = Vecd::Zero();
    Real contact_stiffness_i = solid_.ContactStiffness(index_i);

    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &contact_density_k = *(contact_repulsion_factor_[k]);
        StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);

        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];

            Real contact_stiffness_j = contact_solids_[k]->ContactStiffness(index_j);
            Real contact_stiffness_ave = 2 * contact_stiffness_i * contact_stiffness_j / (contact_stiffness_i + contact_stiffness_j);

            Vecd e_ij = contact_neighborhood.e_ij_[n];

            Real sigma_star = (sigma_i + contact_density_k[index_j]) * contact_stiffness_ave;
            // force due to pressure
            force -= sigma_star * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    repulsion_force_[index_i] = force * particles_->ParticleVolume(index_i);
}
//=================================================================================================//
RepulsionForce<Contact<Wall>>::RepulsionForce(BaseContactRelation &solid_body_contact_relation)
    : RepulsionForce<Base, DataDelegateContact>(solid_body_contact_relation, "RepulsionForce"),
      ForcePrior(particles_, "RepulsionForce"),
      solid_(DynamicCast<Solid>(this, sph_body_.getBaseMaterial())),
      repulsion_factor_(*particles_->getVariableDataByName<Real>("RepulsionFactor"))
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        contact_Vol_.push_back(this->contact_particles_[k]->template registerSharedVariable<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
void RepulsionForce<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real p_i = repulsion_factor_[index_i] * solid_.ContactStiffness(index_i);
    /** Contact interaction. */
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];

            // force due to pressure
            force -= 2.0 * p_i * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    repulsion_force_[index_i] = force * particles_->ParticleVolume(index_i);
}
//=================================================================================================//
RepulsionForce<Wall, Contact<>>::RepulsionForce(BaseContactRelation &solid_body_contact_relation)
    : RepulsionForce<Base, DataDelegateContact>(solid_body_contact_relation, "RepulsionForce"),
      ForcePrior(particles_, "RepulsionForce")
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_solids_.push_back(&DynamicCast<Solid>(this, contact_bodies_[k]->getBaseMaterial()));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_repulsion_factor_.push_back(contact_particles_[k]->getVariableDataByName<Real>("RepulsionFactor"));
    }
}
//=================================================================================================//
void RepulsionForce<Wall, Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
        StdLargeVec<Real> &contact_density_k = *(contact_repulsion_factor_[k]);
        Solid *solid_k = contact_solids_[k];

        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];

            Real p_star = contact_density_k[index_j] * solid_k->ContactStiffness(index_j);
            // force due to pressure
            force -= 2.0 * p_star * e_ij * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
        }
    }
    repulsion_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
