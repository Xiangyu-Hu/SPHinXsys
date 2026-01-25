#ifndef REPULSION_FORCE_HPP
#define REPULSION_FORCE_HPP

#include "base_particles.hpp"
#include "repulsion_force.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <typename... Parameters>
RepulsionForceCK<Base, Contact<Parameters...>>::
    RepulsionForceCK(Contact<Parameters...> &contact_relation, Real numerical_damping)
    : Interaction<Contact<Parameters...>>(contact_relation),
      ForcePriorCK(this->particles_, "RepulsionForce"),
      solid_contact_(DynamicCast<SolidContact>(this, this->particles_->getBaseMaterial())),
      numerical_damping_(numerical_damping),
      stiffness_(solid_contact_.ContactStiffness()),
      impedance_(sqrt(solid_contact_.ContactReferenceDensity() * stiffness_)),
      dv_repulsion_factor_(this->particles_->template getVariableByName<Real>("RepulsionFactor")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_repulsion_force_(ForcePriorCK::getCurrentForce()) {}
//=================================================================================================//
template <typename... Parameters>
template <typename... Args>
RepulsionForceCK<Contact<WithUpdate, Parameters...>>::
    RepulsionForceCK(Contact<Parameters...> &contact_relation, Args &&...args)
    : BaseInteractionType(contact_relation, std::forward<Args>(args)...),
      dv_n_(this->particles_->template getVariableByName<Vecd>("NormalDirection"))
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Solid &solid = DynamicCast<Solid>(this, this->contact_bodies_[k]->getBaseMaterial());
        contact_stiffness_.push_back(solid.ContactStiffness());
        contact_impedance_.push_back(sqrt(solid.ReferenceDensity() * solid.ContactStiffness()));
        dv_contact_repulsion_factor_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("RepulsionFactor"));
        dv_contact_vel_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
        dv_contact_n_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>("NormalDirection"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RepulsionForceCK<Contact<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteractionType::InteractKernel(ex_policy, encloser, contact_index),
      numerical_damping_(encloser.numerical_damping_),
      stiffness_ave_(2.0 * encloser.stiffness_ * encloser.contact_stiffness_[contact_index] /
                     (encloser.stiffness_ + encloser.contact_stiffness_[contact_index])),
      impedance_(encloser.impedance_), contact_impedance_(encloser.contact_impedance_[contact_index]),
      repulsion_factor_(encloser.dv_repulsion_factor_->DelegatedData(ex_policy)),
      contact_repulsion_factor_(encloser.dv_contact_repulsion_factor_[contact_index]->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      contact_vel_(encloser.dv_contact_vel_[contact_index]->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      contact_n_(encloser.dv_contact_n_[contact_index]->DelegatedData(ex_policy)),
      repulsion_force_(encloser.dv_repulsion_force_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void RepulsionForceCK<Contact<WithUpdate, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd n_ave = (impedance_ * n_[index_i] - contact_impedance_ * contact_n_[index_j]) /
                     (impedance_ + contact_impedance_);

        Real p_star = 0.5 * stiffness_ave_ * (repulsion_factor_[index_i] + contact_repulsion_factor_[index_j]);
        Real impedance_p = 0.5 * numerical_damping_ * (impedance_ + contact_impedance_) *
                           (vel_[index_i] - contact_vel_[index_j]).dot(-e_ij);
        // contact force to mimic pressure but pointing to the average normal direction
        force -= 2.0 * (p_star + impedance_p) * n_ave * n_ave.dot(e_ij) *
                 this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
    }
    repulsion_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <typename... Parameters>
template <typename... Args>
RepulsionForceCK<Contact<WithUpdate, Wall, Parameters...>>::
    RepulsionForceCK(Contact<Parameters...> &contact_relation, Args &&...args)
    : BaseInteractionType(contact_relation, std::forward<Args>(args)...),
      Interaction<Wall>(contact_relation) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RepulsionForceCK<Contact<WithUpdate, Wall, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteractionType::InteractKernel(ex_policy, encloser, contact_index),
      numerical_damping_(encloser.numerical_damping_),
      stiffness_(encloser.stiffness_), impedance_(encloser.impedance_),
      repulsion_factor_(encloser.dv_repulsion_factor_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      contact_vel_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedData(ex_policy)),
      contact_n_(encloser.dv_wall_n_[contact_index]->DelegatedData(ex_policy)),
      repulsion_force_(encloser.dv_repulsion_force_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void RepulsionForceCK<Contact<WithUpdate, Wall, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);

        Real p = stiffness_ * repulsion_factor_[index_i];
        Real impedance_p = 0.5 * numerical_damping_ * impedance_ *
                           (vel_[index_i] - contact_vel_[index_j]).dot(-e_ij);
        Vecd projection = contact_n_[index_j] * contact_n_[index_j].dot(e_ij);
        // contact force to mimic pressure but pointing to the wall normal direction
        force -= 2.0 * (p + impedance_p) * projection *
                 this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
    }
    repulsion_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FORCE_HPP