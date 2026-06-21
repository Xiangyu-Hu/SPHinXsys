/**
 * @file 	diffusion_dynamics.hpp
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_CK_HPP
#define DIFFUSION_DYNAMICS_CK_HPP

#include "diffusion_dynamics_ck.h"

#include "base_body.hpp"

namespace SPH
{
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion)
    : BaseInteractionType(identifier),
      diffusions_(this->obtainConcreteDiffusions(*abstract_diffusion)),
      species_names_(this->obtainSpeciesNames(diffusions_)),
      dv_species_(registerSpecies(this->particles_, species_names_)),
      dv_species_dt_(registerSpecies(this->particles_, species_names_, "ChangeRate"))
{
    this->particles_->template addEvolvingVariable<Real>(dv_species_);
    this->particles_->template addVariableToWrite<Real>(dv_species_);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
DiscreteVariable<Real> *DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::registerSpecies(
    BaseParticles *particles, StdVec<std::string> &species_names, std::string suffix)
{
    if (species_names.size() == 1)
    {
        return particles->template registerStateVariable<Real>(species_names[0] + suffix);
    }
    return particles->template registerStateVariable<Real>("Species" + suffix, species_names);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
DiscreteVariable<Real> *DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::getSpeciesByName(
    BaseParticles *particles, StdVec<std::string> &species_names, std::string suffix)
{
    if (species_names.size() == 1)
    {
        return particles->template getVariableByName<Real>(species_names[0] + suffix);
    }
    return particles->template getVariableByName<Real>("Species" + suffix);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier)
    : DiffusionRelaxationCK(
          identifier, &identifier.getSPHBody().template getMaterialProperty<AbstractDiffusion>()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<DiffusionType *> DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainConcreteDiffusions(AbstractDiffusion &abstract_diffusion)
{
    StdVec<AbstractDiffusion *> all_diffusions = abstract_diffusion.AllDiffusions();
    StdVec<DiffusionType *> diffusions;
    for (auto &diffusion : all_diffusions)
    {
        diffusions.push_back(DynamicCast<DiffusionType>(this, diffusion));
    }
    return diffusions;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<std::string> DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainSpeciesNames(StdVec<DiffusionType *> &diffusions)
{
    StdVec<std::string> diffusion_species_names;
    for (auto &diffusion : diffusions)
    {
        diffusion_species_names.push_back(diffusion->DiffusionSpeciesName());
    }
    return diffusion_species_names;
}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
template <typename... Args>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...),
      kernel_correction_method_(this->particles_),
      cka_inter_particle_diffusion_coeff_(this->diffusions_),
      smoothing_length_sq_(pow(this->sph_adaptation_->ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_method_),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      species_dt_(encloser.dv_species_dt_->DelegatedMultiEntryView(ex_policy)),
      inter_particle_diffusion_coeff_(encloser.cka_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
void DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

        Vecd surface_area = dW_ijV_j * (correction_(index_i) + correction_(index_j)) * e_ij;
        Vecd rij_1 = vec_r_ij / (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);

        for (UnsignedInt m = 0; m < this->species_dt_.Width(); ++m)
        {
            Vecd derivative = (species_[index_i][m] - species_[index_j][m]) * rij_1;
            species_dt_[index_i][m] +=
                (inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * derivative).dot(surface_area);
        }
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType,
          class KernelCorrectionType, class... Parameters>
template <typename... Args>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType, Parameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...), kernel_correction_method_(this->particles_)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_dv_transfer_.push_back(
            this->registerSpecies( // registered in source body particles
                this->particles_, this->species_names_, "TransferWith" + this->sph_body_->Name()));
        contact_boundary_method_.push_back(
            boundaries_keeper_.template createPtr<BoundaryType<DiffusionType>>(
                *this, this->contact_particles_[k]));
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType,
          class KernelCorrectionType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_method_),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      species_dt_(encloser.dv_species_dt_->DelegatedMultiEntryView(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_transfer_(encloser.contact_dv_transfer_[contact_index]->DelegatedMultiEntryView(ex_policy)),
      boundary_flux_(ex_policy, *encloser.contact_boundary_method_[contact_index]) {}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType,
          class KernelCorrectionType, class... Parameters>
void DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < contact_transfer_.Width(); ++m)
    {
        contact_transfer_[index_i][m] = 0.0;
    }

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
        Vecd surface_area = 2.0 * dW_ijV_j * correction_(index_i) * e_ij;

        for (UnsignedInt m = 0; m < contact_transfer_.Width(); ++m)
        {
            Real flux = boundary_flux_(m, index_i, index_j, e_ij, vec_r_ij).dot(surface_area);
            contact_transfer_[index_i][m] += flux;
            species_dt_[index_i][m] += flux;
        }
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...), cka_inverse_volume_capacity_(this->diffusions_) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : species_dt_(encloser.dv_species_dt_->DelegatedMultiEntryView(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < species_dt_.Width(); ++m)
    {
        species_dt_[index_i][m] = 0.0;
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      species_dt_(encloser.dv_species_dt_->DelegatedMultiEntryView(ex_policy)),
      cv1_(encloser.cka_inverse_volume_capacity_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < species_.Width(); ++m)
    {
        species_[index_i][m] += cv1_[m](index_i) * dt * species_dt_[index_i][m];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_species_s_(this->registerSpecies(this->particles_, this->species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      species_s_(encloser.dv_species_s_->DelegatedMultiEntryView(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::InitializeKernel::initialize(index_i, dt);

    for (UnsignedInt m = 0; m < this->species_.Width(); ++m)
    {
        species_s_[index_i][m] = species_[index_i][m];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_species_s_(this->getSpeciesByName(this->particles_, this->species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      species_s_(encloser.dv_species_s_->DelegatedMultiEntryView(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::UpdateKernel::update(index_i, dt);
    for (UnsignedInt m = 0; m < this->species_.Width(); ++m)
    {
        species_[index_i][m] = 0.5 * species_s_[index_i][m] + 0.5 * species_[index_i][m];
    }
}
//=================================================================================================//
template <class DiffusionType>
template <class DiffusionDynamics>
Dirichlet<DiffusionType>::Dirichlet(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles)
    : smoothing_length_sq_(pow(diffusion_dynamics.getSPHAdaptation().ReferenceSmoothingLength(), 2)),
      dv_species_(diffusion_dynamics.getSpeciesByName(
          &diffusion_dynamics.getBaseParticles(), diffusion_dynamics.getSpeciesNames())),
      dv_contact_species_(diffusion_dynamics.registerSpecies(
          contact_particles, diffusion_dynamics.getSpeciesNames())),
      cka_inter_particle_diffusion_coeff_(diffusion_dynamics.getDiffusions())
{
    contact_particles->template addVariableToWrite<Real>(dv_contact_species_);
}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
Dirichlet<DiffusionType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : smoothing_length_sq_(encloser.smoothing_length_sq_),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      contact_species_(encloser.dv_contact_species_->DelegatedMultiEntryView(ex_policy)),
      inter_particle_diffusion_coeff_(
          encloser.cka_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)){};
//=================================================================================================//
template <class DiffusionType>
Vecd Dirichlet<DiffusionType>::ComputingKernel::operator()(
    UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j, const Vecd &e_ij, const Vecd &vec_r_ij)
{
    Real species_diff = species_[index_i][m] - contact_species_[index_j][m];
    return vec_r_ij * species_diff * inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) /
           (vec_r_ij.squaredNorm() + 0.01 * smoothing_length_sq_);
}
//=================================================================================================//
template <class DiffusionType>
template <class DiffusionDynamics>
Neumann<DiffusionType>::Neumann(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles)
    : dv_contact_n_(contact_particles->getVariableByName<Vecd>("NormalDirection")),
      dv_contact_species_flux_(diffusion_dynamics.registerSpecies(
          contact_particles, diffusion_dynamics.getSpeciesNames(), "Flux")) {}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
Neumann<DiffusionType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : contact_n_(encloser.dv_contact_n_->DelegatedData(ex_policy)),
      contact_species_flux_(encloser.dv_contact_species_flux_->DelegatedMultiEntryView(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType>
Vecd Neumann<DiffusionType>::ComputingKernel::operator()(
    UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j,
    const Vecd &e_ij, const Vecd &vec_r_ij)
{
    return -contact_species_flux_[index_j][m] * contact_n_[index_j];
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_HPP
