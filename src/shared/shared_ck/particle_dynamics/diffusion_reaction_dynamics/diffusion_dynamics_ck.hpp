/**
 * @file 	diffusion_dynamics.hpp
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_CK_HPP
#define DIFFUSION_DYNAMICS_CK_HPP

#include "diffusion_dynamics_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion)
    : BaseInteractionType(identifier),
      diffusions_(this->obtainConcreteDiffusions(*abstract_diffusion)),
      diffusion_species_names_(this->obtainDiffusionSpeciesNames(diffusions_)),
      gradient_species_names_(this->obtainGradientSpeciesNames(diffusions_)),
      dv_diffusion_species_array_(this->particles_->template registerStateVariables<Real>(
          diffusion_species_names_, "")),
      dv_gradient_species_array_(this->particles_->template registerStateVariables<Real>(
          gradient_species_names_, "")),
      dv_diffusion_dt_array_(this->particles_->template registerStateVariables<Real>(
          diffusion_species_names_, "ChangeRate")),
      ca_inverse_volume_capacity_(this->diffusions_)
{
    this->particles_->template addVariableToWrite<Real>(&dv_diffusion_species_array_);
    this->particles_->template addVariableToWrite<Real>(&dv_gradient_species_array_);
    this->particles_->template addEvolvingVariable<Real>(&dv_diffusion_species_array_);
    this->particles_->template addEvolvingVariable<Real>(&dv_gradient_species_array_);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class DynamicsIdentifier>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier)
    : DiffusionRelaxationCK(identifier, DynamicCast<AbstractDiffusion>(
                                            this, &identifier.getSPHBody().getBaseMaterial())) {}
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
    obtainDiffusionSpeciesNames(StdVec<DiffusionType *> &diffusions)
{
    StdVec<std::string> diffusion_species_names;
    for (auto &diffusion : diffusions)
    {
        diffusion_species_names.push_back(diffusion->DiffusionSpeciesName());
    }
    return diffusion_species_names;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<std::string> DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    obtainGradientSpeciesNames(StdVec<DiffusionType *> &diffusions)
{
    StdVec<std::string> gradient_species_names;
    for (auto &diffusion : diffusions)
    {
        gradient_species_names.push_back(diffusion->GradientSpeciesName());
    }
    return gradient_species_names;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType, typename... Args>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseInteractionType::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      gradient_species_(encloser.dv_gradient_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : cv1_(encloser.ca_inverse_volume_capacity_.DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
template <typename... Args>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...),
      kernel_correction_method_(this->particles_),
      ca_inter_particle_diffusion_coeff_(this->diffusions_),
      smoothing_length_sq_(pow(this->sph_adaptation_->ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_method_),
      inter_particle_diffusion_coeff_(encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class DiffusionType, class KernelCorrectionType, class... Parameters>
void DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        Real d_species = 0.0;
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

            Vecd surface_area = dW_ijV_j * (correction_(index_i) + correction_(index_j)) * e_ij;
            Vecd derivative = (this->gradient_species_[m][index_i] - this->gradient_species_[m][index_j]) *
                              vec_r_ij / (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);
            d_species += (inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * derivative).dot(surface_area);
        }
        this->diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelCorrectionType>
template <typename... Args>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...), kernel_correction_method_(this->particles_)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_dv_transfer_array_.push_back(
            contact_transfer_array_ptrs_keeper_.createPtr<DiscreteVariableArray<Real>>(
                this->particles_->template registerStateVariables<Real>(
                    this->diffusion_species_names_, "TransferWith" + this->sph_body_->getName())));
        contact_boundary_method_.push_back(
            boundary_ptrs_keeper_.template createPtr<BoundaryType<DiffusionType>>(
                *this, this->contact_particles_[k]));
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelCorrectionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_method_),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_transfer_(encloser.contact_dv_transfer_array_[contact_index]->DelegatedDataArray(ex_policy)),
      boundary_flux_(ex_policy, *encloser.contact_boundary_method_[contact_index]) {}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelCorrectionType>
void DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelCorrectionType>>::
    InteractKernel::interact(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        contact_transfer_[m][index_i] = 0.0;
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->contact_Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

            Vecd surface_area = 2.0 * dW_ijV_j * correction_(index_i) * e_ij;
            contact_transfer_[m][index_i] += boundary_flux_(m, index_i, index_j, e_ij, vec_r_ij).dot(surface_area);
        }
        this->diffusion_dt_[m][index_i] += contact_transfer_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < number_of_species_; ++m)
    {
        diffusion_dt_[m][index_i] = 0;
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < number_of_species_; ++m)
    {
        diffusion_species_[m][index_i] += this->cv1_[m](index_i) * dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->particles_->template registerStateVariables<Real>(
          this->diffusion_species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    InitializeKernel::initialize(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::InitializeKernel::initialize(index_i, dt);

    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        diffusion_species_s_[m][index_i] = diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->particles_->template getVariablesByName<Real>(
          this->diffusion_species_names_, "Intermediate")) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    BaseDynamicsType::UpdateKernel::update(index_i, dt);
    for (UnsignedInt m = 0; m < this->number_of_species_; ++m)
    {
        this->diffusion_species_[m][index_i] = 0.5 * diffusion_species_s_[m][index_i] +
                                               0.5 * this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType>
template <class DiffusionDynamics>
Dirichlet<DiffusionType>::Dirichlet(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles)
    : smoothing_length_sq_(
          pow(diffusion_dynamics.getSPHAdaptation().ReferenceSmoothingLength(), 2)),
      dv_gradient_species_array_(diffusion_dynamics.dvGradientSpeciesArray()),
      contact_dv_gradient_species_array_(contact_particles->template registerStateVariables<Real>(
          diffusion_dynamics.getGradientSpeciesNames(), "")),
      ca_inter_particle_diffusion_coeff_(diffusion_dynamics.getDiffusions())
{
    contact_particles->template addVariableToWrite<Real>(&contact_dv_gradient_species_array_);
}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
Dirichlet<DiffusionType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : smoothing_length_sq_(encloser.smoothing_length_sq_),
      gradient_species_(
          encloser.dv_gradient_species_array_.DelegatedDataArray(ex_policy)),
      contact_gradient_species_(
          encloser.contact_dv_gradient_species_array_.DelegatedDataArray(ex_policy)),
      inter_particle_diffusion_coeff_(
          encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)){};
//=================================================================================================//
template <class DiffusionType>
Vecd Dirichlet<DiffusionType>::ComputingKernel::operator()(
    UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j,
    const Vecd &e_ij, const Vecd &vec_r_ij)
{
    Real phi_ij = gradient_species_[m][index_i] - contact_gradient_species_[m][index_j];
    return vec_r_ij * phi_ij * inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) /
           (vec_r_ij.squaredNorm() + 0.01 * smoothing_length_sq_);
}
//=================================================================================================//
template <class DiffusionType>
template <class DiffusionDynamics>
Neumann<DiffusionType>::Neumann(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles)
    : dv_contact_n_(contact_particles->getVariableByName<Vecd>("NormalDirection")),
      contact_dv_species_flux_array_(contact_particles->template getVariablesByName<Real>(
          diffusion_dynamics.getDiffusionSpeciesNames(), "Flux")) {}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
Neumann<DiffusionType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : contact_n_(encloser.dv_contact_n_->DelegatedData(ex_policy)),
      contact_species_flux_(
          encloser.contact_dv_species_flux_array_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType>
Vecd Neumann<DiffusionType>::ComputingKernel::operator()(
    UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j,
    const Vecd &e_ij, const Vecd &vec_r_ij)
{
    return -contact_species_flux_[m][index_j] * contact_n_[index_j];
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_HPP
