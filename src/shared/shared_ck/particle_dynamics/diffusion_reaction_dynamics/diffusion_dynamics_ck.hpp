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
      diffusions_(this->getConcreteDiffusions(*abstract_diffusion)),
      dv_diffusion_species_array_(this->getDiffusionVariables(this->particles_, "")),
      dv_gradient_species_array_(this->getGradientVariables(this->particles_, "")),
      dv_diffusion_dt_array_(this->getDiffusionVariables(this->particles_, "ChangeRate"))
{
    this->particles_->template addVariableToWrite<Real>(&dv_diffusion_species_array_);
    this->particles_->template addVariableToWrite<Real>(&dv_gradient_species_array_);
    this->particles_->template addVariableToSort<Real>(&dv_diffusion_species_array_);
    this->particles_->template addVariableToSort<Real>(&dv_gradient_species_array_);
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
    getConcreteDiffusions(AbstractDiffusion &abstract_diffusion)
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
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    getDiffusionVariables(BaseParticles *particles, const std::string &suffix)
{
    StdVec<DiscreteVariable<Real> *> diffusion_variables;
    for (auto &diffusion : diffusions_)
    {
        std::string variable_name = diffusion->DiffusionSpeciesName() + suffix;
        diffusion_variables.push_back(
            particles->template registerStateVariableOnly<Real>(variable_name));
    }
    return diffusion_variables;
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<DiffusionType, BaseInteractionType>::
    getGradientVariables(BaseParticles *particles, const std::string &suffix)
{
    StdVec<DiscreteVariable<Real> *> diffusion_variables;
    for (auto &diffusion : diffusions_)
    {
        std::string variable_name = diffusion->GradientSpeciesName() + suffix;
        diffusion_variables.push_back(
            particles->template registerStateVariableOnly<Real>(variable_name));
    }
    return diffusion_variables;
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
template <class DiffusionType, class KernelGradientType, class... Parameters>
template <typename... Args>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...),
      kernel_gradient_(this->particles_),
      ca_inter_particle_diffusion_coeff_(this->diffusions_),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      smoothing_length_sq_(pow(this->sph_adaptation_->ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, class... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      gradient_(ex_policy, encloser.kernel_gradient_),
      inter_particle_diffusion_coeff_(encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, class... Parameters>
void DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>::
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

            Real surface_area_ij = 2.0 * gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(vec_r_ij) /
                                   (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);
            Real phi_ij = this->gradient_species_[m][index_i] - this->gradient_species_[m][index_j];
            d_species += inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * phi_ij * surface_area_ij;
        }
        this->diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
template <typename... Args>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseInteraction(std::forward<Args>(args)...)
{
    for (UnsignedInt k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_dv_transfer_array_.push_back(
            contact_transfer_array_ptrs_keeper_.createPtr<DiscreteVariableArray<Real>>(
                this->getDiffusionVariables(this->particles_, "TransferWith" + this->sph_body_.getName())));
        contact_kernel_gradient_method_.push_back(
            kernel_gradient_ptrs_keeper_.template createPtr<KernelGradientType>(
                this->particles_, this->contact_particles_[k]));
        contact_boundary_method_.push_back(
            boundary_ptrs_keeper_.template createPtr<BoundaryType<DiffusionType>>(
                *this, this->contact_particles_[k]));
    }
}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_transfer_(encloser.contact_dv_transfer_array_[contact_index]->DelegatedDataArray(ex_policy)),
      gradient_(ex_policy, *encloser.contact_kernel_gradient_method_[contact_index]),
      boundary_flux_(ex_policy, *encloser.contact_boundary_method_[contact_index]) {}
//=================================================================================================//
template <class DiffusionType, template <typename...> class BoundaryType, class KernelGradientType>
void DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>::
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

            Vecd surface_area_ij = 2.0 * gradient_(index_i, index_j, dW_ijV_j, e_ij);
            contact_transfer_[m][index_i] += boundary_flux_(m, index_i, index_j, e_ij, vec_r_ij).dot(surface_area_ij);
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
    : diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
void DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>::
    UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    for (UnsignedInt m = 0; m < number_of_species_; ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <template <typename...> class RelationType, class... InteractionParameters>
template <typename... Args>
DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->getDiffusionVariables(this->particles_, "Intermediate")) {}
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
      dv_diffusion_species_array_s_(this->getDiffusionVariables(this->particles_, "Intermediate")) {}
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
          pow(diffusion_dynamics.getSPHAdaptation()->ReferenceSmoothingLength(), 2)),
      dv_gradient_species_array_(diffusion_dynamics.dvGradientSpeciesArray()),
      contact_dv_gradient_species_array_(diffusion_dynamics.getGradientVariables(contact_particles, "")),
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
      contact_dv_species_flux_array_(diffusion_dynamics.getDiffusionVariables(contact_particles, "Flux")) {}
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
