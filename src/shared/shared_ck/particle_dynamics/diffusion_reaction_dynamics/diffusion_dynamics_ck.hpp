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
DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier), dv_diffusion_species_array_(nullptr),
      dv_gradient_species_array_(nullptr), dv_diffusion_dt_array_(nullptr),
      smoothing_length_sq_(pow(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2))
{
    getDiffusions();

    StdVec<DiscreteVariable<Real> *> diffusion_species;
    StdVec<DiscreteVariable<Real> *> diffusion_dt;
    StdVec<DiscreteVariable<Real> *> gradient_species;
    for (auto &diffusion : diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species.push_back(this->particles_->template registerStateVariableOnly<Real>(diffusion_species_name));
        this->particles_->template addVariableToSort<Real>(diffusion_species_name);
        this->particles_->template addVariableToWrite<Real>(diffusion_species_name);
        diffusion_dt.push_back(this->particles_->template registerStateVariableOnly<Real>(diffusion_species_name + "ChangeRate"));

        std::string gradient_species_name = diffusion->GradientSpeciesName();
        gradient_species.push_back(this->particles_->template registerStateVariableOnly<Real>(gradient_species_name));
        this->particles_->template addVariableToSort<Real>(gradient_species_name);
        this->particles_->template addVariableToWrite<Real>(gradient_species_name);
    }
    dv_diffusion_species_array_ = variable_array_ptrs_
                                      .template createPtr<DiscreteVariableArray<Real>>(diffusion_species);
    dv_gradient_species_array_ = variable_array_ptrs_
                                     .template createPtr<DiscreteVariableArray<Real>>(gradient_species);
    dv_diffusion_dt_array_ = variable_array_ptrs_
                                 .template createPtr<DiscreteVariableArray<Real>>(diffusion_dt);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::getDiffusions()
{
    AbstractDiffusion &abstract_diffusion = DynamicCast<AbstractDiffusion>(this, this->sph_body_.getBaseMaterial());
    StdVec<AbstractDiffusion *> all_diffusions = abstract_diffusion.AllDiffusions();
    for (auto &diffusion : all_diffusions)
    {
        diffusions_.push_back(DynamicCast<DiffusionType>(this, diffusion));
    }
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_dt_(encloser.dv_diffusion_species_array_->DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        diffusion_dt_[m][index_i] = 0;
    }
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_species_(encloser.dv_diffusion_species_array_->DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_->DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(nullptr)
{
    StdVec<DiscreteVariable<Real> *> diffusion_species_s;
    for (auto &diffusion : this->diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_s.push_back(
            this->particles_->template registerStateVariableOnly<Real>(diffusion_species_name + "Intermediate"));
    }
    dv_diffusion_species_array_s_ = this->variable_array_ptrs_
                                        .template createPtr<DiscreteVariableArray<Real>>(diffusion_species_s);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType, BaseInteractionType>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_->DelegatedDataArray(ex_policy)),
      diffusion_species_(encloser.dv_diffusion_species_array_->DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType, BaseInteractionType>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    BaseDynamicsType::InitializeKernel::initialize(index_i, dt);

    for (size_t m = 0; m < this->number_of_species_; ++m)
    {
        diffusion_species_s_[m][index_i] = this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(Args &&...args)
    : BaseDynamicsType(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType, BaseInteractionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_->DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType, BaseInteractionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    BaseDynamicsType::UpdateKernel::update(index_i, dt);
    for (size_t m = 0; m < this->number_of_species_; ++m)
    {
        this->diffusion_species_[m][index_i] = 0.5 * this->diffusion_species_s_[m][index_i] +
                                               0.5 * this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    DiffusionRelaxationCK(Relation<Inner<Parameters...>> &inner_relation)
    : BaseInteraction(inner_relation),
      kernel_gradient_(this->particles_),
      ca_inter_particle_diffusion_coeff_(this->diffusions_),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      inter_particle_diffusion_coeff_(encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      diffusion_species_(encloser.dv_diffusion_species_array_->DelegatedDataArray(ex_policy)),
      gradient_species_(encloser.dv_gradient_species_array_->DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_->DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()),
      gradient_(ex_policy, encloser.kernel_gradient_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
void DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        Real d_species = 0.0;
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

            Real surface_area_ij = 2.0 * gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(vec_r_ij) /
                                   (vec_r_ij.squaredNorm() + 0.01 * smoothing_length_sq_);
            Real phi_ij = gradient_species_[m][index_i] - gradient_species_[m][index_j];
            d_species += inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * phi_ij * surface_area_ij;
        }
        diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_HPP
