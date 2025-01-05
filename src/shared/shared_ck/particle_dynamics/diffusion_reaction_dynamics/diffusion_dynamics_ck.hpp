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
DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::
    DiffusionRelaxationCK(DynamicsIdentifier &identifier) : BaseInteractionType(identifier)
{
    getDiffusions();

    for (auto &diffusion : diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name));
        this->particles_->template addVariableToSort<Real>(diffusion_species_name);
        this->particles_->template addVariableToWrite<Real>(diffusion_species_name);
        diffusion_dt_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name + "ChangeRate"));

        std::string gradient_species_name = diffusion->GradientSpeciesName();
        gradient_species_.push_back(this->particles_->template registerStateVariable<Real>(gradient_species_name));
        this->particles_->template addVariableToSort<Real>(gradient_species_name);
        this->particles_->template addVariableToWrite<Real>(gradient_species_name);
    }

    dv_diffusion_species_array_ = VariableArray<DiscreteVariable<Real>>("DiffusionSpecies", diffusion_species_);
    dv_gradient_species_array_ = VariableArray<DiscreteVariable<Real>>("GradientSpecies", gradient_species_);
    dv_diffusion_dt_array_ = VariableArray<DiscreteVariable<Real>>("DiffusionChangeRate", diffusion_dt_);
}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::getDiffusions()
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
DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_dt_(encloser.dv_diffusion_species_array_.DelegatedVariableArrayData(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::
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
DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedVariableArrayData(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedVariableArrayData(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType, class BaseInteractionType>
void DiffusionRelaxationCK<Base, DiffusionType, BaseInteractionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < diffusions_.size(); ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, typename... Parameters>
DiffusionRelaxationCK<Inner<OneLevel, DiffusionType, KernelGradientType, Parameters...>>::
    DiffusionRelaxationCK(Relation<Inner<Parameters...>> &inner_relation)
    : BaseInteraction(inner_relation),
      kernel_gradient_(this->particles_),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<OneLevel, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      inter_particle_diffusion_coeff_(encloser.dv_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedVariableArrayData(ex_policy)),
      gradient_species_(encloser.dv_gradient_species_array_.DelegatedVariableArrayData(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedVariableArrayData(ex_policy)),
      number_of_species_(encloser.diffusions_.size()),
      gradient_(ex_policy, encloser.kernel_gradient_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType, class KernelGradientType, typename... Parameters>
void DiffusionRelaxationCK<Inner<OneLevel, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        Real d_species = 0.0;
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            Real r_ij = this->vec_r_ij(index_i, index_j).norm();
            Vecd e_ij = this->e_ij(index_i, index_j);

            const Vecd grad_ijV_j = gradient_(index_i, index_j, dW_ijV_j, e_ij);
            Real surface_area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij;
            Real phi_ij = gradient_species_[m][index_i] - gradient_species_[m][index_j];
            d_species += inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * phi_ij * surface_area_ij;
        }
        diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_HPP