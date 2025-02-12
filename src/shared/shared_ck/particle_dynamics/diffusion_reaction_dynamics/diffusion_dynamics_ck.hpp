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
template <class DiffusionType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    DiffusionRelaxationCK(SPHBody &sph_body, AbstractDiffusion &abstract_diffusion)
    : diffusions_(this->getConcreteDiffusions(abstract_diffusion)),
      dv_diffusion_species_array_(this->getDiffusionSpecies(sph_body.getBaseParticles())),
      dv_gradient_species_array_(this->getGradientSpecies(sph_body.getBaseParticles())),
      dv_diffusion_dt_array_(this->getSpeciesChangeRates(sph_body.getBaseParticles())),
      smoothing_length_sq_(pow(sph_body.getSPHAdaptation().ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <class DiffusionType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType>::DiffusionRelaxationCK(SPHBody &sph_body)
    : DiffusionRelaxationCK(
          sph_body, DynamicCast<AbstractDiffusion>(this, sph_body.getBaseMaterial())) {}
//=================================================================================================//
template <class DiffusionType>
StdVec<DiffusionType *> DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
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
template <class DiffusionType>
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    getDiffusionSpecies(BaseParticles &particles)
{
    StdVec<DiscreteVariable<Real> *> diffusion_species;
    for (auto &diffusion : diffusions_)
    {
        std::string species_name = diffusion->DiffusionSpeciesName();
        diffusion_species.push_back(
            particles.template registerStateVariableOnly<Real>(species_name));
        particles.template addVariableToSort<Real>(species_name);
        particles.template addVariableToWrite<Real>(species_name);
    }
    return diffusion_species;
}
//=================================================================================================//
template <class DiffusionType>
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    getGradientSpecies(BaseParticles &particles)
{
    StdVec<DiscreteVariable<Real> *> gradient_species;
    for (auto &diffusion : diffusions_)
    {
        std::string species_name = diffusion->GradientSpeciesName();
        gradient_species.push_back(
            particles.template registerStateVariableOnly<Real>(species_name));
        particles.template addVariableToSort<Real>(species_name);
        particles.template addVariableToWrite<Real>(species_name);
    }
    return gradient_species;
}
//=================================================================================================//
template <class DiffusionType>
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    getSpeciesChangeRates(BaseParticles &particles)
{
    StdVec<DiscreteVariable<Real> *> diffusion_dt;
    for (auto &diffusion : diffusions_)
    {
        std::string species_name = diffusion->DiffusionSpeciesName();
        diffusion_dt.push_back(
            particles.template registerStateVariableOnly<Real>(species_name + "ChangeRate"));
    }
    return diffusion_dt;
}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType>
void DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        diffusion_dt_[m][index_i] = 0;
    }
}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()) {}
//=================================================================================================//
template <class DiffusionType>
void DiffusionRelaxationCK<ForwardEuler, DiffusionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < number_of_species_; ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta, DiffusionType>::
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args)
    : BaseDynamicsType(sph_body, std::forward<Args>(args)...),
      dv_diffusion_species_array_s_(this->getIntermediateSpecies(sph_body.getBaseParticles())) {}
//=================================================================================================//
template <class DiffusionType>
StdVec<DiscreteVariable<Real> *> DiffusionRelaxationCK<RungeKutta, DiffusionType>::
    getIntermediateSpecies(BaseParticles &particles)
{
    StdVec<DiscreteVariable<Real> *> diffusion_species_s;
    for (auto &diffusion : this->diffusions_)
    {
        std::string species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_s.push_back(
            particles.template registerStateVariableOnly<Real>(species_name + "Intermediate"));
    }
    return diffusion_species_s;
}
//=================================================================================================//
template <class DiffusionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType>::
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args)
    : BaseDynamicsType(sph_body, std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InitializeKernel(ex_policy, encloser),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType>
void DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    BaseDynamicsType::InitializeKernel::initialize(index_i, dt);

    for (size_t m = 0; m < this->number_of_species_; ++m)
    {
        diffusion_species_s_[m][index_i] = diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionType>
template <typename... Args>
DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType>::
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args)
    : BaseDynamicsType(sph_body, std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::UpdateKernel(ex_policy, encloser),
      diffusion_species_s_(encloser.dv_diffusion_species_array_s_.DelegatedDataArray(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType>
void DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    BaseDynamicsType::UpdateKernel::update(index_i, dt);
    for (size_t m = 0; m < this->number_of_species_; ++m)
    {
        this->diffusion_species_[m][index_i] = 0.5 * diffusion_species_s_[m][index_i] +
                                               0.5 * this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType>
template <typename... Args>
DiffusionRelaxationCK<TimeSteppingType, DiffusionType, KernelGradientType>::
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args)
    : BaseRelaxation(sph_body, std::forward<Args>(args)...),
      kernel_gradient_(&sph_body.getBaseParticles()),
      ca_inter_particle_diffusion_coeff_(this->diffusions_),
      dv_Vol_(sph_body.getBaseParticles().template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<TimeSteppingType, DiffusionType, KernelGradientType>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : inter_particle_diffusion_coeff_(encloser.ca_inter_particle_diffusion_coeff_.DelegatedData(ex_policy)),
      diffusion_species_(encloser.dv_diffusion_species_array_.DelegatedDataArray(ex_policy)),
      gradient_species_(encloser.dv_gradient_species_array_.DelegatedDataArray(ex_policy)),
      diffusion_dt_(encloser.dv_diffusion_dt_array_.DelegatedDataArray(ex_policy)),
      number_of_species_(encloser.diffusions_.size()),
      gradient_(ex_policy, encloser.kernel_gradient_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
template <typename... Args>
DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    DiffusionRelaxationCK(Relation<Inner<Parameters...>> &relation, Args &&...args)
    : BaseRelaxation(relation.getSPHBody(), std::forward<Args>(args)...),
      Interaction<Inner<Parameters...>>(relation) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseRelaxation::InteractKernel(ex_policy, encloser),
      BaseInteraction::InteractKernel(ex_policy, encloser) {}
//=================================================================================================//
template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
void DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (size_t m = 0; m < this->number_of_species_; ++m)
    {
        Real d_species = 0.0;
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            Real dW_ijV_j = this->dW_ij(index_i, index_j) * this->Vol_[index_j];
            Vecd e_ij = this->e_ij(index_i, index_j);
            Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

            Real surface_area_ij = 2.0 * this->gradient_(index_i, index_j, dW_ijV_j, e_ij).dot(vec_r_ij) /
                                   (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);
            Real phi_ij = this->gradient_species_[m][index_i] - this->gradient_species_[m][index_j];
            d_species += this->inter_particle_diffusion_coeff_[m](index_i, index_j, e_ij) * phi_ij * surface_area_ij;
        }
        this->diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_HPP
