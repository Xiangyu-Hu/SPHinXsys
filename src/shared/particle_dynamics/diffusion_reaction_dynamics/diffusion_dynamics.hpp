/**
 * @file 	diffusion_dynamics.hpp
 * @brief 	This is the particle dynamics applicable for all type bodies
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_HPP
#define DIFFUSION_DYNAMICS_HPP

#include "diffusion_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
template <class BodyRelationType>
DiffusionRelaxation<DataDelegationType, DiffusionType>::
    DiffusionRelaxation(BodyRelationType &body_relation)
    : LocalDynamics(body_relation.getSPHBody()), DataDelegationType(body_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure"))
{
    getDiffusions();

    for (auto &diffusion : diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name));
        this->particles_->template addEvolvingVariable<Real>(diffusion_species_name);
        this->particles_->template addVariableToWrite<Real>(diffusion_species_name);
        diffusion_dt_.push_back(this->particles_->template registerStateVariable<Real>(diffusion_species_name + "ChangeRate"));

        std::string gradient_species_name = diffusion->GradientSpeciesName();
        gradient_species_.push_back(this->particles_->template registerStateVariable<Real>(gradient_species_name));
        this->particles_->template addEvolvingVariable<Real>(gradient_species_name);
        this->particles_->template addVariableToWrite<Real>(gradient_species_name);
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void DiffusionRelaxation<DataDelegationType, DiffusionType>::getDiffusions()
{
    AbstractDiffusion &abstract_diffusion = DynamicCast<AbstractDiffusion>(this, this->sph_body_.getBaseMaterial());
    StdVec<AbstractDiffusion *> all_diffusions = abstract_diffusion.AllDiffusions();
    for (auto &diffusion : all_diffusions)
    {
        diffusions_.push_back(DynamicCast<DiffusionType>(this, diffusion));
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void DiffusionRelaxation<DataDelegationType, DiffusionType>::initialization(size_t index_i, Real dt)
{
    for (size_t m = 0; m < diffusions_.size(); ++m)
    {
        diffusion_dt_[m][index_i] = 0;
    }
}
//=================================================================================================//
template <class DataDelegationType, class DiffusionType>
void DiffusionRelaxation<DataDelegationType, DiffusionType>::update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < diffusions_.size(); ++m)
    {
        diffusion_species_[m][index_i] += dt * diffusion_dt_[m][index_i];
    }
}
//=================================================================================================//
template <class KernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<DataDelegateInner, DiffusionType>(std::forward<Args>(args)...),
      kernel_gradient_(this->particles_) {}
//=================================================================================================//
template <class KernelGradientType, class DiffusionType>
void DiffusionRelaxation<Inner<KernelGradientType>, DiffusionType>::interaction(size_t index_i, Real dt)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        auto diffusion_m = this->diffusions_[m];
        Real *gradient_species = this->gradient_species_[m];
        Real d_species = 0.0;
        Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
            Real r_ij_ = inner_neighborhood.r_ij_[n];
            Vecd &e_ij = inner_neighborhood.e_ij_[n];

            Real diff_coeff_ij = diffusion_m->getInterParticleDiffusionCoeff(index_i, index_j, e_ij);
            const Vecd &grad_ijV_j = this->kernel_gradient_(index_i, index_j, dW_ijV_j, e_ij);
            Real surface_area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            Real phi_ij = gradient_species[index_i] - gradient_species[index_j];
            d_species += diff_coeff_ij * phi_ij * surface_area_ij;
        }
        this->diffusion_dt_[m][index_i] += d_species;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<DataDelegateContact, DiffusionType>(
          std::forward<Args>(args)...)
{
    contact_transfer_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_kernel_gradients_.push_back(ContactKernelGradientType(this->particles_, contact_particles_k));
        contact_Vol_.push_back(contact_particles_k->template registerStateVariable<Real>("VolumetricMeasure"));

        std::string diffusion_direction = "From" + this->contact_bodies_[k]->getName();
        for (auto &diffusion : this->diffusions_)
        {
            std::string variable_name = diffusion->GradientSpeciesName() + "Transfer" + diffusion_direction;
            contact_transfer_[k].push_back(
                this->particles_->template registerStateVariable<Real>(variable_name));
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    resetContactTransfer(size_t index_i)
{
    for (size_t k = 0; k < this->contact_particles_.size(); ++k)
    {
        for (size_t m = 0; m < this->diffusions_.size(); ++m)
        {
            this->contact_transfer_[k][m][index_i] = 0.0;
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>::
    accumulateDiffusionRate(size_t index_i)
{
    for (size_t k = 0; k < this->contact_particles_.size(); ++k)
    {
        for (size_t m = 0; m < this->diffusions_.size(); ++m)
        {
            this->diffusion_dt_[m][index_i] += this->contact_transfer_[k][m][index_i];
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(std::forward<Args>(args)...)
{
    contact_gradient_species_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        for (auto &diffusion : this->diffusions_)
        {
            std::string gradient_species_name = diffusion->GradientSpeciesName();
            contact_gradient_species_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(gradient_species_name));
            contact_particles_k->template addVariableToWrite<Real>(gradient_species_name);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateDirichlet(size_t particle_i, size_t particle_j, Vecd &e_ij,
                                    Real surface_area_ij, const StdVec<Real *> &gradient_species_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real diff_coeff_ij =
            this->diffusions_[m]->getInterParticleDiffusionCoeff(particle_i, particle_i, e_ij);
        Real phi_ij = 2.0 * (this->gradient_species_[m][particle_i] - gradient_species_k[m][particle_j]);
        this->diffusion_dt_[m][particle_i] += diff_coeff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Dirichlet<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &gradient_species_k = this->contact_gradient_species_[k];
        Real *wall_Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij_ = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            getDiffusionChangeRateDirichlet(index_i, index_j, e_ij, area_ij, gradient_species_k);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(
          std::forward<Args>(args)...),
      n_(this->particles_->template getVariableDataByName<Vecd>("NormalDirection"))
{
    contact_diffusive_flux_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_n_.push_back(this->contact_particles_[k]->template getVariableDataByName<Vecd>("NormalDirection"));

        for (auto &diffusion : this->diffusions_)
        {
            std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
            contact_diffusive_flux_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(diffusion_species_name + "Flux"));
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j,
                                  Real surface_area_ij_Neumann,
                                  const StdVec<Real *> &diffusive_flux_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_dt_[m][particle_i] += surface_area_ij_Neumann * diffusive_flux_k[m][particle_j];
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Neumann<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<Real *> &diffusive_flux_k = contact_diffusive_flux_[k];
        Vecd *n_k = contact_n_[k];
        Real *Vol_k = this->contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Vecd n_ij = n_[index_i] - n_k[index_j];
            Real area_ij_Neumann = grad_ijV_j.dot(n_ij);
            getDiffusionChangeRateNeumann(index_i, index_j, area_ij_Neumann, diffusive_flux_k);
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
template <typename... Args>
DiffusionRelaxation<Robin<ContactKernelGradientType>, DiffusionType>::
    DiffusionRelaxation(Args &&...args)
    : DiffusionRelaxation<Contact<ContactKernelGradientType>, DiffusionType>(
          std::forward<Args>(args)...),
      n_(this->particles_->template getVariableDataByName<Vecd>("NormalDirection"))
{
    contact_convection_.resize(this->contact_particles_.size());
    contact_species_infinity_.resize(this->contact_particles_.size());

    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles_k = this->contact_particles_[k];
        contact_n_.push_back(contact_particles_k->template getVariableDataByName<Vecd>("NormalDirection"));

        for (auto &diffusion : this->diffusions_)
        {
            std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
            contact_convection_[k].push_back(
                contact_particles_k->template registerStateVariable<Real>(diffusion_species_name + "Convection"));
            contact_species_infinity_[k].push_back(
                contact_particles_k->template registerSingularVariable<Real>(diffusion_species_name + "Infinity")->Data());
        }
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Robin<ContactKernelGradientType>, DiffusionType>::
    getTransferRateRobin(size_t particle_i, size_t particle_j,
                         Real surface_area_ij_Robin,
                         StdVec<Real *> &transfer_k,
                         StdVec<Real *> &convection_k,
                         StdVec<Real *> &species_infinity_k)
{
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        Real phi_ij = *species_infinity_k[m] - this->diffusion_species_[m][particle_i];
        Real transfer_ij = convection_k[m][particle_j] * phi_ij * surface_area_ij_Robin;
        transfer_k[m][particle_i] += transfer_ij;
    }
}
//=================================================================================================//
template <class ContactKernelGradientType, class DiffusionType>
void DiffusionRelaxation<Robin<ContactKernelGradientType>, DiffusionType>::
    interaction(size_t index_i, Real dt)
{
    this->resetContactTransfer(index_i);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Vecd *n_k = contact_n_[k];
        Real *Vol_k = this->contact_Vol_[k];
        StdVec<Real *> &transfer_k = this->contact_transfer_[k];
        StdVec<Real *> &convection_k = contact_convection_[k];
        StdVec<Real *> &species_infinity_k = contact_species_infinity_[k];

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j, e_ij);
            Vecd n_ij = n_[index_i] - n_k[index_j];
            Real area_ij_Robin = grad_ijV_j.dot(n_ij);
            getTransferRateRobin(index_i, index_j, area_ij_Robin, transfer_k, convection_k, species_infinity_k);
        }
    }
    this->accumulateDiffusionRate(index_i);
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
RungeKuttaStep<DiffusionRelaxationType>::RungeKuttaStep(Args &&...args)
    : DiffusionRelaxationType(std::forward<Args>(args)...)
{
    for (auto &diffusion : this->diffusions_)
    {
        std::string diffusion_species_name = diffusion->DiffusionSpeciesName();
        diffusion_species_s_.push_back(
            this->particles_->template registerStateVariable<Real>(diffusion_species_name + "Intermediate"));
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
FirstStageRK2<DiffusionRelaxationType>::FirstStageRK2(Args &&...args)
    : RungeKuttaStep<DiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void FirstStageRK2<DiffusionRelaxationType>::initialization(size_t index_i, Real dt)
{
    DiffusionRelaxationType::initialization(index_i, dt);

    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_species_s_[m][index_i] = this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename... Args>
SecondStageRK2<DiffusionRelaxationType>::SecondStageRK2(Args &&...args)
    : RungeKuttaStep<DiffusionRelaxationType>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void SecondStageRK2<DiffusionRelaxationType>::update(size_t index_i, Real dt)
{
    DiffusionRelaxationType::update(index_i, dt);
    for (size_t m = 0; m < this->diffusions_.size(); ++m)
    {
        this->diffusion_species_[m][index_i] = 0.5 * this->diffusion_species_s_[m][index_i] +
                                               0.5 * this->diffusion_species_[m][index_i];
    }
}
//=================================================================================================//
template <class DiffusionRelaxationType>
template <typename FirstArg, typename... OtherArgs>
DiffusionRelaxationRK2<DiffusionRelaxationType>::
    DiffusionRelaxationRK2(FirstArg &first_arg, OtherArgs &&...other_args)
    : BaseDynamics<void>(),
      rk2_1st_stage_(first_arg, std::forward<OtherArgs>(other_args)...),
      rk2_2nd_stage_(first_arg, std::forward<OtherArgs>(other_args)...) {}
//=================================================================================================//
template <class DiffusionRelaxationType>
void DiffusionRelaxationRK2<DiffusionRelaxationType>::exec(Real dt)
{
    rk2_1st_stage_.exec(dt);
    rk2_2nd_stage_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_HPP