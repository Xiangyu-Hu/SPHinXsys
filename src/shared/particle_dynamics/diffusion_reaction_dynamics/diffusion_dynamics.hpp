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
template <class ParticlesType>
GetDiffusionTimeStepSize<ParticlesType>::
    GetDiffusionTimeStepSize(SPHBody &sph_body)
    : BaseDynamics<Real>(sph_body),
      DiffusionReactionSimpleData<ParticlesType>(sph_body)
{
    Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
    diff_time_step_ = this->particles_->diffusion_reaction_material_
                          .getDiffusionTimeStepSize(smoothing_length);
}
//=================================================================================================//
template <class ParticlesType>
BaseDiffusionRelaxation<ParticlesType>::
    BaseDiffusionRelaxation(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      DiffusionReactionSimpleData<ParticlesType>(sph_body),
      material_(this->particles_->diffusion_reaction_material_),
      all_diffusions_(material_.AllDiffusions()),
      diffusion_species_(this->particles_->DiffusionSpecies()),
      gradient_species_(this->particles_->GradientSpecies())
{
    diffusion_dt_.resize(all_diffusions_.size());
    StdVec<std::string> &all_species_names = this->particles_->AllSpeciesNames();
    IndexVector &diffusion_species_indexes = material_.DiffusionSpeciesIndexes();
    for (size_t i = 0; i != all_diffusions_.size(); ++i)
    {
        // Register specie change rate as shared variable
        std::string &diffusion_species_name = all_species_names[diffusion_species_indexes[i]];
        diffusion_dt_[i] = this->particles_->template registerSharedVariable<Real>(diffusion_species_name + "ChangeRate");
    }
} //=================================================================================================//
template <class ParticlesType, class KernelGradientType>
DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    DiffusionRelaxationInner(BaseInnerRelation &inner_relation)
    : BaseDiffusionRelaxation<ParticlesType>(inner_relation.getSPHBody()),
      DataDelegateInner<ParticlesType, DataDelegateEmptyBase>(inner_relation),
      kernel_gradient_(this->particles_) {}
//=================================================================================================//
template <class ParticlesType, class KernelGradientType>
void DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    initializeDiffusionChangeRate(size_t particle_i)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        (*this->diffusion_dt_[m])[particle_i] = 0;
    }
}
//=================================================================================================//
template <class ParticlesType, class KernelGradientType>
void DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        Real diff_coff_ij =
            this->all_diffusions_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
        StdLargeVec<Real> &gradient_species = *this->gradient_species_[m];
        Real phi_ij = gradient_species[particle_i] - gradient_species[particle_j];
        (*this->diffusion_dt_[m])[particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ParticlesType, class KernelGradientType>
void DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    updateSpeciesDiffusion(size_t particle_i, Real dt)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        (*this->diffusion_species_[m])[particle_i] += dt * (*this->diffusion_dt_[m])[particle_i];
    }
}
//=================================================================================================//
template <class ParticlesType, class KernelGradientType>
void DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    interaction(size_t index_i, Real dt)
{
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];

    initializeDiffusionChangeRate(index_i);
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
        Real r_ij_ = inner_neighborhood.r_ij_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        const Vecd &grad_ijV_j = this->kernel_gradient_(index_i, index_j, dW_ijV_j_, e_ij);
        Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
        getDiffusionChangeRate(index_i, index_j, e_ij, area_ij);
    }
}
//=================================================================================================//
template <class ParticlesType, class KernelGradientType>
void DiffusionRelaxationInner<ParticlesType, KernelGradientType>::
    update(size_t index_i, Real dt)
{
    updateSpeciesDiffusion(index_i, dt);
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>::
    BaseDiffusionRelaxationContact(BaseContactRelation &contact_relation)
    : BaseDiffusionRelaxation<ParticlesType>(contact_relation.getSPHBody()),
      DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>(contact_relation)
{
    StdVec<std::string> &all_species_names = this->particles_->AllSpeciesNames();
    contact_gradient_species_names_.resize(this->contact_particles_.size());

    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        size_t l = this->all_diffusions_[m]->gradient_species_index_;
        std::string &inner_species_name_m = all_species_names[l];
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            auto all_species_map_k = this->contact_particles_[k]->AllSpeciesIndexMap();
            if (all_species_map_k.find(inner_species_name_m) != all_species_map_k.end())
            {
                contact_gradient_species_names_[k].push_back(inner_species_name_m);
            }
            else
            {
                std::cout << "\n Error: inner species '" << inner_species_name_m
                          << "' is not found in contact particles" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
        }
    }

    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_kernel_gradients_.push_back(KernelGradientType(this->particles_, this->contact_particles_[k]));
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
DiffusionRelaxationDirichlet<ParticlesType, ContactParticlesType, KernelGradientType>::
    DiffusionRelaxationDirichlet(BaseContactRelation &contact_relation)
    : BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>(contact_relation)
{
    contact_gradient_species_.resize(this->contact_particles_.size());

    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            std::string contact_species_names_k_m = this->contact_gradient_species_names_[k][m];
            auto all_species_map_k = this->contact_particles_[k]->AllSpeciesIndexMap();
            size_t contact_species_index_k_m = all_species_map_k[contact_species_names_k_m];
            StdVec<StdLargeVec<Real>> &all_contact_species_k = this->contact_particles_[k]->all_species_;
            contact_gradient_species_[k].push_back(&all_contact_species_k[contact_species_index_k_m]);
        }
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationDirichlet<ParticlesType, ContactParticlesType, KernelGradientType>::
    getDiffusionChangeRateDirichletContact(size_t particle_i, size_t particle_j, Vecd &e_ij,
                                           Real surface_area_ij, const StdVec<StdLargeVec<Real> *> &gradient_species_k)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        Real diff_coff_ij =
            this->all_diffusions_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
        Real phi_ij = (*this->gradient_species_[m])[particle_i] - (*gradient_species_k[m])[particle_j];
        (*this->diffusion_dt_[m])[particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationDirichlet<ParticlesType, ContactParticlesType, KernelGradientType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdVec<StdLargeVec<Real> *> &gradient_species_k = this->contact_gradient_species_[k];

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij_ = contact_neighborhood.r_ij_[n];
            Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j_, e_ij);
            Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
            getDiffusionChangeRateDirichletContact(index_i, index_j, e_ij, area_ij, gradient_species_k);
        }
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
DiffusionRelaxationNeumann<ParticlesType, ContactParticlesType, KernelGradientType>::
    DiffusionRelaxationNeumann(BaseContactRelation &contact_relation)
    : BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>(contact_relation),
      n_(this->particles_->n_)
{
    contact_heat_flux_.resize(this->contact_particles_.size());
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_n_.push_back(&(this->contact_particles_[k]->n_));
            contact_heat_flux_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HeatFlux");
        }
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationNeumann<ParticlesType, ContactParticlesType, KernelGradientType>::
    getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j,
                                         Real surface_area_ij_Neumann, StdLargeVec<Real> &heat_flux_k)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        (*this->diffusion_dt_[m])[particle_i] += surface_area_ij_Neumann * heat_flux_k[particle_j];
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationNeumann<ParticlesType, ContactParticlesType, KernelGradientType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &heat_flux_k = *(contact_heat_flux_[k]);
        StdLargeVec<Vecd> &n_k = *(contact_n_[k]);

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j_, e_ij);
            Vecd n_ij = n_[index_i] - n_k[index_j];
            Real area_ij_Neumann = grad_ijV_j.dot(n_ij);
            getDiffusionChangeRateNeumannContact(index_i, index_j, area_ij_Neumann, heat_flux_k);
        }
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
DiffusionRelaxationRobin<ParticlesType, ContactParticlesType, KernelGradientType>::
    DiffusionRelaxationRobin(BaseContactRelation &contact_relation)
    : BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType>(contact_relation),
      n_(this->particles_->n_)
{
    contact_convection_.resize(this->contact_particles_.size());
    contact_T_infinity_.resize(this->all_diffusions_.size());

    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_n_.push_back(&(this->contact_particles_[k]->n_));
            contact_convection_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("Convection");
            contact_T_infinity_[m] = this->contact_particles_[k]->template registerGlobalVariable<Real>("T_infinity");
        }
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationRobin<ParticlesType, ContactParticlesType, KernelGradientType>::
    getDiffusionChangeRateRobinContact(size_t particle_i, size_t particle_j,
                                       Real surface_area_ij_Robin, StdLargeVec<Real> &convection_k, Real &T_infinity_k)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        Real phi_ij = T_infinity_k - (*this->diffusion_species_[m])[particle_i];
        (*this->diffusion_dt_[m])[particle_i] += convection_k[particle_j] * phi_ij * surface_area_ij_Robin;
    }
}
//=================================================================================================//
template <class ParticlesType, class ContactParticlesType, class KernelGradientType>
void DiffusionRelaxationRobin<ParticlesType, ContactParticlesType, KernelGradientType>::
    interaction(size_t index_i, Real dt)
{
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(contact_n_[k]);

        StdLargeVec<Real> &convection_k = *(contact_convection_[k]);
        Real &T_infinity_k = *(contact_T_infinity_[k]);

        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];

            const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j_, e_ij);
            Vecd n_ij = n_[index_i] - n_k[index_j];
            Real area_ij_Robin = grad_ijV_j.dot(n_ij);
            getDiffusionChangeRateRobinContact(index_i, index_j, area_ij_Robin, convection_k, T_infinity_k);
        }
    }
}
//=================================================================================================//
template <class ParticlesType>
InitializationRK<ParticlesType>::
    InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &diffusion_species_s)
    : BaseDiffusionRelaxation<ParticlesType>(sph_body),
      diffusion_species_s_(diffusion_species_s) {}
//=================================================================================================//
template <class ParticlesType>
void InitializationRK<ParticlesType>::
    update(size_t index_i, Real dt)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        diffusion_species_s_[m][index_i] = (*this->diffusion_species_[m])[index_i];
    }
}
//=================================================================================================//
template <class FirstStageType>
void SecondStageRK2<FirstStageType>::
    updateSpeciesDiffusion(size_t particle_i, Real dt)
{
    for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
    {
        (*this->diffusion_species_[m])[particle_i] =
            0.5 * diffusion_species_s_[m][particle_i] +
            0.5 * ((*this->diffusion_species_[m])[particle_i] + dt * (*this->diffusion_dt_[m])[particle_i]);
    }
}
//=================================================================================================//
template <class FirstStageType>
template <typename... ContactArgsType>
DiffusionRelaxationRK2<FirstStageType>::
    DiffusionRelaxationRK2(typename FirstStageType::BodyRelationType &body_relation, ContactArgsType &&...contact_args)
    : BaseDynamics<void>(body_relation.getSPHBody()),
      rk2_initialization_(body_relation.getSPHBody(), diffusion_species_s_),
      rk2_1st_stage_(body_relation, std::forward<ContactArgsType>(contact_args)...),
      rk2_2nd_stage_(body_relation, diffusion_species_s_, std::forward<ContactArgsType>(contact_args)...),
      all_diffusions_(rk2_1st_stage_.AllDiffusions())
{
    diffusion_species_s_.resize(all_diffusions_.size());
    StdVec<std::string> &all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
    for (size_t i = 0; i != all_diffusions_.size(); ++i)
    {
        // Register diffusion species intermediate
        size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
        std::string &diffusion_species_name = all_species_names[diffusion_species_index];
        rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
    }
}
//=================================================================================================//
template <class FirstStageType>
void DiffusionRelaxationRK2<FirstStageType>::exec(Real dt)
{
    rk2_initialization_.exec();
    rk2_1st_stage_.exec(dt);
    rk2_2nd_stage_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_HPP