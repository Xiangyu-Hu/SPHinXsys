#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_HPP
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_HPP

#include "particle_dynamics_diffusion_reaction_with_boundary.h"

namespace SPH
{
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	DiffusionReactionInitialConditionWithNeumann<DiffusionReactionParticlesType>::
		DiffusionReactionInitialConditionWithNeumann(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<DiffusionReactionParticlesType>(sph_body),
		heat_flux_(*(this->particles_->template getVariableByName<Real>("HeatFlux"))) {}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	RelaxationOfAllDiffusionSpeciesNeumann<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		RelaxationOfAllDiffusionSpeciesNeumann(BaseContactRelation& contact_relation)
		: RelaxationOfAllDiffusionSpeciesContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>(contact_relation),
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
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesNeumann<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j,
			Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k)
	{
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			(*this->diffusion_dt_[m])[particle_i] += surface_area_ij_Neumann * heat_flux_k[particle_j];
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesNeumann<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		interaction(size_t index_i, Real dt)
	{
		DiffusionReactionParticlesType* particles = this->particles_;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& heat_flux_k = *(contact_heat_flux_[k]);
			StdLargeVec<Vecd>& n_k = *(contact_n_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];

				const Vecd& grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Vecd n_ij = n_[index_i] - n_k[index_j];
				Real area_ij_Neumann = grad_ijV_j.dot(n_ij);
				getDiffusionChangeRateNeumannContact(index_i, index_j, area_ij_Neumann, heat_flux_k);
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	DiffusionReactionInitialConditionWithRobin<DiffusionReactionParticlesType>::
		DiffusionReactionInitialConditionWithRobin(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<DiffusionReactionParticlesType>(sph_body)
	{
		this->particles_->registerVariable(convection_, "Convection");
		this->particles_->registerGlobalVariable("T_infinity", &T_infinity_);
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	RelaxationOfAllDiffusionSpeciesRobin<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		RelaxationOfAllDiffusionSpeciesRobin(BaseContactRelation& contact_relation)
		: RelaxationOfAllDiffusionSpeciesContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>(contact_relation),
		n_(this->particles_->n_)
	{
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				contact_n_.push_back(&(this->contact_particles_[k]->n_));
				contact_convection_.push_back(this->contact_particles_[k]->getVariableByName<Real>("Convection"));
				contact_T_infinity_[m] = this->contact_particles_[k]->getGlobalVariableByName<Real>("T_infinity");
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesRobin<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		getDiffusionChangeRateRobinContact(size_t particle_i, size_t particle_j,
			Real surface_area_ij_Robin, StdLargeVec<Real>& convection_k, Real& T_infinity_k)
	{
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			Real phi_ij = T_infinity_k - (*this->diffusion_species_[m])[particle_i];
			(*this->diffusion_dt_[m])[particle_i] += convection_k[particle_j] * phi_ij * surface_area_ij_Robin;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesRobin<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		interaction(size_t index_i, Real dt)
	{
		DiffusionReactionParticlesType* particles = this->particles_;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Vecd>& n_k = *(contact_n_[k]);

			StdLargeVec<Real>& convection_k = *(contact_convection_[k]);
			Real& T_infinity_k = *(contact_T_infinity_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];

				const Vecd& grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Vecd n_ij = n_[index_i] - n_k[index_j];
				Real area_ij_Robin = grad_ijV_j.dot(n_ij);
				getDiffusionChangeRateRobinContact(index_i, index_j, area_ij_Robin, convection_k, T_infinity_k);
			}
		}
	}
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_HPP