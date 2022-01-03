/**
 * @file 	fluid_dynamics_multi_phase.hpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "fluid_dynamics_multi_phase.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template<class RelaxationInnerType>
        RelaxationMultiPhase<RelaxationInnerType>::
            RelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation) :
				RelaxationInnerType(inner_relation), MultiPhaseContactData(contact_relation)
		{
			if (inner_relation.sph_body_ != contact_relation.sph_body_)
			{
				std::cout << "\n Error: the two body_realtions do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_p_.push_back(&(contact_particles_[k]->p_));
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));
			}
		}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>::
            BasePressureRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation) :
				RelaxationMultiPhase<PressureRelaxationInnerType>(inner_relation,
					contact_relation)
		{
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				riemann_solvers_.push_back(CurrentRiemannSolver(*this->material_, *this->contact_material_[k]));
			}
		}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>::
            BasePressureRelaxationMultiPhase(ComplexBodyRelation &complex_relation) :
				BasePressureRelaxationMultiPhase(complex_relation.inner_relation_, complex_relation.contact_relation_) {}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        void BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxationInnerType::Interaction(index_i, dt);

			FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i]);
			Vecd acceleration(0.0);
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->contact_Vol_[k]);
				StdLargeVec<Real>& rho_k = *(this->contact_rho_n_[k]);
				StdLargeVec<Real>& p_k = *(this->contact_p_[k]);
				StdLargeVec<Vecd>& vel_k = *(this->contact_vel_n_[k]);
				CurrentRiemannSolver& riemann_solver_k = riemann_solvers_[k];
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];

					FluidState state_j(rho_k[index_j], vel_k[index_j], p_k[index_j]);
					Real p_star = riemann_solver_k.getPStar(state_i, state_j, e_ij);
					acceleration -= 2.0 * p_star * e_ij * Vol_k[index_j] * dW_ij / state_i.rho_;
				}
			}
			this->dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        Vecd BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>::
            computeNonConservativeAcceleration(size_t index_i)
		{
			Vecd acceleration = PressureRelaxationInnerType::computeNonConservativeAcceleration(index_i);

			Real rho_i = this->rho_n_[index_i];
			Real p_i = this->p_[index_i];
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& rho_k = *(this->contact_rho_n_[k]);
				StdLargeVec<Real>& p_k = *(this->contact_p_[k]);
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];

					Real rho_j = rho_k[index_j];
					Real p_j = p_k[index_j];

					Real p_star = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
					acceleration += (p_i - p_star) * this->Vol_[index_j] * dW_ij * e_ij / rho_i;
				}
			}
			return acceleration;
		}
		//=================================================================================================//
		template<class DensityRelaxationInnerType>
        BaseDensityRelaxationMultiPhase<DensityRelaxationInnerType>::
            BaseDensityRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation) :
				RelaxationMultiPhase<DensityRelaxationInnerType>(inner_relation,
					contact_relation)
			{
				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					riemann_solvers_.push_back(CurrentRiemannSolver(*this->material_, *this->contact_material_[k]));
				}
			}
 		//=================================================================================================//
		template<class DensityRelaxationInnerType>
        BaseDensityRelaxationMultiPhase<DensityRelaxationInnerType>::
            BaseDensityRelaxationMultiPhase(ComplexBodyRelation &complex_relation) :
				BaseDensityRelaxationMultiPhase(complex_relation.inner_relation_, complex_relation.contact_relation_) {}
 		//=================================================================================================//
		template<class DensityRelaxationInnerType>
        void BaseDensityRelaxationMultiPhase<DensityRelaxationInnerType>::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxationInnerType::Interaction(index_i, dt);

			FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i]);
			Real density_change_rate = 0.0;
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->contact_Vol_[k]);
				StdLargeVec<Real>& rho_k = *(this->contact_rho_n_[k]);
				StdLargeVec<Real>& p_k = *(this->contact_p_[k]);
				StdLargeVec<Vecd>& vel_k = *(this->contact_vel_n_[k]);
				CurrentRiemannSolver& riemann_solver_k = riemann_solvers_[k];
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];

					FluidState state_j(rho_k[index_j], vel_k[index_j], p_k[index_j]);
					Vecd vel_star = riemann_solver_k.getVStar(state_i, state_j, e_ij);
					density_change_rate += 2.0 * state_i.rho_ * Vol_k[index_j] * dot(state_i.vel_ - vel_star, e_ij) * dW_ij;
				}
			}
			this->drho_dt_[index_i] += density_change_rate;
		}
       //=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//