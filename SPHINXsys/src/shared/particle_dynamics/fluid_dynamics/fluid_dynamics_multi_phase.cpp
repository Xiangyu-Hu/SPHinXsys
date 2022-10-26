/**
 * @file 	fluid_dynamics_multi_phase.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "fluid_dynamics_multi_phase.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		ViscousAccelerationMultiPhase::ViscousAccelerationMultiPhase(BaseInnerRelation &inner_relation,
																	 BaseContactRelation &contact_relation)
			: ViscousAccelerationInner(inner_relation), MultiPhaseContactData(contact_relation)
		{
			if (&inner_relation.sph_body_ != &contact_relation.sph_body_)
			{
				std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_fluids_.push_back(&contact_particles_[k]->fluid_);
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
			}
		}
		//=================================================================================================//
		ViscousAccelerationMultiPhase::
			ViscousAccelerationMultiPhase(ComplexRelation &complex_relation)
			: ViscousAccelerationMultiPhase(complex_relation.inner_relation_,
											complex_relation.contact_relation_) {}
		//=================================================================================================//
		void ViscousAccelerationMultiPhase::interaction(size_t index_i, Real dt)
		{
			ViscousAccelerationInner::interaction(index_i, dt);

			Real rho_i = this->rho_[index_i];
			const Vecd &vel_i = this->vel_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				Real mu_j = this->contact_fluids_[k]->ReferenceViscosity();
				StdLargeVec<Vecd> &vel_k = *(this->contact_vel_n_[k]);
				Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0 * (vel_i - vel_k[index_j]) /
									 (r_ij + 0.01 * this->smoothing_length_);
					Real mu_ij = 2.0 * this->mu_ * mu_j / (this->mu_ + mu_j);
					acceleration += 2.0 * mu_ij * vel_derivative *
									contact_neighborhood.dW_ijV_j_[n] / rho_i;
				}
			}

			acc_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		MultiPhaseColorFunctionGradient::
			MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation)
			: LocalDynamics(contact_relation.sph_body_), MultiPhaseData(contact_relation),
			  rho0_(sph_body_.base_material_->ReferenceDensity()), Vol_(particles_->Vol_),
			  pos_div_(*particles_->getVariableByName<Real>("PositionDivergence")),
			  surface_indicator_(particles_->surface_indicator_)
		{
			particles_->registerVariable(color_grad_, "ColorGradient");
			particles_->registerVariable(surface_norm_, "SurfaceNormal");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
				contact_rho0_.push_back(rho0_k);
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void MultiPhaseColorFunctionGradient::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			Vecd gradient(0.0);
			if (surface_indicator_[index_i])
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					Real rho0_k = contact_rho0_[k];
					StdLargeVec<Real> &contact_Vol_k = *(contact_Vol_[k]);
					Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						/** Norm of interface.*/
						Real rho_ij = rho0_ / (rho0_ + rho0_k);
						Real area_ij = (Vol_i * Vol_i + contact_Vol_k[index_j] * contact_Vol_k[index_j]) *
									   contact_neighborhood.dW_ijV_j_[n] / contact_Vol_k[index_j];
						gradient += rho_ij * area_ij * contact_neighborhood.e_ij_[n] / Vol_i;
					}
				}
			}
			color_grad_[index_i] = gradient;
			surface_norm_[index_i] = gradient / (gradient.norm() + TinyReal);
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//