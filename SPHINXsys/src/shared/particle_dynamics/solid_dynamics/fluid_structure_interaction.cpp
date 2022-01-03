/**
 * @file 	fluid_structure_interaction.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "fluid_structure_interaction.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		FluidViscousForceOnSolid::
			FluidViscousForceOnSolid(BaseBodyRelationContact &contact_relation)
			: InteractionDynamics(*contact_relation.sph_body_),
			  FSIContactData(contact_relation),
			  Vol_(particles_->Vol_), vel_ave_(particles_->vel_ave_)
		{
			particles_->registerAVariable<indexVector, Vecd>(viscous_force_from_fluid_, "ViscousForceFromFluid");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_n_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_n_));

				mu_.push_back(contact_material_[k]->ReferenceViscosity());
				smoothing_length_.push_back(contact_bodies_[k]->sph_adaptation_->ReferenceSmoothingLength());
			}
		}
		//=================================================================================================//
		void FluidViscousForceOnSolid::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			const Vecd &vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					Vecd vel_derivative = 2.0 * (vel_ave_i - vel_n_k[index_j]) /
										  (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);

					force += 2.0 * mu_k * vel_derivative *
							 Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnSolid::Interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			const Vecd &vel_ave_i = vel_ave_[index_i];

			Vecd force(0);
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.) */
					Real v_r_ij = dot(vel_ave_i - vel_n_k[index_j],
									  contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n]);
					Real vel_difference = 0.0 * (vel_ave_i - vel_n_k[index_j]).norm() * contact_neighborhood.r_ij_[n];
					Real eta_ij = 8.0 * SMAX(mu_k, rho_n_k[index_j] * vel_difference) * v_r_ij /
								  (contact_neighborhood.r_ij_[n] * contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);
					force += eta_ij * Vol_i * Vol_k[index_j] * contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		TotalViscousForceOnSolid ::TotalViscousForceOnSolid(SolidBody &solid_body)
			: ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>(solid_body),
			  SolidDataSimple(solid_body),
			  viscous_force_from_fluid_(*particles_->getVariableByName<indexVector, Vecd>("ViscousForceFromFluid"))
		{
			quantity_name_ = "TotalViscousForceOnSolid";
			initial_reference_ = Vecd(0);
		}
		//=================================================================================================//
		Vecd TotalViscousForceOnSolid::ReduceFunction(size_t index_i, Real dt)
		{
			return viscous_force_from_fluid_[index_i];
		}
		//=================================================================================================//
		TotalForceOnSolid::TotalForceOnSolid(SolidBody &solid_body)
			: ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>(solid_body),
			  SolidDataSimple(solid_body),
			  force_from_fluid_(particles_->force_from_fluid_)
		{
			quantity_name_ = "TotalForceOnSolid";
			initial_reference_ = Vecd(0);
		}
		//=================================================================================================//
		Vecd TotalForceOnSolid::ReduceFunction(size_t index_i, Real dt)
		{
			return force_from_fluid_[index_i];
		}
		//=================================================================================================//
		InitializeDisplacement::
			InitializeDisplacement(SolidBody &solid_body, StdLargeVec<Vecd> &pos_temp)
			: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body),
			  pos_temp_(pos_temp), pos_n_(particles_->pos_n_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_)
		{
		}
		//=================================================================================================//
		void InitializeDisplacement::Update(size_t index_i, Real dt)
		{
			pos_temp_[index_i] = pos_n_[index_i];
		}
		//=================================================================================================//
		void UpdateAverageVelocityAndAcceleration::Update(size_t index_i, Real dt)
		{
			Vecd updated_vel_ave = (pos_n_[index_i] - pos_temp_[index_i]) / (dt + Eps);
			dvel_dt_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
			vel_ave_[index_i] = updated_vel_ave;
		}
		//=================================================================================================//
		AverageVelocityAndAcceleration::
			AverageVelocityAndAcceleration(SolidBody &solid_body)
			: initialize_displacement_(solid_body, pos_temp_),
			  update_averages_(solid_body, pos_temp_)
		{
			solid_body.base_particles_->registerAVariable<indexVector, Vecd>(pos_temp_, "TemporaryPosition");
		}
		//=================================================================================================//
	}
}
