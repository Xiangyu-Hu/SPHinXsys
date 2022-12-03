#include "fluid_structure_interaction.h"

namespace SPH
{
	//=====================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		FluidViscousForceOnSolid::FluidViscousForceOnSolid(BaseContactRelation &contact_relation)
			: LocalDynamics(contact_relation.sph_body_), FSIContactData(contact_relation),
			  Vol_(particles_->Vol_), vel_ave_(*particles_->AverageVelocity())
		{
			particles_->registerVariable(viscous_force_from_fluid_, "ViscousForceFromFluid");
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_fluids_.push_back(&contact_particles_[k]->fluid_);
				contact_rho_n_.push_back(&(contact_particles_[k]->rho_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));

				mu_.push_back(contact_fluids_[k]->ReferenceViscosity());
				smoothing_length_.push_back(contact_bodies_[k]->sph_adaptation_->ReferenceSmoothingLength());
			}
		}
		//=================================================================================================//
		void FluidViscousForceOnSolid::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			const Vecd &vel_ave_i = vel_ave_[index_i];

			Vecd force = Vecd::Zero();
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					Vecd vel_derivative = 2.0 * (vel_ave_i - vel_n_k[index_j]) /
										  (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);

					force += 2.0 * mu_k * vel_derivative *
							 Vol_i * contact_neighborhood.dW_ijV_j_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnSolid::interaction(size_t index_i, Real dt)
		{
			Real Vol_i = Vol_[index_i];
			const Vecd &vel_ave_i = vel_ave_[index_i];

			Vecd force = Vecd::Zero();
			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Real mu_k = mu_[k];
				Real smoothing_length_k = smoothing_length_[k];
				StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
				StdLargeVec<Vecd> &vel_n_k = *(contact_vel_n_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.) */
					Real v_r_ij =  contact_neighborhood.r_ij_[n] * (vel_ave_i - vel_n_k[index_j]).dot(contact_neighborhood.e_ij_[n]);
					Real vel_difference = 0.0 * (vel_ave_i - vel_n_k[index_j]).norm() * contact_neighborhood.r_ij_[n];
					Real eta_ij = 8.0 * SMAX(mu_k, rho_n_k[index_j] * vel_difference) * v_r_ij /
								  (contact_neighborhood.r_ij_[n] * contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_k);
					force += eta_ij * Vol_i * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		TotalViscousForceOnSolid ::TotalViscousForceOnSolid(SPHBody &sph_body)
			: LocalDynamicsReduce<Vecd, ReduceSum<Vecd>>(sph_body, Vecd::Zero()),
			  SolidDataSimple(sph_body),
			  viscous_force_from_fluid_(*particles_->getVariableByName<Vecd>("ViscousForceFromFluid"))
		{
			quantity_name_ = "TotalViscousForceOnSolid";
		}
		//=================================================================================================//
		Vecd TotalViscousForceOnSolid::reduce(size_t index_i, Real dt)
		{
			return viscous_force_from_fluid_[index_i];
		}
		//=================================================================================================//
		TotalForceOnSolid::TotalForceOnSolid(SPHBody &sph_body)
			: LocalDynamicsReduce<Vecd, ReduceSum<Vecd>>(sph_body, Vecd::Zero()),
			  SolidDataSimple(sph_body),
			  force_from_fluid_(*particles_->getVariableByName<Vecd>("ForceFromFluid"))
		{
			quantity_name_ = "TotalForceOnSolid";
		}
		//=================================================================================================//
		Vecd TotalForceOnSolid::reduce(size_t index_i, Real dt)
		{
			return force_from_fluid_[index_i];
		}
		//=================================================================================================//
		InitializeDisplacement::
			InitializeDisplacement(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp)
			: LocalDynamics(sph_body), ElasticSolidDataSimple(sph_body),
			  pos_temp_(pos_temp), pos_(particles_->pos_) {}
		//=================================================================================================//
		void InitializeDisplacement::update(size_t index_i, Real dt)
		{
			pos_temp_[index_i] = pos_[index_i];
		}
		//=================================================================================================//
		UpdateAverageVelocityAndAcceleration::
			UpdateAverageVelocityAndAcceleration(SPHBody &sph_body, StdLargeVec<Vecd> &pos_temp)
			: LocalDynamics(sph_body), ElasticSolidDataSimple(sph_body),
			  pos_temp_(pos_temp), pos_(particles_->pos_),
			  vel_ave_(particles_->vel_ave_),
			  acc_ave_(particles_->acc_ave_) {}
		//=================================================================================================//
		void UpdateAverageVelocityAndAcceleration::update(size_t index_i, Real dt)
		{
			Vecd updated_vel_ave = (pos_[index_i] - pos_temp_[index_i]) / (dt + Eps);
			acc_ave_[index_i] = (updated_vel_ave - vel_ave_[index_i]) / (dt + Eps);
			vel_ave_[index_i] = updated_vel_ave;
		}
		//=================================================================================================//
		AverageVelocityAndAcceleration::
			AverageVelocityAndAcceleration(SolidBody &solid_body)
			: initialize_displacement_(solid_body, pos_temp_),
			  update_averages_(solid_body, pos_temp_)
		{
			solid_body.getBaseParticles().registerVariable(pos_temp_, "TemporaryPosition");
		}
		//=================================================================================================//
	}
}
