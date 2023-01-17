#include "shell_fluid_interaction.h"

namespace SPH
{
	//=====================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		FluidViscousForceOnShell::FluidViscousForceOnShell(BaseContactRelation &contact_relation)
			: LocalDynamics(contact_relation.sph_body_)
            , FluidShellContactData(contact_relation)
            , Vol_(particles_->Vol_)
            , vel_ave_(*particles_->AverageVelocity())
		{
			particles_->registerVariable(viscous_force_from_fluid_, "ViscousForceFromFluid");

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_fluids_.push_back(&contact_particles_[k]->fluid_);
				mu_.push_back(contact_fluids_[k]->ReferenceViscosity());
				smoothing_length_.push_back(contact_bodies_[k]->sph_adaptation_->ReferenceSmoothingLength());

				contact_rho_n_.push_back(&(contact_particles_[k]->rho_));
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
			}
		}
		//=================================================================================================//
		void FluidViscousForceOnShell::interaction(size_t index_i, Real dt)
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

					force += 2.0 * mu_k * vel_derivative * contact_neighborhood.dW_ijV_j_[n] * Vol_i * particles_->DegeneratedSpacing(index_i);
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
		void FluidAngularConservativeViscousForceOnShell::interaction(size_t index_i, Real dt)
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
					force += eta_ij * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n] * Vol_i * particles_->DegeneratedSpacing(index_i);
				}
			}

			viscous_force_from_fluid_[index_i] = force;
		}
		//=================================================================================================//
	}
}
