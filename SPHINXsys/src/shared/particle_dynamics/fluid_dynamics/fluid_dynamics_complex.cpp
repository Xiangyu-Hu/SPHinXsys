#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		void DensitySummationComplex::interaction(size_t index_i, Real dt)
		{
			BaseDensitySummationComplex<DensitySummationInner>::interaction(index_i, dt);
			Real sigma = BaseDensitySummationComplex<DensitySummationInner>::ContactSummation(index_i);
			rho_sum_[index_i] += sigma * rho0_ * rho0_ * inv_sigma0_ / mass_[index_i];
		}
		//=================================================================================================//
		void DensitySummationComplexAdaptive::interaction(size_t index_i, Real dt)
		{
			BaseDensitySummationComplex<DensitySummationInnerAdaptive>::interaction(index_i, dt);
			Real sigma = BaseDensitySummationComplex<DensitySummationInnerAdaptive>::ContactSummation(index_i);
			rho_sum_[index_i] += sigma * rho0_ * rho0_ / mass_[index_i] /
								 sph_adaptation_.ReferenceNumberDensity(h_ratio_[index_i]);
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplex::interaction(size_t index_i, Real dt)
		{
			TransportVelocityCorrectionInner::interaction(index_i, dt);

			Vecd acceleration_trans = Vecd::Zero();
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

					// acceleration for transport velocity
					acceleration_trans -= 2.0 * nablaW_ijV_j;
				}
			}

			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
				pos_[index_i] += coefficient_ * smoothing_length_sqr_ * acceleration_trans;
		}
		//=================================================================================================//
		void TransportVelocityCorrectionComplexAdaptive::interaction(size_t index_i, Real dt)
		{
			TransportVelocityCorrectionInnerAdaptive::interaction(index_i, dt);

			Vecd acceleration_trans = Vecd::Zero();
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					Vecd nablaW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

					// acceleration for transport velocity
					acceleration_trans -= 2.0 * nablaW_ijV_j;
				}
			}

			/** correcting particle position */
			if (surface_indicator_[index_i] == 0)
			{
				Real inv_h_ratio = 1.0 / sph_adaptation_.SmoothingLengthRatio(index_i);
				pos_[index_i] += coefficient_ * smoothing_length_sqr_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
			}
		}
		//=================================================================================================//
		void Oldroyd_BIntegration1stHalfWithWall::interaction(size_t index_i, Real dt)
		{
			BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>::interaction(index_i, dt);

			Real rho_i = rho_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration = Vecd::Zero();
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
					/** stress boundary condition. */
					acceleration += 2.0 * tau_i * nablaW_ijV_j / rho_i;
				}
			}

			acc_[index_i] += acceleration;
		}
		//=================================================================================================//
		void Oldroyd_BIntegration2ndHalfWithWall::interaction(size_t index_i, Real dt)
		{
			BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>::interaction(index_i, dt);

			Vecd vel_i = vel_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate = Matd::Zero();
			for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
				Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];

					Matd velocity_gradient = -2.0 * (vel_i - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
					stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
								   (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
				}
			}
			dtau_dt_[index_i] += stress_rate;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//