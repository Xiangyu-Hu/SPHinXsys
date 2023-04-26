#include "general_dynamics_wkgc.h"

namespace SPH
{
	//=================================================================================================//
	void GlobalCorrectionMatrix::interaction(size_t index_i, Real dt)
	{
		Matd local_configuration = Eps * Matd::Identity();

		/** Inner interaction. */

		const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
			local_configuration -= r_ji * gradW_ij.transpose();
		}
		/** Contact interaction. */

		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				Vecd gradW_ij = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
				Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
				local_configuration -= r_ji * gradW_ij.transpose();
			}
		}
		B_[index_i] = (pow(local_configuration.determinant(), beta_) * local_configuration.inverse() + alpha_ * Matd::Identity())
			           / (alpha_ * Matd::Identity().determinant() + pow(local_configuration.determinant(), beta_));
	}
	//=================================================================================================//
}
//=================================================================================================//