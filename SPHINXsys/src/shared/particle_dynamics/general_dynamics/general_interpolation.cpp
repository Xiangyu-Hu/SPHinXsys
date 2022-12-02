/**
 * @file 	general_interpolation.cpp
 * @brief
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "general_interpolation.h"

namespace SPH
{
	//=================================================================================================//
	CorrectInterpolationKernelWeights::
		CorrectInterpolationKernelWeights(BaseContactRelation &contact_relation) : LocalDynamics(contact_relation.sph_body_),
																				   InterpolationContactData(contact_relation)
	{
		for (size_t k = 0; k != contact_particles_.size(); ++k)
		{
			contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
		}
	}
	//=================================================================================================//
	void CorrectInterpolationKernelWeights::interaction(size_t index_i, Real dt)
	{
		Vecd weight_correction = Vecd::Zero();
		Matd local_configuration = Eps * Matd::Identity();

		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
			Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
				Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
				Vecd gradW_ijV_j = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

				weight_correction += weight_j * r_ji;
				local_configuration += r_ji * gradW_ijV_j.transpose();
			}
		}

		// correction matrix for interacting configuration
		Matd B_ = local_configuration.inverse();
		// Add the kernel weight correction to W_ij_ of neighboring particles.
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				Vecd normalized_weight_correction = B_ * weight_correction;
				contact_neighborhood.W_ij_[n] -= contact_neighborhood.dW_ijV_j_[n] *
												 normalized_weight_correction.dot(contact_neighborhood.e_ij_[n]);
			}
		}
	}
	//=================================================================================================//
}
