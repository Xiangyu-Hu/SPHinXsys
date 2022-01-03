/**
 * @file 	observer_dynamics.cpp
 * @brief 	Here, Functions defined in observer_dyanmcis.h are detailed.
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "observer_dynamics.h"
//=================================================================================================//
using namespace SimTK;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace observer_dynamics
	{
		//=================================================================================================//
		CorrectInterpolationKernelWeights::
			CorrectInterpolationKernelWeights(BaseBodyRelationContact &contact_relation) : 
			InteractionDynamics(*contact_relation.sph_body_),
			InterpolationContactData(contact_relation)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void CorrectInterpolationKernelWeights::Interaction(size_t index_i, Real dt)
		{
			Vecd weight_correction(0.0);
			Matd local_configuration(Eps); // small number added to diagonal to avoid divide zero
			// Compute the first order consistent kernel weights
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
					Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
					Vecd gradw_ij = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];

					weight_correction += r_ji * weight_j;
					local_configuration += Vol_k[index_j] * SimTK::outer(r_ji, gradw_ij);
				}
			}

			// correction matrix for interacting configuration
			Matd B_ = SimTK::inverse(local_configuration);

			// Add the kernel weight correction to W_ij_ of neighboring particles.
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					Vecd normalized_weight_correction = B_ * weight_correction;
					contact_neighborhood.W_ij_[n] 
						-= dot(normalized_weight_correction, contact_neighborhood.e_ij_[n])
						 * contact_neighborhood.dW_ij_[n];
				}
			}
		}
		//=================================================================================================//
	}
//=================================================================================================//
}

