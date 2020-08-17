/**
 * @file 	observer_dynamics.cpp
 * @brief 	Here, Functions defined in observer_dyanmcis.h are detailed.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
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
		void CorrectInterpolationKernelWeights::ContactInteraction(size_t index_particle_i, Real dt)
		{
			Vecd weight_correction(0.0);
			/** a small number added to diagnal to avoid divide zero */
			Matd local_configuration(Eps);
			/** Compute the first order consistent kernel weights */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_particle_i];
				KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
				CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					CommonRelation& common_relation = contact_common_relations[n];
					size_t index_particle_j = common_relation.j_;
					BaseParticleData& base_particle_data_j = contact_particles_[k]->base_particle_data_[index_particle_j];

					Vecd r_ji = -common_relation.r_ij_ * common_relation.e_ij_;
					weight_correction += r_ji * kernel_value_list[n] * base_particle_data_j.Vol_;
					Vecd gradw_ij = common_relation.dW_ij_ * common_relation.e_ij_;
					local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ji, gradw_ij);
				}
			}

			/** correction matrix for interacting configuration */
			Matd B_ = SimTK::inverse(local_configuration);

			/** Add the kernel weight correction to W_ij_ of neighboring particles. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_particle_i];
				KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
				CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					CommonRelation& common_relation = contact_common_relations[n];
					size_t index_particle_j = common_relation.j_;
					BaseParticleData& base_particle_data_j = contact_particles_[k]->base_particle_data_[index_particle_j];

					Vecd normalized_weight_correction = B_ * weight_correction;
					kernel_value_list[n] -= dot(normalized_weight_correction, common_relation.e_ij_) * common_relation.dW_ij_;
				}
			}
		}
		//=================================================================================================//
	}
//=================================================================================================//
}

