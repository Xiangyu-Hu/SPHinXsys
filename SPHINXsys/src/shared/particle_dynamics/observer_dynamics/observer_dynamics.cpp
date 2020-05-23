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
		void CorrectKenelWeightsForInterpolation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Vecd weight_correction(0.0);
			Matd local_configuration(0.0);
			/** Compute the first order consistent kernel weights */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j = (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					Vecd r_ji = -neighboring_particle->r_ij_ * neighboring_particle->e_ij_;
					weight_correction += r_ji * neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
					Vecd gradw_ij = neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
					local_configuration += base_particle_data_j.Vol_ * SimTK::outer(r_ji, gradw_ij);
				}
			}

			/** correction matrix for interacting configuration */
			Matd B_ = GeneralizedInverse(local_configuration);

			/** Add the kernel weight correction to W_ij_ of neighboring particles. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j = (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					Vecd normalized_weight_correction = B_ * weight_correction;
					neighboring_particle->W_ij_ 
						-= dot(normalized_weight_correction, neighboring_particle->e_ij_) * neighboring_particle->dW_ij_;
				}
			}
		}
		//=================================================================================================//
	}
//=================================================================================================//
}

