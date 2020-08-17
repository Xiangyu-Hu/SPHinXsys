/**
 * @file 	fluid_dynamics_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "fluid_dynamics.h"

using namespace std;

namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
//=================================================================================================//
		void ComputingVorticityInFluidField::InnerInteraction(size_t index_particle_i, Real dt)
		{
			FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real vort_temp = 0.0;
			Neighborhood& inner_neighborhood = inner_configuration_[index_particle_i];
			CommonRelationList& inner_common_relations = inner_neighborhood.common_relation_list_;
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				CommonRelation& common_relation = inner_common_relations[n];
				size_t index_particle_j = common_relation.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];

				Vecd vel_diff = base_particle_data_j.vel_n_ - base_particle_data_i.vel_n_;
				Vecd r_ij = common_relation.r_ij_ * common_relation.e_ij_;
				Real vort = vel_diff[0] * r_ij[1] - vel_diff[1] * r_ij[0];
				vort_temp -= vort * base_particle_data_j.Vol_ * common_relation.dW_ij_;
			}

			fluid_data_i.vorticity_ = upgradeToVector3D(vort_temp);
		}
//=================================================================================================//
	}
//=================================================================================================//
}