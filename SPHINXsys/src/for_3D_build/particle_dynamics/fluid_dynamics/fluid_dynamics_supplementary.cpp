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

			Vecd vort_temp(0);
			Vecd vort(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				FluidParticleData &fluid_data_j = particles_->fluid_particle_data_[index_particle_j];


				//low dissipation Riemann problem
				Vecd vel_diff = base_particle_data_j.vel_n_ - base_particle_data_i.vel_n_;
				Vecd r_ij = neighboring_particle.r_ij_ * neighboring_particle.e_ij_;
				vort[0] = vel_diff[1] * r_ij[2] - vel_diff[2] * r_ij[1];
				vort[1] = vel_diff[2] * r_ij[0] - vel_diff[0] * r_ij[2];
				vort[2] = vel_diff[0] * r_ij[1] - vel_diff[1] * r_ij[0];
				vort_temp += vort * base_particle_data_j.Vol_ * neighboring_particle.dW_ij_;
			}

			fluid_data_i.vorticity_ = vort_temp;
		}
//=================================================================================================//
	}
//=================================================================================================//
}