#include "fluid_dynamics.h"

using namespace std;

namespace SPH
{
	namespace fluid_dynamics
	{
		//===========================================================//
		void ComputingVorticityInFluidField::InnerInteraction(size_t index_particle_i, Real dt)
		{
			WeaklyCompressibleFluidParticleData &fluid_data_i = particles_->fluid_data_[index_particle_i];
			BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Vecd vort_temp(0);
			Vecd vort(0);
			StdVec<NeighboringParticle>  &neighors = (*current_inner_configuration_)[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				WeaklyCompressibleFluidParticleData &fluid_data_j = particles_->fluid_data_[index_particle_j];


				//low dissipation Riemann problem
				Vecd vel_diff = base_particle_data_j.vel_n_ - base_particle_data_i.vel_n_;
				vort[0] = vel_diff[1] * neighboring_particle.r_ij_[2] - vel_diff[2] * neighboring_particle.r_ij_[1];
				vort[1] = vel_diff[2] * neighboring_particle.r_ij_[0] - vel_diff[0] * neighboring_particle.r_ij_[2];
				vort[2] = vel_diff[0] * neighboring_particle.r_ij_[1] - vel_diff[1] * neighboring_particle.r_ij_[0];
				vort_temp += vort * base_particle_data_j.Vol_ * neighboring_particle.dW_ij_;
			}

			fluid_data_i.vorticity_ = vort_temp;
		}
		//===========================================================//
	}
}