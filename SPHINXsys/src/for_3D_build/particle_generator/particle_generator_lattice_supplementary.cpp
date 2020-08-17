//common functions used by 3d buildings only

#include "particle_generator_lattice.h"
#include "base_body.h"
#include "base_particles.h"
#include "array_allocation.h"
#include "base_kernel.h"

namespace SPH {
	//===============================================================//
	void ParticleGeneratorLattice::CreateBaseParticles(BaseParticles* base_particles)
	{
		size_t number_of_particles = 0;
		Real vol = lattice_spacing_ * lattice_spacing_*lattice_spacing_;
		Real sigma = ComputeReferenceNumberDensity();
		for (int i = 0; i < number_of_lattices_[0]; ++i)
			for (int j = 0; j < number_of_lattices_[1]; ++j) 
				for (int k = 0; k < number_of_lattices_[2]; ++k){
				Point particle_location(lower_bound_[0] + (i + 0.5)*lattice_spacing_,
					lower_bound_[1] + (j + 0.5)*lattice_spacing_,
					lower_bound_[2] + (k + 0.5)*lattice_spacing_);
				if (body_shape_.checkContain(particle_location))
				{
					base_particles->InitializeABaseParticle(particle_location, vol, sigma);
					number_of_particles++;
				}
			}

		body_.number_of_particles_ = number_of_particles;
	}
	//===============================================================//
	Real ParticleGeneratorLattice::ComputeReferenceNumberDensity()
	{
		Real sigma(0);
		Real cutoff_radius = body_.kernel_->GetCutOffRadius();
		int search_range = int(cutoff_radius / body_.particle_spacing_) + 1;
		for (int k = -search_range; k <= search_range; ++k)
			for (int j = -search_range; j <= search_range; ++j)
				for (int i = -search_range; i <= search_range; ++i)
				{
					Point particle_location(Real(i) * lattice_spacing_, 
						Real(j) * lattice_spacing_, Real(k) * lattice_spacing_);
					if (particle_location.norm() < cutoff_radius)
						sigma += body_.kernel_->W(particle_location);
				}
		return sigma;
	}
	//===============================================================//
}