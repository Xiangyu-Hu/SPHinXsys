//common functions used by 2d buildings only

#include "particle_generator_lattice.h"
#include "base_body.h"
#include "base_particles.h"
#include "array_allocation.h"

namespace SPH {
	//===============================================================//
	void ParticleGeneratorLattice
		::CreateParticles(SPHBody &sph_body)
	{
		size_t number_of_particles = 0;
		Real vol = lattice_spacing_ * lattice_spacing_;
		for (int j = 0; j < number_of_lattices_[1]; ++j)
			for (int i = 0; i < number_of_lattices_[0]; ++i) 
			{
				Point particle_location(lower_bound_[0] + (Real(i) + 0.5)*lattice_spacing_,
					lower_bound_[1] + (Real(j) + 0.5)*lattice_spacing_);
				if(sph_body.BodyContain(particle_location))
				{
						sph_body.GenerateAParticle(particle_location, vol);
					number_of_particles++;
				}
			}

		sph_body.number_of_real_particles_ = number_of_particles;
	}
	//===============================================================//
}