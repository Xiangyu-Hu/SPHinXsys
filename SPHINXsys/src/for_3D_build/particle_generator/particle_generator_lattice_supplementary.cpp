//common functions used by 3d buildings only

#include "particle_generator_lattice.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	void ParticleGeneratorLattice::CreateBaseParticles(BaseParticles* base_particles)
	{
		size_t number_of_particles = 0;
		Real particle_volume = lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
		Vecu number_of_lattices = mesh_->NumberOfCells();
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) 
				for (size_t k = 0; k < number_of_lattices[2]; ++k) {
					Point particle_position = mesh_->CellPositionFromIndexes(Vecu(i, j, k));
					if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
					{
						if (body_shape_->checkContain(particle_position))
						{
								base_particles->initializeABaseParticle(particle_position, particle_volume);
							number_of_particles++;
						}
					}
				}
		sph_body_->number_of_particles_ = number_of_particles;
	}
	//=================================================================================================//
}
