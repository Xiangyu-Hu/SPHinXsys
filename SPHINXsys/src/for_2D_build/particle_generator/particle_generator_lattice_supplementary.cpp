/**
 * @file 	particle_generator_lattic_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "particle_generator_lattice.h"

#include "geometry.h"
#include "base_mesh.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	void ParticleGeneratorLattice::createBaseParticles(BaseParticles* base_particles)
	{
		std::unique_ptr<BaseMesh> mesh(new BaseMesh(domain_bounds_, lattice_spacing_, 0));
		Real particle_volume = lattice_spacing_ * lattice_spacing_;
		Vecu number_of_lattices = mesh->NumberOfCellsFromNumberOfGridPoints(mesh->NumberOfGridPoints());
			for (size_t i = 0; i < number_of_lattices[0]; ++i)
				for (size_t j = 0; j < number_of_lattices[1]; ++j)
				{
				Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i,j));
				if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
				{
					if (body_shape_->checkContain(particle_position))
					{
						createABaseParticle(base_particles, particle_position, particle_volume);
					}
				}
			}
	}
	//=================================================================================================//
}
