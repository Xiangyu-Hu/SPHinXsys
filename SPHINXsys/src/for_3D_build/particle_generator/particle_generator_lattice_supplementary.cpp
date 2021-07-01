//common functions used by 3d buildings only

#include "particle_generator_lattice.h"
#include "particle_generator_lattice_shell.h"

#include "geometry.h"
#include "base_mesh.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	void ParticleGeneratorLattice::createBaseParticles(BaseParticles* base_particles)
	{
		std::unique_ptr<Mesh> mesh(new Mesh(domain_bounds_, lattice_spacing_, 0));
		size_t total_real_particles = 0;
		Real particle_volume = lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
		Vecu number_of_lattices = mesh->NumberOfCells();
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) 
				for (size_t k = 0; k < number_of_lattices[2]; ++k) {
					Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j, k));
					if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
					{
						if (body_shape_->checkContain(particle_position))
						{
							createABaseParticle(base_particles, particle_position, particle_volume, total_real_particles);
						}
					}
				}
		base_particles->total_real_particles_ = total_real_particles;
	}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::createBaseParticles(BaseParticles* base_particles)
	{
		int count = 0;
		std::unique_ptr<Mesh> mesh(new Mesh(domain_bounds_, lattice_spacing_, 0));
		size_t total_real_particles = 0;
		Real particle_volume = lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
		Vecu number_of_lattices = mesh->NumberOfCells();
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) 
				for (size_t k = 0; k < number_of_lattices[2]; ++k) {
					Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j, k));
					if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
					{
						if (body_shape_->checkContain(particle_position))
						{
							if (count % 10 == 0)
							{
								createABaseParticle(base_particles, particle_position, particle_volume, total_real_particles);
							}
							count++;
						}
					}
				}
		base_particles->total_real_particles_ = total_real_particles;
	}
	//=================================================================================================//
}
