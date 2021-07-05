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
		// count the number of cells
		int number_of_cells = 0;
		std::unique_ptr<Mesh> mesh(new Mesh(domain_bounds_, lattice_spacing_, 0));
		Vecu number_of_lattices = mesh->NumberOfCells();
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) 
				for (size_t k = 0; k < number_of_lattices[2]; ++k) {
					Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j, k));
					if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
					{
						if (body_shape_->checkContain(particle_position))
						{
							number_of_cells++;
						}
					}
				}
		std::cout << "number of cells: " << number_of_cells << std::endl;

		// calculate the number of particles needed
		size_t total_real_particles = 0;
		int count = 0;
		int interval = int(number_of_cells / number_of_particles_ + 0.5);
		if (interval <= 0) interval = 1;
		int random_int = generateRandomIntInRange(0, interval - 1);

		// add a particle in an idex interval, randomly for each interval
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) 
				for (size_t k = 0; k < number_of_lattices[2]; ++k) {
					Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j, k));
					//std::cout << body_shape_->checkNotFar(particle_position, 0.0);
					if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
					{
						if (body_shape_->checkContain(particle_position))
						{
							// if we hit the random number in the interval, add a particle
							if ( count % interval == random_int )
							{
								createABaseParticle(base_particles, particle_position, avg_particle_volume_, total_real_particles);
							}
							count++;
							// new random number in new interval
							if ( count % interval == 0 )
							{
								count = 0; // reset the count if we hit a new interval
								random_int = generateRandomIntInRange(0, interval - 1);
							}
						}
					}
				}
		base_particles->total_real_particles_ = total_real_particles;
	}
	//=================================================================================================//
}
