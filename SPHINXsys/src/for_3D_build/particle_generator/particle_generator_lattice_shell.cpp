/**
 * @file 	particle_generator_lattice.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */
#include "particle_generator_lattice_shell.h"

#include "geometry.h"
#include "base_body.h"
#include "base_particles.h"
#include "particle_adaptation.h"

namespace SPH {
	//=================================================================================================//
	ShellParticleGeneratorLattice::ShellParticleGeneratorLattice()
		: ParticleGenerator(), lattice_spacing_(0), 
		domain_bounds_(0, 0), body_shape_(NULL)
	{
		// total_volume_ = 1000;
		// global_avg_thickness_ = 4;
		// particle_spacing_ = 4;
		// avg_particle_size_ = 4 * 10 * 10; // calculated as: particle spacing^2 * global avg thickness
		number_of_particles_ = 25;
		number_of_cells_ = 10000;
	}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = sph_body_->body_shape_;
	}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::
		createABaseParticle(BaseParticles* base_particles,
			Vecd& particle_position, Real particle_volume, size_t& total_real_particles)
	{
		base_particles->initializeABaseParticle(particle_position, particle_volume);
		total_real_particles++;
	}
}
