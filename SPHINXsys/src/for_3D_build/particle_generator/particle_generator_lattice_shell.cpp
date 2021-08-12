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
	ShellParticleGeneratorLattice::ShellParticleGeneratorLattice(Real global_avg_thickness)
		: ParticleGenerator(), lattice_spacing_(0), 
		domain_bounds_(0, 0), body_shape_(NULL), total_volume_(0), global_avg_thickness_(global_avg_thickness)
	{}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		number_of_cells_ = 0;
		total_real_particles_ = 0;
		
		sph_body_ = sph_body;
		lattice_spacing_ = global_avg_thickness_ * 0.25;
		particle_spacing_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = sph_body_->body_shape_;

		avg_particle_volume_ = global_avg_thickness_ * particle_spacing_ * particle_spacing_;
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
