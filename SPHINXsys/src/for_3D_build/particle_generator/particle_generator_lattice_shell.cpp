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
	{}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		total_volume_plate_ = 10000;
		total_volume_tube_ = 149225.65104552;
		global_avg_thickness_ = 4.0; // 4 plate, 5 tube

		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		particle_spacing_ = 5;
		//lattice_spacing_ = 1;

		avg_particle_volume_ = global_avg_thickness_ * particle_spacing_ * particle_spacing_;
		number_of_particles_ = total_volume_plate_ / avg_particle_volume_;

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
