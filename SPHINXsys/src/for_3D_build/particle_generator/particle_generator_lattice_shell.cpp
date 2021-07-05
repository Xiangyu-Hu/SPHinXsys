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
		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->particle_adaptation_->ReferenceSpacing();

		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = sph_body_->body_shape_;
		// BoundingBox bb = body_shape_->findBounds();
		// std::cout << "bb: " << bb.first[0] << bb.first[1] << bb.first[2];

		total_volume_ = 149225.65104552; // plate 10000, tube 149225.65104552
		global_avg_thickness_ = 5.0; // 4 plate, 5 tube
		avg_particle_volume_ = global_avg_thickness_ * 5 * 5;
		Real number_of_particles = total_volume_ / avg_particle_volume_ + 0.5;
		number_of_particles_ = int(number_of_particles);
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
