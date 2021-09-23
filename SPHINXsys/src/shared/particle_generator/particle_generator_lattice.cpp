/**
 * @file 	particle_generator_lattice.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */
#include "particle_generator_lattice.h"

#include "geometry.h"
#include "base_body.h"
#include "base_particles.h"
#include "particle_adaptation.h"

namespace SPH {
	//=================================================================================================//
	ParticleGeneratorLattice::ParticleGeneratorLattice()
		: ParticleGenerator(), lattice_spacing_(0), 
		domain_bounds_(0, 0), body_shape_(nullptr)
	{
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->particle_adaptation_->ReferenceSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = sph_body_->body_shape_;
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::
		createABaseParticle(BaseParticles* base_particles, Vecd& particle_position, Real particle_volume)
	{
		base_particles->initializeABaseParticle(particle_position, particle_volume);
	}
	//=================================================================================================//
	ParticleGeneratorMultiResolution::ParticleGeneratorMultiResolution()
		: ParticleGeneratorLattice(), particle_adapation_(nullptr) {}
	//=================================================================================================//
	void ParticleGeneratorMultiResolution::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		particle_adapation_ =
			dynamic_cast<ParticleSpacingByBodyShape*>(sph_body->particle_adaptation_);
		lattice_spacing_ = particle_adapation_->MinimumSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = sph_body_->body_shape_;
	}
	//=================================================================================================//
	void ParticleGeneratorMultiResolution::
		createABaseParticle(BaseParticles* base_particles, Vecd& particle_position, Real particle_volume)
	{
		Real local_particle_spacing 
			= particle_adapation_->getLocalSpacing(*body_shape_, particle_position);
		Real local_particle_volume_ratio = powerN(lattice_spacing_ / local_particle_spacing, Dimensions);
		if ((double)rand() / (RAND_MAX) < local_particle_volume_ratio)
		{
			base_particles->initializeABaseParticle(particle_position, particle_volume / local_particle_volume_ratio);
		}
	}
	//=================================================================================================//
}
