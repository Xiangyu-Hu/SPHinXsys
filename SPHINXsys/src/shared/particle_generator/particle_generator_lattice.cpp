/**
 * @file 	particle_generator_lattice.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */
#include "particle_generator_lattice.h"

#include "complex_shape.h"
#include "base_body.h"
#include "base_particles.h"
#include "adaptation.h"

namespace SPH
{
	//=================================================================================================//
	ParticleGeneratorLattice::ParticleGeneratorLattice()
		: ParticleGenerator(), lattice_spacing_(0),
		  domain_bounds_(0, 0), body_shape_(nullptr)
	{
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::initialize(SPHBody *sph_body)
	{
		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->sph_adaptation_->ReferenceSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = &sph_body_->body_shape_;
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::
		createABaseParticle(BaseParticles *base_particles, Vecd &particle_position, Real particle_volume)
	{
		base_particles->initializeABaseParticle(particle_position, particle_volume);
	}
	//=================================================================================================//
	ParticleGeneratorMultiResolution::ParticleGeneratorMultiResolution()
		: ParticleGeneratorLattice(), particle_adapation_(nullptr) {}
	//=================================================================================================//
	void ParticleGeneratorMultiResolution::initialize(SPHBody *sph_body)
	{
		sph_body_ = sph_body;
		particle_adapation_ =
			DynamicCast<ParticleSpacingByBodyShape>(this, sph_body->sph_adaptation_);
		lattice_spacing_ = particle_adapation_->MinimumSpacing();
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = &sph_body_->body_shape_;
	}
	//=================================================================================================//
	void ParticleGeneratorMultiResolution::
		createABaseParticle(BaseParticles *base_particles, Vecd &particle_position, Real particle_volume)
	{
		Real local_particle_spacing = particle_adapation_->getLocalSpacing(*body_shape_, particle_position);
		Real local_particle_volume_ratio = powerN(lattice_spacing_ / local_particle_spacing, Dimensions);
		if ((double)rand() / (RAND_MAX) < local_particle_volume_ratio)
		{
			base_particles->initializeABaseParticle(particle_position, particle_volume / local_particle_volume_ratio);
		}
	}
	//=================================================================================================//
	ShellParticleGeneratorLattice::ShellParticleGeneratorLattice(Real global_avg_thickness)
		: ParticleGeneratorLattice(), total_volume_(0), global_avg_thickness_(global_avg_thickness)
	{
	}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		domain_bounds_ = sph_body_->getSPHSystemBounds();
		body_shape_ = &sph_body_->body_shape_;

		number_of_cells_ = 0;
		particle_spacing_ = sph_body_->sph_adaptation_->ReferenceSpacing();
		lattice_spacing_ = 0.5 * global_avg_thickness_;
		avg_particle_volume_ = powerN(particle_spacing_, Dimensions - 1) * global_avg_thickness_;
	}
	//=================================================================================================//
}
