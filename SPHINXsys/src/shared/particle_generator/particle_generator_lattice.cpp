/**
 * @file 	particle_generator_lattice.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "particle_generator_lattice.h"
#include "base_mesh.h"
#include "base_body.h"

namespace SPH {
	//=================================================================================================//
	ParticleGeneratorLattice::ParticleGeneratorLattice()
		: ParticleGenerator(), lattice_spacing_(0), 
		lower_bound_(0), upper_bound_(0), body_shape_(NULL)
	{
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		lattice_spacing_ = sph_body_->particle_spacing_;
		sph_body_->getSPHSystemBound(lower_bound_, upper_bound_);
		mesh_ = std::unique_ptr<Mesh>(new Mesh(lower_bound_, upper_bound_, lattice_spacing_));
		body_shape_ = sph_body_->body_shape_;
	}
	//=================================================================================================//
}
