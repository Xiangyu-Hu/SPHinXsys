/**
 * @file 	particle_generator_lattice.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "particle_generator_lattice.h"
#include "mesh_cell_linked_list.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	ParticleGeneratorLattice::ParticleGeneratorLattice()
		: ParticleGenerator(), lower_bound_(0), upper_bound_(0),
		body_shape_(NULL), lattice_spacing_(0), number_of_lattices_(0)
	{
	}
	//=================================================================================================//
	void ParticleGeneratorLattice::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		sph_body_->getSPHSystemBound(lower_bound_, upper_bound_);
		body_shape_ = sph_body_->body_shape_;
		lattice_spacing_ = sph_body_->particle_spacing_;
		CalcNumberOfLattices(lower_bound_, upper_bound_, lattice_spacing_);
	}
	//=================================================================================================//
	void ParticleGeneratorLattice
		::CalcNumberOfLattices(Vecd lower_bound, Vecd upper_bound,
			Real lattice_spacing)
	{
		Vecd zero(0);
		for (int i = 0; i < zero.size(); ++i) {
			number_of_lattices_[i] = static_cast<int>(ceil((upper_bound[i] 
				- lower_bound[i]) / lattice_spacing));
		}
	}
	//=================================================================================================//
	ParticleGeneratorRegularized::ParticleGeneratorRegularized()
		: ParticleGeneratorLattice()
	{

	}
	//=================================================================================================//
}