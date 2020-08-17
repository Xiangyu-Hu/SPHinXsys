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
	ParticleGeneratorLattice
		::ParticleGeneratorLattice(SPHBody &sph_body)
		: ParticleGenerator(sph_body), body_shape_(sph_body.getBodyShape())
	{
		sph_body.getSPHSystemBound(lower_bound_, upper_bound_);
		lattice_spacing_ = sph_body.particle_spacing_;
		CalcNumberOfLattices(lower_bound_, upper_bound_, lattice_spacing_);
		Vecd body_lower_bound, bound_upper_bound;
		sph_body.findBodyShapeBounds(body_lower_bound, bound_upper_bound);
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
	ParticleGeneratorRegularized::ParticleGeneratorRegularized(SPHBody& sph_body)
		: ParticleGeneratorLattice(sph_body)
	{

	}
	//=================================================================================================//
}