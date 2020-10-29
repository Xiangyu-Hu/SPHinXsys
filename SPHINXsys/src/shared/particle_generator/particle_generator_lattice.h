/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	particle_generator_lattice.h
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume. The lattice generator generate
 * 			at lattice position by check whether the poision is contained by a SPH body.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_particle_generator.h"
#include "geometry.h"
#include "base_mesh.h"

namespace SPH {
	/**
	 * @class ParticleGeneratorLattice
	 * @brief generate particles from lattice positions for a body.
	 */
	class ParticleGeneratorLattice : public ParticleGenerator
	{
	public:
		ParticleGeneratorLattice();
		virtual ~ParticleGeneratorLattice() {};

		virtual void initialize(SPHBody* sph_body) override;
		virtual void CreateBaseParticles(BaseParticles* base_particles) override;
	protected:
		Real lattice_spacing_;		/**< Lattice size. */
		Vecd lower_bound_, upper_bound_;	/**< Domain bounds. */
		std::unique_ptr<Mesh> mesh_;
		ComplexShape* body_shape_;
	};
}
