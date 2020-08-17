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

namespace SPH {
	/**
	 * @class ParticleGeneratorLattice
	 * @brief generate particles from lattice positions for a body.
	 */
	class ParticleGeneratorLattice : public ParticleGenerator
	{
	public:
		ParticleGeneratorLattice(SPHBody& sph_body);
		virtual ~ParticleGeneratorLattice() {};

		/** Compute reference number density*/
		virtual Real ComputeReferenceNumberDensity();
		/** Create lattice particle for a body. */
		virtual void CreateBaseParticles(BaseParticles* base_particles) override;

	protected:
		Vecd lower_bound_, upper_bound_;	/**< Domain bounds. */
		ComplexShape &body_shape_;
		Real lattice_spacing_;		/**< Lattice size. */
		Vecu number_of_lattices_;	/**< Number of lattice. */ 
		/**
		 * @brief Calculate the number of Lattices.
		 * @param[in] lower_bound Lower bound of lattice size.
		 * @param[in] upuper_bound Upper bound of lattice size.
		 * @param[in] lattice_spacing Lattice size.
		 */
		void CalcNumberOfLattices(Vecd lower_bound, Vecd upper_bound, Real lattice_spacing);
	};

	/**
	 * @class ParticleGeneratorRegularized
	 * @brief generate particles from lattice positions for a body.
	 */
	class ParticleGeneratorRegularized : public ParticleGeneratorLattice
	{
	public:
		ParticleGeneratorRegularized(SPHBody& sph_body);
		virtual ~ParticleGeneratorRegularized() {};

		/** Create lattice particle for a body. */
		virtual void CreateBaseParticles(BaseParticles* base_particles) override {};
	};
}