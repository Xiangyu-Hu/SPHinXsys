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
	 * @brief generate particles from lattice poistions for a body.
	 */
	class ParticleGeneratorLattice : public ParticleGenerator
	{
		Region &region_;
		Real lattice_spacing_;		/**< Lattice size. */
		Vecu number_of_lattices_;	/**< Number of lattice. */ 
		Real sigma_ref_;
		/**
		 * @brief Claculate the number of Lattices.
		 * @param[in] lower_bound Lower bound of lattice size.
		 * @param[in] upuper_bound Upper bound of lattice size.
		 * @param[in] lattice_spacing Lattice size.
		 */
		void CalcNumberOfLattices(Vecd lower_bound, Vecd upper_bound, Real lattice_spacing);
	protected:
		Vecd lower_bound_, upper_bound_;	/**< Domain bounds. */
	public:
		ParticleGeneratorLattice(SPHBody &sph_body, Region &region);
		virtual ~ParticleGeneratorLattice() {};

		/** Compute reference number density*/
		virtual Real ComputeReferenceNumberDensity();
		/** Create lattice particle for a body. */
		virtual void CreateBaseParticles() override;
	};
}