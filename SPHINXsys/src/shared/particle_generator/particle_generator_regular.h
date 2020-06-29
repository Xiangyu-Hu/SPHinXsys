/**
 * @file 	particle_generator_regular.h
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The regular generator generate
 * 			at lattice position by check whether the poision is contained by a 0 levelset.
 * @author	Yongchuan Yu and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_particle_generator.h"
#include "geometry.h"
#include "base_mesh.h"

namespace SPH {
	/**
	 * @class ParticleGeneratorLattice
	 * @brief generate particles from lattice poistions for a body.
	 */
	class ParticleGeneratorRegular : public ParticleGenerator
	{
		Region &region_;
		Real lattice_spacing_;		/**< Lattice size. */
		Vecu number_of_lattices_;	/**< Number of lattice. */ 
		Real sigma_ref_;
		Real phi_input_;
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
		ParticleGeneratorRegular(SPHBody &sph_body);
		virtual ~ParticleGeneratorRegular() {};

		/** Compute reference number density*/

		/** Compute reference number density*/
		virtual Real ComputeReferenceNumberDensity();
		/** Create lattice particle for a body. */
		virtual void CreateBaseParticles(BaseParticles *base_particles) override;
	};
}