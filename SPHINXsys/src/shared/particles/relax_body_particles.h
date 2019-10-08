/**
 * @file 	relax_body_particle.h
 * @brief 	This is the derived class of base particle.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_particles.h"
#include "xml_engine.h"

#include <fstream>
using namespace std;

namespace SPH {
	/**
	 * @class RelaxBodyParticleData
	 * @brief Data for relax body particles.
	 */
	class RelaxBodyParticleData
	{

	public:
		RelaxBodyParticleData(Vecd position);
		virtual ~RelaxBodyParticleData() {};
		
		/** Paticle postion after mapping particles to surface. */
		Vecd pos_0_;	
	};
	/**
	 * @class RelaxBodyParticles
	 * @brief A group of particles with relax body particle data.
	 */
	class RelaxBodyParticles : public Particles
	{
	public:
		RelaxBodyParticles(SPHBody *body);
		virtual ~RelaxBodyParticles() {};

		/** Vector of particle data. */
		StdLargeVec<RelaxBodyParticleData> relax_body_data_;	

		IndexVector lists_of_singularity_particles_;

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) override;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) override;

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) override;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) override;

		/** Pointer to this object. */
		virtual RelaxBodyParticles* PointToThisObject() override;
	};
}