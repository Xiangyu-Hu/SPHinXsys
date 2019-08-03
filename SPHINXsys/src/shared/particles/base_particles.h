/**
 * @file 	base_particle.h
 * @brief 	This is the base class of SPH particles. The basic data of the particles
 *			is saved in a large vector. Each derived class will introduce an extra     
 * 			vector for the new data. Note that there is no class of particle.
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"
#include "xml_engine.h"

#include <fstream>
using namespace std;

namespace SPH {

	/** preclaimed classes*/
	class SPHBody;
	class ParticleGenerator;

	/**
	 * @class BaseParticleData
	 * @brief A based particle with essential data.
	  */
	class BaseParticleData
	{
	protected:
		/** Maximum ID number for all particles. */
		static int ID_max_;		

	public:
		BaseParticleData(Vecd position, Real volume);
		virtual ~BaseParticleData() {};

		/** Particle ID. */
		int particle_ID_;			
		/** Cell location of a particle, used building inner configuration. */
		Vecu cell_location_;		
		/** Current position. */
		Point pos_n_;	
		/** Current particle velocity and stress induced and other accelerations. */
		Vecd  vel_n_, dvel_dt_, dvel_dt_others_;	
		/** Particle volume and initial volume. */
		Real Vol_, Vol_0_;							
	};
	/**
	 * @class Particles
	 * @brief Agroup of based particles with essential data.
	 */
	class Particles
	{
	protected:
		/** The body in which the particles belongs to. */
		SPHBody *body_;
		string body_name_;

		ParticleGenerator *particle_generator_;
	public:
		Particles(SPHBody *body);
		virtual ~Particles() {};
		/** Vector of base particle data. */
		StdLargeVec<BaseParticleData> base_particle_data_;	
		/** Initialize a base prticle by input a postion and volume. */
		void InitializeABaseParticle(Vecd pnt, Real particle_volume);

		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) = 0;
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream &output_file) = 0;

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) = 0;
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string &filefullpath) = 0;

		/** Output particle position and volume in XML file for reloading particles. */
		virtual void WriteToXmlForReloadParticle(std::string &filefullpath);
		/** Reload particle position and volume from XML files. */
		virtual void ReadFromXmlForReloadParticle(std::string &filefullpath);

		/** Pointer to this object. */
		virtual Particles* PointToThisObject();
	};
}
