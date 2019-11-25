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
#include "sph_data_conainers.h"
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
		/** Total number for all real particles 
		  * which excludes the buffer particles. */
		static int total_number_of_particles_;		

	public:
		/** default constructor. */
		BaseParticleData();
		/** in this constructor, the particle is set at rest. */
		BaseParticleData(Vecd position, Real Vol_0, Real sigma_0);
		virtual ~BaseParticleData() {};

		/** Cell location of a particle, used for building inner configuration. */
		Vecu cell_location_;		
		/** Current position. */
		Point pos_n_;	
		/** Current particle velocity and stress induced and other accelerations. */
		Vecd  vel_n_, dvel_dt_, dvel_dt_others_;
		/** Particle volume and its reference volume. */
		Real Vol_, Vol_0_;
		/** Particle number density and its reference value. */
		Real sigma_, sigma_0_;
	};
	/**
	 * @class Particles
	 * @brief A group of based particles with essential (geometric and kinematic) data.
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
		
		//----------------------------------------------------------------------
		//Global data
		//----------------------------------------------------------------------
		Real speed_max_;		/**< Maxium particle speed. */
		
		/** Initialize a base prticle by input a postion, volume
		  * and reference number density. */
		void InitializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0);
		/** add buffer particles which latter may be realized for particle dynamics */
		virtual void AddABufferParticle();
		/** Realize a buffer particle from a real particle */
		virtual void RealizeABufferParticle(size_t buffer_particle_index, size_t real_particle_index);

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
