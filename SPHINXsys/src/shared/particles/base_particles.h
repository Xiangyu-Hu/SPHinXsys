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
#include <fstream>

using namespace std;

namespace SPH {
	/**
	 * @class BaseParticleData
	 * @brief A based particle with essential data.
	  */
	class BaseParticleData
	{
	protected:
		static int ID_max_;		/**< Maximum ID number for all particles. */

	public:
		/**
		 * @brief Default constrcutor.
		 * @detail Create a group of particle refered to a body.
		 */
		BaseParticleData() {};
		/**
		 * @brief Constrcutor.
		 * @detail Create a particle.
		 * @param[in] positioni Particle positiion.
		 * @param[in] volum Particle volume.
		 */
		BaseParticleData(Vecd position, Real volume);
		/**
		 * @brief Default destructor.
		 */
		virtual ~BaseParticleData() {};

		int particle_ID_;			/**< Particle ID. */
		Vecu cell_location_;		/**< Cell location of a particle, */
		Point pos_n_;	/**< Current position. */
		/** Current and transport particle velocity and acceleration. **/
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
		/** The name of the body in which the particles belongs to. */
		string body_name_;		
	public:
		/**
		 * @brief Default Constructor.
		 * @detail Create a group of particles referred to a body.
		 * @param[in] body_name Name of a body.
		 */
		Particles(string body_name);
		/**
		 * @brief Destructor.
		 */
		virtual ~Particles() {};

		StdLargeVec<BaseParticleData> base_particle_data_;	/**< Vector of base particle data. */
		size_t number_of_particles_;						/**< Number of particles. */
		/**
		 * @brief Initialize a prticle by input a postion and volume. 
		 * @param[in] pnt Vecotor of particle position.
		 * @param[in] particle_volume Volume of particle.
		 */
		virtual void InitializeAParticle(Vecd pnt, Real particle_volume) = 0;
		/**
		 * @brief Write particle data in VTU format for Paraview.
		 * @param[out] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToVtuFile(ofstream &output_file) = 0;
		/**
		 * @brief Write particle data in PLT format for Tecplot.
		 * @param[out] output_file Ofstream of particle data.
		 */
		virtual void WriteParticlesToPltFile(ofstream &output_file) = 0;
		/**
		 * @brief Write particle data in XML format.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlFile(std::string &filefullpath) = 0;
		/**
		 * @brief Write particle data in XML format for restart.
		 * @param[out] filefullpath Full path to file being write.
		 */
		virtual void WriteParticlesToXmlForRestart(std::string &filefullpath) = 0;
		/**
		 * @brief Initialize particle data from restart xml file.
		 */
		virtual void InitialParticleFromRestartXmlFile(std::string &filefullpath) = 0;
		/**
		 * @brief Pointer to this object. 
		 */
		virtual Particles* PointToThisObject();
	};
}
