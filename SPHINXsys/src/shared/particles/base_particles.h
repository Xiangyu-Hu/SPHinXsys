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
		/** Total number for all, including the buffer, particles. */
		static size_t total_number_of_particles_;		
	public:
		/** Default constructor. */
		BaseParticleData();
		/** In this constructor, the particle state is set at rest. */
		BaseParticleData(Vecd position, Real Vol_0, Real sigma_0);
		virtual ~BaseParticleData() {};

		/** Particle ID 
		 *	@brief For a real particle, it is the particle index.
		 *	For a ghost particle, it is the index of its corresponding real particle.
		 */
		size_t particle_id_;
		/** Current position. */
		Point pos_n_, pos_0_;	
		/** Current particle velocity and stress-induced and other accelerations. */
		Vecd  vel_n_, dvel_dt_, dvel_dt_others_;
		/** Particle volume and its reference volume. */
		Real Vol_, Vol_0_;
		/** Particle reference number density. */
		Real sigma_0_;
		/** smoothing length of the particle. */
		Real smoothing_length_;
		/** particle with fixed index, not subject to sorting. */
		bool is_sortable_;
	};

	/**
	 * @class BaseParticles
	 * @brief A group of base particles with essential (geometric and kinematic) data.
	 * There are three types of particles. 
	 * One is real particles whose states are updated by particle dynmaics. 
	 * One is buffer particles whose state are not updated by particle dynmaics. 
	 * They may be switched from real particles or switch to real particles. 
	 * The other is ghost particles whose states are updated according to 
	 * boundary condition if their indices are inlcuded in the neigbor particle list.   
	 */
	class BaseParticles
	{
	protected:
		/** The body in which the particles belongs to. */
		SPHBody *body_;
		string body_name_;
	public:
		/** Base material corresponding to base particles*/
		BaseMaterial* base_material_;

		BaseParticles(SPHBody *body, BaseMaterial* base_material);
		BaseParticles(SPHBody* body);
		virtual ~BaseParticles() {};
	
		/** Vector of base particle data. */
		StdLargeVec<BaseParticleData> base_particle_data_;	
		
		//----------------------------------------------------------------------
		//Global information for all partiles
		//----------------------------------------------------------------------
		Real speed_max_;		/**< Maxium particle speed. */
		/** Maxium possible number of real particles. 
		  * Also the start index of ghost particles. */
		size_t real_particles_bound_;
		size_t number_of_ghost_particles_;
		
		/** Initialize a base prticle by input a postion, volume
		  * and reference number density. */
		void InitializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0);
		/** Add buffer particles which latter may be realized for particle dynamics,
		  * or used as ghost particle. */
		virtual void AddABufferParticle();
		/** Copy state, except particle id, from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index);
		/** Update the state of a particle from another particle */
		virtual void UpdateFromAnotherParticle(size_t this_particle_index, size_t another_particle_index);
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index);
		/** Check whether partcles allowed for swaping*/
		bool allowSwapping(size_t this_particle_index, size_t that_particle_index);
		/** Getinsert a ghost particle. */
		size_t insertAGhostParticle(size_t index_particle_i);

		/** acess the sph body*/
		SPHBody* getSPHBody();
		/** Write particle data in VTU format for Paraview. */
		virtual void WriteParticlesToVtuFile(ofstream &output_file);
		/** Write particle data in PLT format for Tecplot. */
		virtual void WriteParticlesToPltFile(ofstream& output_file) {};

		/** Write particle data in XML format for restart. */
		virtual void WriteParticlesToXmlForRestart(std::string& filefullpath) {};
		/** Initialize particle data from restart xml file. */
		virtual void ReadParticleFromXmlForRestart(std::string& filefullpath) {};

		/** Output particle position and volume in XML file for reloading particles. */
		virtual void WriteToXmlForReloadParticle(std::string &filefullpath);
		/** Reload particle position and volume from XML files. */
		virtual void ReadFromXmlForReloadParticle(std::string &filefullpath);

		/** Pointer to this object. */
		virtual BaseParticles* PointToThisObject();

		/** Access a real data*/
		virtual Real accessAParticleDataTypeReal(size_t particel_index) { return 0.0; };
		/** Access a matrix data*/
		virtual Matd accessAParticleDataTypeMatd(size_t particel_index) { return Matd(0.0); };

		/** Normalize a gradient. */
		virtual Vecd normalizeGradient(size_t particle_index_i, Vecd& gradient) { return gradient; };
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j,
			Real dW_ij, Vecd& e_ij) { return dW_ij * e_ij; };
		/** Get mirror a particle along an axis direaction. */
		virtual void mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction);

	};
}
