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
	 * @brief A based particle with essential data for all types of particles.
	  */
	class BaseParticleData
	{
	protected:
		/** All, including the buffer, particles. */
		static size_t total_number_of_particles_;		
	public:
		BaseParticleData();
		/** In this constructor, the particle state is set at rest. */
		BaseParticleData(Vecd position, Real Vol_0, Real sigma_0);
		virtual ~BaseParticleData() {};

		/** For a real particle, it is the index.
		 *	For a ghost particle, it is the index of its corresponding real particle. */
		size_t particle_id_;
		/** Current and initial position. */
		Point pos_n_, pos_0_;	
		/** Current particle velocity and stress-induced and other accelerations. */
		Vecd  vel_n_, dvel_dt_, dvel_dt_others_;
		/** Particle volume and its initial value. */
		Real Vol_, Vol_0_;
		/** Particle reference number density. */
		Real sigma_0_;
		/** Smoothing length of the particle. */
		Real smoothing_length_;
		/** Particle with fixed index, not subject to sorting. */
		bool is_sortable_;
	};

	/**
	 * @class BaseParticles
	 * @brief Particles with essential (geometric and kinematic) data.
	 * There are three types of particles. 
	 * One is real particles whose states are updated by particle dynamics. 
	 * One is buffer particles whose state are not updated by particle dynamics. 
	 * They may be switched from real particles or switch to real particles. 
	 * The other is ghost particles whose states are updated according to 
	 * boundary condition if their indices are included in the neighbor particle list.   
	 */
	class BaseParticles
	{
	protected:
		SPHBody *body_; /**< The body in which the particles belongs to. */
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
		Real speed_max_;		/**< Maximum particle speed. */
		Real signal_speed_max_; /**< Maximum signal speed.*/
		/** Maximum possible number of real particles. Also the start index of ghost particles. */
		size_t real_particles_bound_;
		size_t number_of_ghost_particles_;
		
		/** Initialize a base particle by input a postion, volume and reference number density. */
		void InitializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0);
		/** Add buffer particles which latter may be realized for particle dynamics, or used as ghost particle. */
		virtual void AddABufferParticle();
		/** Copy state, except particle id, from another particle */
		virtual void CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index);
		/** Update the state of a particle from another particle */
		virtual void UpdateFromAnotherParticle(size_t this_particle_index, size_t another_particle_index);
		/** Swapping particles. */
		virtual void swapParticles(size_t this_particle_index, size_t that_particle_index);
		/** Check whether particles allowed for swaping*/
		bool isSwappingAllowed(size_t this_particle_index, size_t that_particle_index);
		/** Insert a ghost particle into the particle list. */
		size_t insertAGhostParticle(size_t index_particle_i);

		/** access the sph body*/
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

		/** Normalize the kernel gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& kernel_gradient) { return kernel_gradient; };
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j,
			Real dW_ij, Vecd& e_ij) { return dW_ij * e_ij; };
		/** Get mirror a particle along an axis direaction. */
		virtual void mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction);

	};
}
