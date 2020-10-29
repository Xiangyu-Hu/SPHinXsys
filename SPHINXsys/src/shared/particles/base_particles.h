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
 * @file 	base_particles.h
 * @brief 	This is the base class of SPH particles. The basic data of the particles
 *			is saved in separated large vectors. Each derived class will introduce several extra     
 * 			vectors for the new data. Note that there is no class of single particle.
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

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class SPHBody;
	class ParticleGenerator;

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
	public:
		/** Base material corresponding to base particles*/
		BaseMaterial* base_material_;

		BaseParticles(SPHBody *body, BaseMaterial* base_material);
		BaseParticles(SPHBody* body);
		virtual ~BaseParticles() {};
	
		StdLargeVec<Vecd> pos_n_;	/**< current position */
		StdLargeVec<Vecd> vel_n_;	/**< current particle velocity */
		StdLargeVec<Vecd> dvel_dt_;	/**< inner pressure- or stress-induced acceleration */
		StdLargeVec<Vecd> dvel_dt_others_; /**<  other, such as gravity and viscous acceleration */

		StdLargeVec<Real> Vol_;		/**< particle volume */
		StdLargeVec<Real> rho_n_;	/**< current particle density */
		StdLargeVec<Real> mass_;	/**< particle mass */
		StdLargeVec<Real> smoothing_length_;

		//----------------------------------------------------------------------
		//Global information for all partiles
		//----------------------------------------------------------------------
		Real rho_0_;			/**< global reference density*/
		Real sigma_0_;			/**< global reference number density. */
		Real speed_max_;		/**< Maximum particle speed. */
		Real signal_speed_max_; /**< Maximum signal speed.*/
		/** Maximum possible number of real particles. Also the start index of ghost particles. */
		size_t real_particles_bound_;
		size_t number_of_ghost_particles_;
		
		//----------------------------------------------------------------------
		//		Registered particle data
		//----------------------------------------------------------------------
		StdVec<StdLargeVec<Matd>*> registered_matrices_;
		StdVec<StdLargeVec<Vecd>*> registered_vectors_;
		StdVec<StdLargeVec<Real>*> registered_scalars_;
		map<string, size_t> matrices_map_; 	/**< Map from matrix names to indexes. */
		map<string, size_t> vectors_map_;	/**< Map from vector names to indexes. */
		map<string, size_t> scalars_map_;	/**< Map from scalar names to indexes. */
		StdVec<string> matrices_to_write_;		/**< matrix variables to be written in file */
		StdVec<string> vectors_to_write_;		/**< vactor variables to be written in file */
		StdVec<string> scalars_to_write_;		/**< scalar variables to be written in file */
		/** register a variable into particles */
		template<typename VariableType>
		void registerAVariable(StdLargeVec<VariableType>& variable_addrs,
			StdVec<StdLargeVec<VariableType>*>& registered_variables,
			map<string, size_t>& name_map, StdVec<string>& variable_to_write, 	
			string variable_name, bool is_to_write, VariableType initial_value = VariableType(0))
		{
			variable_addrs.resize(real_particles_bound_, initial_value);
			registered_variables.push_back(&variable_addrs);
			name_map.insert(make_pair(variable_name, registered_variables.size() - 1));
			if (is_to_write) variable_to_write.push_back(variable_name);
		};

		//----------------------------------------------------------------------
		//		Particle data for sorting
		//----------------------------------------------------------------------
		StdLargeVec<size_t> sequence_;
		StdLargeVec<size_t> sorted_id_;
		StdLargeVec<size_t> unsorted_id_;
		StdVec<StdLargeVec<Matd>*> sortable_matrices_;
		StdVec<StdLargeVec<Vecd>*> sortable_vectors_;
		StdVec<StdLargeVec<Real>*> sortable_scalars_;

		/** access the sph body*/
		SPHBody* getSPHBody() { return body_; };
		/** Initialize a base particle by input a postion, volume and reference number density. */
		void initializeABaseParticle(Vecd pnt, Real Vol_0);
		/** Add buffer particles which latter may be realized for particle dynamics, or used as ghost particle. */
		void addABufferParticle();
		/** Copy physical state from another particle */
		void copyFromAnotherParticle(size_t this_index, size_t another_index);
		/** Update physical state of a particle from another particle */
		void updateFromAnotherParticle(size_t this_index, size_t another_index);
		/** Insert a ghost particle into the particle list. */
		size_t insertAGhostParticle(size_t index_i);

		/** Write particle data in VTU format for Paraview. */
		virtual void writeParticlesToVtuFile(ofstream &output_file);
		/** Write particle data in PLT format for Tecplot. */
		virtual void writeParticlesToPltFile(ofstream& output_file) {};

		/** Write particle data in XML format for restart. */
		virtual void writeParticlesToXmlForRestart(std::string& filefullpath) {};
		/** Initialize particle data from restart xml file. */
		virtual void readParticleFromXmlForRestart(std::string& filefullpath) {};

		/** Output particle position and volume in XML file for reloading particles. */
		virtual void writeToXmlForReloadParticle(std::string &filefullpath);
		/** Reload particle position and volume from XML files. */
		virtual void readFromXmlForReloadParticle(std::string &filefullpath);

		/** Pointer to this object. */
		virtual BaseParticles* pointToThisObject();

		/** Normalize the kernel gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& kernel_gradient) { return kernel_gradient; };
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j,
			Real dW_ij, Vecd& e_ij) { return dW_ij * e_ij; };
	protected:
		SPHBody* body_; /**< The body in which the particles belongs to. */
		string body_name_;
	};
}
