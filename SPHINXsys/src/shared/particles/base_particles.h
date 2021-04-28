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
 */
#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"

#include <fstream>
using namespace std;

namespace SPH {

	class SPHBody;
	class BaseMaterial;

	/**
	 * @class BaseParticles
	 * @brief Particles with essential (geometric and kinematic) data.
	 * There are three types of particles， all particles of a same type are saved with continuous memory segments. 
	 * The first type is real particles whose states are updated by particle dynamics. 
	 * One is buffer particles whose state are not updated by particle dynamics. 
	 * Buffer particles are saved behind real particles.
	 * The global value of total_real_particles_ separate the real and buffer particles.
	 * They may be switched from real particles or switch to real particles. 
	 * As the memory for both particles are continuous, such switch is achieved at the memory boundary sequentially.
	 * The basic idea is swap the data of the last real particle with the one will be switched particle, 
	 * and then switch this swapped last particle as buffer particle by decrease the total_real_particles_ by one.
	 * Switch from buffer particle to real particle is easy. One just need to assign expect state to 
	 * the first buffer particle and increase total_real_particles_ by one.    
	 * The other is ghost particles whose states are updated according to 
	 * boundary condition if their indices are included in the neighbor particle list.
	 * The ghost particles are saved behind the buffer particles.
	 * The global value of real_particles_bound_ separate the sum of real and buffer particles with ghost particles.
	 * The global value of total_ghost_particles_ indicates the total number of ghost particles in use.
	 * It will be initialized to zero before a time step. 
	 */
	class BaseParticles
	{
	public:
		BaseParticles(SPHBody *body, BaseMaterial* base_material);
		BaseParticles(SPHBody* body);
		virtual ~BaseParticles() {};
	
		BaseMaterial* base_material_; /**< for dynamic cast in particle data delegation */

		StdLargeVec<Vecd> pos_n_;	/**< current position */
		StdLargeVec<Vecd> vel_n_;	/**< current particle velocity */
		StdLargeVec<Vecd> dvel_dt_;	/**< inner pressure- or stress-induced acceleration */
		StdLargeVec<Vecd> dvel_dt_others_; /**<  other, such as gravity and viscous, accelerations */

		StdLargeVec<Real> Vol_;		/**< particle volume */
		StdLargeVec<Real> rho_n_;	/**< current particle density */
		StdLargeVec<Real> mass_;	/**< particle mass */
		StdLargeVec<Real> h_ratio_;	/**< the ratio between reference smoothing length to variable smoothing length */
		//----------------------------------------------------------------------
		//Global information for all particles
		//----------------------------------------------------------------------
		Real rho_0_;			/**< reference density*/
		Real sigma_0_;			/**< reference number density. */
		Real speed_max_;		/**< Maximum particle speed. */
		Real signal_speed_max_; /**< Maximum signal speed.*/
		//----------------------------------------------------------------------
		//Global information for defining partilce groups
		//----------------------------------------------------------------------
		size_t total_real_particles_;
		size_t real_particles_bound_; /**< Maximum possible number of real particles. Also the start index of ghost particles. */
		size_t total_ghost_particles_;
		//----------------------------------------------------------------------
		//		Generalized particle data for parameterized management
		//----------------------------------------------------------------------
		std::tuple<StdVec<StdLargeVec<Real>*>, StdVec<StdLargeVec<Vecd>*>, StdVec<StdLargeVec<Matd>*>,
			StdVec<StdLargeVec<int>*>, StdVec<StdLargeVec<bool>*>> all_particle_data_;
		std::array<map<string, size_t>, 5> name_index_maps_;
		std::array<StdVec<string>, 5> variables_to_write_;

		/** register a defined variable into particles */
		template<int DataTypeIndex, typename VariableType>
		void registerAVariable(StdLargeVec<VariableType>& variable_addrs,
			string variable_name, bool is_to_write = false, VariableType initial_value = VariableType(0))
		{
			if (name_index_maps_[DataTypeIndex].find(variable_name) == name_index_maps_[DataTypeIndex].end()) 
			{
				variable_addrs.resize(real_particles_bound_, initial_value);
				std::get<DataTypeIndex>(all_particle_data_).push_back(&variable_addrs);
				name_index_maps_[DataTypeIndex].insert(make_pair(variable_name, std::get<DataTypeIndex>(all_particle_data_).size() - 1));
				if (is_to_write) variables_to_write_[DataTypeIndex].push_back(variable_name);
			}
			else
			{
				std::cout << "\n Error: the variable '" << variable_name << "' has already been registered!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		};
		/** create and register a new variable, which has not been defined yet, into particles */
		template<int DataTypeIndex, typename VariableType>
		StdLargeVec<VariableType>* createAVariable(string variable_name, bool is_to_write = false, VariableType initial_value = VariableType(0))
		{
			if (name_index_maps_[DataTypeIndex].find(variable_name) == name_index_maps_[DataTypeIndex].end()) {
				StdLargeVec<VariableType>* variable = new StdLargeVec<VariableType>;
				variable->resize(real_particles_bound_, initial_value);
				std::get<DataTypeIndex>(all_particle_data_).push_back(variable);
				name_index_maps_[DataTypeIndex].insert(make_pair(variable_name, std::get<DataTypeIndex>(all_particle_data_).size() - 1));
				if (is_to_write) variables_to_write_[DataTypeIndex].push_back(variable_name);
				return variable;
			}
			std::cout << "\n Warning: the variable '" << variable_name << "' has already registered!" << std::endl;
			std::cout << "\n So, the previously registered variable is assigned!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			return getVariableByName<DataTypeIndex, VariableType>(variable_name);
		};

		/** get a registered variable from particles by its name */
		template<int DataTypeIndex, typename VariableType>
		StdLargeVec<VariableType>* getVariableByName(string variable_name) 
		{
			if (name_index_maps_[DataTypeIndex].find(variable_name) != name_index_maps_[DataTypeIndex].end())
				return std::get<DataTypeIndex>(all_particle_data_)[name_index_maps_[DataTypeIndex][variable_name]];

			std::cout << "\n Error: the variable '" << variable_name << "' is not registered!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
			return NULL;
		};

		/** add a variable into the list for state output */
		template<int DataTypeIndex, typename VariableType>
		void addAVariableToWrite(string variable_name)
		{
			if (name_index_maps_[DataTypeIndex].find(variable_name) != name_index_maps_[DataTypeIndex].end())
			{
				bool is_to_write = true;
				for (size_t i = 0; i != variables_to_write_[DataTypeIndex].size(); ++i) {
					if (variables_to_write_[DataTypeIndex][i] == variable_name) is_to_write = false;
				}
				if (is_to_write) variables_to_write_[DataTypeIndex].push_back(variable_name);
			}
			else
			{
				std::cout << "\n Error: the variable '" << variable_name << "' you are going to write is not particle data!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
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

		SPHBody* getSPHBody() { return body_; };
		void initializeABaseParticle(Vecd pnt, Real Vol_0);
		void addABufferParticle();
		void copyFromAnotherParticle(size_t this_index, size_t another_index);
		void updateFromAnotherParticle(size_t this_index, size_t another_index);
		size_t insertAGhostParticle(size_t index_i);
		void switchToBufferParticle(size_t index_i);

		/** Write particle data in VTU format for Paraview. */
		virtual void writeParticlesToVtuFile(ofstream &output_file);
		/** Write particle data in PLT format for Tecplot. */
		virtual void writeParticlesToPltFile(ofstream& output_file);
		virtual void writeParticlesToXmlForRestart(std::string& filefullpath) {};
		virtual void readParticleFromXmlForRestart(std::string& filefullpath) {};
		virtual void writeToXmlForReloadParticle(std::string &filefullpath);
		virtual void readFromXmlForReloadParticle(std::string &filefullpath);
		virtual BaseParticles* pointToThisObject() {return this;};

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
