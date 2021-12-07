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

#ifndef BASE_PARTICLES_H
#define BASE_PARTICLES_H



#include "base_data_package.h"
#include "sph_data_containers.h"
#include "xml_engine.h"

#include <fstream>

namespace SPH
{

	class SPHBody;
	class BaseMaterial;
	class ShapeSurface;

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
		BaseParticles(SPHBody* body, BaseMaterial* base_material);
		BaseParticles(SPHBody* body);
		virtual ~BaseParticles() {};

		BaseMaterial* base_material_; /**< for dynamic cast in particle data delegation */

		StdLargeVec<Vecd> pos_n_;	/**< current position */
		StdLargeVec<Vecd> vel_n_;	/**< current particle velocity */
		StdLargeVec<Vecd> dvel_dt_;	/**< inner pressure- or stress-induced acceleration */
		StdLargeVec<Vecd> dvel_dt_prior_; /**<  other, such as gravity and viscous, accelerations */

		StdLargeVec<Real> Vol_;		/**< particle volume */
		StdLargeVec<Real> rho_n_;	/**< current particle density */
		StdLargeVec<Real> mass_;	/**< particle mass */
		//----------------------------------------------------------------------
		//Global information for all particles
		//----------------------------------------------------------------------
		Real rho0_;			/**< reference density*/
		Real sigma0_;			/**< reference number density. */
		Real speed_max_;		/**< Maximum particle speed. */
		Real signal_speed_max_; /**< Maximum signal speed.*/
		//----------------------------------------------------------------------
		//Global information for defining particle groups
		//----------------------------------------------------------------------
		size_t total_real_particles_;
		size_t real_particles_bound_; /**< Maximum possible number of real particles. Also the start index of ghost particles. */
		size_t total_ghost_particles_;
		//----------------------------------------------------------------------
		//		Generalized particle data for parameterized management
		//----------------------------------------------------------------------
		ParticleData	all_particle_data_;
		ParticleDataMap all_variable_maps_;

		/** register a variable defined in a class (can be non-particle class) to base particles */
		template<int DataTypeIndex, typename VariableType>
		void registerAVariable(StdLargeVec<VariableType>& variable_addrs,
			std::string variable_name, VariableType initial_value = VariableType(0))
		{
			if (all_variable_maps_[DataTypeIndex].find(variable_name) == all_variable_maps_[DataTypeIndex].end())
			{
				variable_addrs.resize(real_particles_bound_, initial_value);
				std::get<DataTypeIndex>(all_particle_data_).push_back(&variable_addrs);
				all_variable_maps_[DataTypeIndex].insert(make_pair(variable_name, std::get<DataTypeIndex>(all_particle_data_).size() - 1));
			}
			else
			{
				std::cout << "\n Error: the variable '" << variable_name << "' has already been registered!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		};

		/** Create and register a new variable which has not been defined yet in particles.
		 * If the variable is registered already, the registered variable will be returned. */
		template<int DataTypeIndex, typename VariableType>
		StdLargeVec<VariableType>* createAVariable(std::string new_variable_name, VariableType initial_value = VariableType(0))
		{
			if (all_variable_maps_[DataTypeIndex].find(new_variable_name) == all_variable_maps_[DataTypeIndex].end()) {
				StdLargeVec<VariableType>* new_variable = new StdLargeVec<VariableType>;
				registerAVariable<DataTypeIndex, VariableType>(*new_variable, new_variable_name, initial_value);
				return new_variable;
			}
			std::cout << "\n Warning: the variable '" << new_variable_name << "' has already registered!";
			std::cout << "\n So, the previously registered variable is assigned!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			return getVariableByName<DataTypeIndex, VariableType>(new_variable_name);
		};

		/** create and register a new variable, which has not been defined yet in particles,
		 *  by copying data from an exist variable */
		template<int DataTypeIndex, typename VariableType>
		StdLargeVec<VariableType>* createAVariable(std::string new_variable_name, std::string old_variable_name)
		{
			if (all_variable_maps_[DataTypeIndex].find(old_variable_name) != all_variable_maps_[DataTypeIndex].end())
			{
				if (all_variable_maps_[DataTypeIndex].find(new_variable_name) == all_variable_maps_[DataTypeIndex].end())
				{
					StdLargeVec<VariableType>* new_variable = new StdLargeVec<VariableType>;
					registerAVariable<DataTypeIndex, VariableType>(*new_variable, new_variable_name);
					StdLargeVec<VariableType>* old_variable =
						std::get<DataTypeIndex>(all_particle_data_)[all_variable_maps_[DataTypeIndex][old_variable_name]];
					for (size_t i = 0; i != real_particles_bound_; ++i) (*new_variable)[i] = (*old_variable)[i];
					return new_variable;
				}
				std::cout << "\n Error: the new variable '" << new_variable_name << "' has already been registered!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			std::cout << "\n Error: the old variable '" << old_variable_name << "' is not registered!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
			return nullptr;
		};

		/** get a registered variable from particles by its name */
		template<int DataTypeIndex, typename VariableType>
		StdLargeVec<VariableType>* getVariableByName(std::string variable_name)
		{
			if (all_variable_maps_[DataTypeIndex].find(variable_name) != all_variable_maps_[DataTypeIndex].end())
				return std::get<DataTypeIndex>(all_particle_data_)[all_variable_maps_[DataTypeIndex][variable_name]];

			std::cout << "\n Error: the variable '" << variable_name << "' is not registered!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
			return nullptr;
		};

		/** add a variable into a particle vairable name list */
		template<int DataTypeIndex, typename VariableType>
		void addAVariableNameToList(ParticleVariableList& variable_name_list, std::string variable_name)
		{
			if (all_variable_maps_[DataTypeIndex].find(variable_name) != all_variable_maps_[DataTypeIndex].end())
			{
				bool is_to_add = true;
				for (size_t i = 0; i != variable_name_list[DataTypeIndex].size(); ++i) {
					if (variable_name_list[DataTypeIndex][i].first == variable_name) is_to_add = false;
				}
				if (is_to_add) {
					size_t variable_index = all_variable_maps_[DataTypeIndex][variable_name];
					variable_name_list[DataTypeIndex].push_back(make_pair(variable_name, variable_index));
				}
			}
			else
			{
				std::cout << "\n Error: the variable '" << variable_name << "' you are going to write is not particle data!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		};

		/** add a variable into the list for state output */
		template<int DataTypeIndex, typename VariableType>
		void addAVariableToWrite(std::string variable_name)
		{
			addAVariableNameToList<DataTypeIndex, VariableType>(variables_to_write_, variable_name);
		};

		/** add a variable into the list for restart */
		template<int DataTypeIndex, typename VariableType>
		void addAVariableToRestart(std::string variable_name)
		{
			addAVariableNameToList<DataTypeIndex, VariableType>(variables_to_restart_, variable_name);
		};

		//----------------------------------------------------------------------
		//		Particle data for sorting
		//----------------------------------------------------------------------
		StdLargeVec<size_t> sequence_;
		StdLargeVec<size_t> sorted_id_;
		StdLargeVec<size_t> unsorted_id_;
		ParticleData	sortable_data_;
		ParticleDataMap sortable_variable_maps_;

		/** register an already defined variable as sortable */
		template<int DataTypeIndex, typename VariableType>
		void registerASortableVariable(std::string variable_name)
		{
			if (sortable_variable_maps_[DataTypeIndex].find(variable_name) == sortable_variable_maps_[DataTypeIndex].end())
			{
				if (all_variable_maps_[DataTypeIndex].find(variable_name) != all_variable_maps_[DataTypeIndex].end())
				{
					StdLargeVec<VariableType>* variable =
						std::get<DataTypeIndex>(all_particle_data_)[all_variable_maps_[DataTypeIndex][variable_name]];
					std::get<DataTypeIndex>(sortable_data_).push_back(variable);
					sortable_variable_maps_[DataTypeIndex].insert(make_pair(variable_name, std::get<DataTypeIndex>(sortable_data_).size() - 1));
				}
				else
				{
					std::cout << "\n Error: the variable '" << variable_name << "' is not registered!" << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
			}
			else
			{
				std::cout << "\n Warning: the variable '" << variable_name << "' has already registered as a sortabele variable!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			}
		};

		SPHBody* getSPHBody() { return body_; };
		void initializeABaseParticle(Vecd pnt, Real Vol_0);
		void addBufferParticles(size_t buffer_size);
		void copyFromAnotherParticle(size_t this_index, size_t another_index);
		void updateFromAnotherParticle(size_t this_index, size_t another_index);
		size_t insertAGhostParticle(size_t index_i);
		void switchToBufferParticle(size_t index_i);

		/** Write particle data in VTU format for Paraview. */
		virtual void writeParticlesToVtuFile(std::ostream& output_file);
		/** Write only surface particle data in VTU format for Paraview. */
		virtual void writeSurfaceParticlesToVtuFile(std::ostream& output_file, ShapeSurface& surface_particles);
		/** Write particle data in PLT format for Tecplot. */
		void writeParticlesToPltFile(std::ofstream& output_file);

		void resizeXmlDocForParticles(XmlEngine& xml_engine);
		void writeParticlesToXmlForRestart(std::string& filefullpath);
		void readParticleFromXmlForRestart(std::string& filefullpath);
		XmlEngine* getReloadXmlEngine() { return &reload_xml_engine_; };
		void writeToXmlForReloadParticle(std::string& filefullpath);
		void readFromXmlForReloadParticle(std::string& filefullpath);

		virtual BaseParticles* ThisObjectPtr() { return this; };

		/** Normalize the kernel gradient. */
		virtual Vecd normalizeKernelGradient(size_t particle_index_i, Vecd& kernel_gradient) { return kernel_gradient; };
		/** Get the kernel gradient in weak form. */
		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j,
			Real dW_ij, Vecd& e_ij) {
			return dW_ij * e_ij;
		};
	protected:
		SPHBody* body_; /**< The body in which the particles belongs to. */
		std::string body_name_;
		XmlEngine restart_xml_engine_;
		XmlEngine reload_xml_engine_;
		ParticleVariableList variables_to_write_;
		ParticleVariableList variables_to_restart_;
		void addAParticleEntry();

		virtual void writePltFileHeader(std::ofstream& output_file);
		virtual void writePltFileParticleData(std::ofstream& output_file, size_t index_i);

		template<int DataTypeIndex, typename VariableType>
		struct addAParticleDataValue
		{
			void operator () (ParticleData& particle_data) const
			{
				for (size_t i = 0; i != std::get<DataTypeIndex>(particle_data).size(); ++i)
					std::get<DataTypeIndex>(particle_data)[i]->push_back(VariableType(0));
			};
		};

		template<int DataTypeIndex, typename VariableType>
		struct copyAParticleDataValue
		{
			void operator () (ParticleData& particle_data, size_t this_index, size_t another_index) const
			{
				for (size_t i = 0; i != std::get<DataTypeIndex>(particle_data).size(); ++i)
					(*std::get<DataTypeIndex>(particle_data)[i])[this_index] =
					(*std::get<DataTypeIndex>(particle_data)[i])[another_index];
			};
		};
	};

	struct WriteAParticleVariableToXml
	{
		XmlEngine& xml_engine_;
		size_t& total_real_particles_;
		WriteAParticleVariableToXml(XmlEngine& xml_engine, size_t& total_real_particles) :
			xml_engine_(xml_engine), total_real_particles_(total_real_particles) {};
		template<typename VariableType>
		void operator () (std::string& variable_name, StdLargeVec<VariableType>& variable)  const
		{
			SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
			for (size_t i = 0; i != total_real_particles_; ++i)
			{
				xml_engine_.setAttributeToElement(ele_ite, variable_name, variable[i]);
				ele_ite++;
			}
		}
	};

	struct ReadAParticleVariableFromXml
	{
		XmlEngine& xml_engine_;
		size_t& total_real_particles_;
		ReadAParticleVariableFromXml(XmlEngine& xml_engine, size_t& total_real_particles) :
			xml_engine_(xml_engine), total_real_particles_(total_real_particles) {};
		template<typename VariableType>
		void operator () (std::string& variable_name, StdLargeVec<VariableType>& variable)  const
		{
			SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
			for (size_t i = 0; i != total_real_particles_; ++i)
			{
				xml_engine_.getRequiredAttributeValue(ele_ite, variable_name, variable[i]);
				ele_ite++;
			}
		}
	};
}
#endif //BASE_PARTICLES_H