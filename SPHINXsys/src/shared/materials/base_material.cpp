/**
 * @file 	base_material.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "base_material.h"
#include "base_particles.h"
#include "base_particles.hpp"

namespace SPH
{
	//=================================================================================================//
	void BaseMaterial::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	void BaseMaterial::writeToXmlForReloadLocalParameters(const std::string &filefullpath)
	{
		std::cout << "\n Material properties writing. " << std::endl;
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleData &all_particle_data = base_particles_->all_particle_data_;
		base_particles_->resizeXmlDocForParticles(reload_material_xml_engine_);
		WriteAParticleVariableToXml
			write_variable_to_xml(reload_material_xml_engine_, total_real_particles);
		ParticleDataOperation<loopVariabaleNameList> loop_variable_namelist;
		loop_variable_namelist(all_particle_data, reload_local_parameters_, write_variable_to_xml);
		reload_material_xml_engine_.writeToXmlFile(filefullpath);
		std::cout << "\n Material properties writing finished. " << std::endl;
	}
	//=================================================================================================//
	void BaseMaterial::readFromXmlForLocalParameters(const std::string &filefullpath)
	{
		reload_material_xml_engine_.loadXmlFile(filefullpath);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleData &all_particle_data = base_particles_->all_particle_data_;
		ReadAParticleVariableFromXml
			read_variable_from_xml(reload_material_xml_engine_, total_real_particles);
		ParticleDataOperation<loopVariabaleNameList> loop_variable_namelist;
		loop_variable_namelist(all_particle_data, reload_local_parameters_, read_variable_from_xml);

		if (total_real_particles != reload_material_xml_engine_.SizeOfXmlDoc())
		{
			std::cout << "\n Error: reload material properties does not match!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		else
		{
			std::cout << "\n Material properties reading finished." << std::endl;
		}
	}
	//=================================================================================================//
}
