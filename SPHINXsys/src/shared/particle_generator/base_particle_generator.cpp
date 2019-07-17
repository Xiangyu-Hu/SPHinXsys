/**
 * @file 	base_particle_generator.cpp
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume. The lattice generator generate
 * 			at lattice position by check whether the poision is contained by a SPH body.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "base_particle_generator.h"
#include "base_body.h"
#include "base_particles.h"
#include "xml_engine.h"

namespace SPH {
	//---------------------------------------------------------------//
	ParticleGenerator
		::ParticleGenerator(SPHBody &sph_body)
	{

	}
	//---------------------------------------------------------------//
	ParticleGeneratorDirect
		::ParticleGeneratorDirect(SPHBody &sph_body)
		: ParticleGenerator(sph_body)
	{

	}
	//---------------------------------------------------------------//
	void ParticleGeneratorDirect
		::CreateParticles(SPHBody &sph_body)
	{
		size_t number_of_particles = 0;
		for (int i = 0; i < sph_body.body_input_points_volumes_.size(); ++i) 
		{
			sph_body.GenerateAParticle(sph_body.body_input_points_volumes_[i].first,
				sph_body.body_input_points_volumes_[i].second);
			number_of_particles++;
		}

		sph_body.number_of_real_particles_ = number_of_particles;
	}
	//---------------------------------------------------------------//
	ReadRelaxedParticlsFromXmlFile
		::ReadRelaxedParticlsFromXmlFile(SPHBody &sph_body)
		: ParticleGenerator(sph_body)
	{

	}
	//---------------------------------------------------------------//
	void ReadRelaxedParticlsFromXmlFile
		::CreateParticles(SPHBody &sph_body)
	{
		std::string input_folder_ = "../input";
		std::string filefullpath = input_folder_ + "/SPHBody_" + sph_body.GetBodyName() + "_" + std::to_string(0) + ".xml";
		if (!fs::exists(filefullpath))
		{
			std::cout << "\n Error: the input file:"<< filefullpath << " is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}else{
			size_t number_of_particles = 0;
			XmlEngine* read_xml = new XmlEngine();
			read_xml->LoadXmlFile(filefullpath);
			SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
			for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
			{
					Vecd pos = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Position");
					Real vol = read_xml->GetRequiredAttributeRealValue(ele_ite_, "Volume");
					sph_body.GenerateAParticle(pos,vol);
					number_of_particles++;
			}
			sph_body.number_of_real_particles_ = number_of_particles;
		}
	}
	//---------------------------------------------------------------//
	ReadRestartParticlsFromXmlFile
		::ReadRestartParticlsFromXmlFile(SPHBody &sph_body)
		: ParticleGenerator(sph_body)
	{

	}
	//---------------------------------------------------------------//
	void ReadRestartParticlsFromXmlFile
		::CreateParticles(SPHBody &sph_body)
	{
		std::string input_folder_ = "./rstfile";
		std::string filefullpath = input_folder_ + "/SPHBody_" + sph_body.GetBodyName() + "_rst_" + std::to_string(sph_body.rst_step_) + ".xml";
		if (!fs::exists(filefullpath))
		{
			std::cout << "\n Error: the input file:"<< filefullpath << " is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}else{
			size_t number_of_particles = 0;
			XmlEngine* read_xml = new XmlEngine();
			read_xml->LoadXmlFile(filefullpath);
			SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
			for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
			{
					Vecd pos = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Position");
					Real vol = read_xml->GetRequiredAttributeRealValue(ele_ite_, "Volume");
					sph_body.GenerateAParticle(pos,vol);
					number_of_particles++;
			}
			sph_body.number_of_real_particles_ = number_of_particles;
		}
	}
	//---------------------------------------------------------------//
}