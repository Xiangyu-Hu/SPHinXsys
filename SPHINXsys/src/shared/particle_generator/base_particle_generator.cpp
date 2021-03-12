/**
 * @file 	base_particle_generator.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "base_particle_generator.h"

#include "base_body.h"
#include "base_particles.h"
#include "in_output.h"


namespace SPH {
	//=================================================================================================//
	void ParticleGenerator::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
	}
	//=================================================================================================//
	void ParticleGeneratorDirect
		::createBaseParticles(BaseParticles* base_particles)
	{
		size_t total_real_particles = 0;
		auto& body_input_points_volumes = sph_body_->body_input_points_volumes_;
		for (size_t i = 0; i < body_input_points_volumes.size(); ++i)
		{
			base_particles->initializeABaseParticle(body_input_points_volumes[i].first,
				body_input_points_volumes[i].second);
			total_real_particles++;
		}

		base_particles->total_real_particles_ = total_real_particles;
	}
	//=================================================================================================//
	ParticleGeneratorReload::
		ParticleGeneratorReload(In_Output* in_output, std::string reload_body_name)
		: ParticleGenerator()
	{
		if (!fs::exists(in_output->reload_folder_))
		{
			std::cout << "\n Error: the particle reload folder:" << in_output->reload_folder_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		file_path_ = in_output->reload_folder_ + "/SPHBody_" + reload_body_name + "_rld.xml";
	}
	//=================================================================================================//
	void ParticleGeneratorReload::createBaseParticles(BaseParticles* base_particles)
	{
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(file_path_);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			base_particles->initializeABaseParticle(position, volume);
			total_real_particles++;
		}
		base_particles->total_real_particles_ = total_real_particles;
	}
	//=================================================================================================//
}
