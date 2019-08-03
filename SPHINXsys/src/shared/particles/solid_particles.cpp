/**
 * @file solid_body_particles.cpp
 * @brief Definition of funcitons declared in solid_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "solid_particles.h"

namespace SPH {
	//===============================================================//
	SolidParticleData::SolidParticleData(Vecd position)
		: pos_0_(position), n_0_(0), n_(0), B_(0), vel_ave_(0), dvel_dt_ave_(0),
		viscous_force_from_fluid_(0), force_from_fluid_(0)
	{

	}
	//===============================================================//
	ElasticSolidParticleData::ElasticSolidParticleData()
		: rho_0_(1.0), rho_n_(1.0), F_(1.0), dF_dt_(0), stress_(0), mass_(1.0), pos_temp_(0)
	{

	}
	//===============================================================//
	SolidParticles::SolidParticles(SPHBody *body)
		: Particles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			Point pnt = base_particle_data_[i].pos_n_;
			solid_body_data_.push_back(SolidParticleData(pnt));
		}
	}
	//===============================================================//
	void SolidParticles::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeVectorValue(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = position;
			solid_body_data_[number_of_particles].pos_0_ = position;
			Real volume = read_xml->GetRequiredAttributeRealValue(ele_ite_, "Volume");
			base_particle_data_[number_of_particles].Vol_ = volume;
			number_of_particles++;
		}

		if (number_of_particles != base_particle_data_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//===============================================================//
	SolidParticles* SolidParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
	ElasticSolidParticles::ElasticSolidParticles(SPHBody *body)
		: SolidParticles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			elastic_body_data_.push_back(ElasticSolidParticleData());
	}
	//===============================================================//
	ElasticSolidParticles* ElasticSolidParticles::PointToThisObject()
	{
		return this;
	}
}
