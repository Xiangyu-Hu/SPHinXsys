/**
 * @file base_particles.cpp
 * @brief Definition of funcitons declared in base_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
 /**
  * @file 	base_particles.cpp
  * @author	Luhui Han, Chi ZHang and Xiangyu Hu
  * @version	0.1
  */

#include "base_particles.h"
#include "base_body.h"
#include "all_particle_generators.h"

namespace SPH
{
	//===============================================================//
	int BaseParticleData::total_number_of_particles_ = 0;
	//===============================================================//
	BaseParticleData::BaseParticleData()
		: vel_n_(0), pos_n_(0), Vol_(0), Vol_0_(0), sigma_(0), sigma_0_(0),
		cell_location_(0), dvel_dt_others_(0), dvel_dt_(0)
	{
		total_number_of_particles_++;
	}
	//===============================================================//
	BaseParticleData::BaseParticleData(Vecd position, Real Vol_0, Real sigma_0) 
		: vel_n_(0), pos_n_(position), Vol_(Vol_0), Vol_0_(Vol_0), sigma_(sigma_0), sigma_0_(sigma_0),
		cell_location_(0), dvel_dt_others_(0), dvel_dt_(0)
	{
		total_number_of_particles_++;
	}
	//===============================================================//
	Particles::Particles(SPHBody *body)
		: body_(body), body_name_(body->GetBodyName()),
		speed_max_(0.0)
	{
		body->base_particles_ = this;

		switch (body->particle_generator_op_)
		{
		case ParticlesGeneratorOps::lattice: {
			particle_generator_ = new ParticleGeneratorLattice(*body, body->body_region_);
			break;
		}

		case ParticlesGeneratorOps::direct: {
			particle_generator_ = new ParticleGeneratorDirect(*body);
			break;
		}

		default: {
			std::cout << "\n FAILURE: the type of particle generator is undefined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
			break;
		}
		}

		particle_generator_->CreateBaseParticles();
	}
	//===============================================================//
	void Particles::InitializeABaseParticle(Vecd pnt, Real Vol_0, Real sigma_0)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, Vol_0, sigma_0));
	}
	//===============================================================//
	void Particles::AddABufferParticle()
	{
		base_particle_data_.push_back(BaseParticleData());
	}
	//===============================================================//
	void Particles::RealizeABufferParticle(size_t buffer_particle_index, size_t real_particle_index)
	{
		base_particle_data_[buffer_particle_index] = base_particle_data_[real_particle_index];
	}
	//===============================================================//
	void Particles::WriteToXmlForReloadParticle(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		XmlEngine* reload_xml = new XmlEngine(xml_name, ele_name);

		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			reload_xml->CreatXmlElement("particle");
			reload_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			reload_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			reload_xml->AddElementToXmlDoc();
		}
		reload_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void Particles::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			base_particle_data_[number_of_particles].Vol_ = volume;
			number_of_particles++;
		}

		if(number_of_particles != base_particle_data_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//===============================================================//
	Particles* Particles::PointToThisObject()
	{
		return this;
	}
}