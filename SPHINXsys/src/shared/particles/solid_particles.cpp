/**
 * @file solid_body_particles.cpp
 * @brief Definition of funcitons declared in solid_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 * @version 0.2.1
 * 			add muscle particles and muscle data.
 */
#include "solid_particles.h"
#include "elastic_solid.h"
//=============================================================================================//
namespace SPH {
//=============================================================================================//
	SolidParticleData::SolidParticleData(Vecd position)
		: pos_0_(position), n_0_(0), n_(0), B_(0), vel_ave_(0), dvel_dt_ave_(0),
		viscous_force_from_fluid_(0), force_from_fluid_(0)
	{

	}
//=============================================================================================//
	ElasticSolidParticleData::ElasticSolidParticleData(BaseParticleData &base_particle_data,
		ElasticSolid *elastic_solid)
		:F_(1.0), dF_dt_(0), stress_(0), mass_(1.0), pos_temp_(0)
	{
		rho_0_ = elastic_solid->rho_0_;
		rho_n_ = rho_0_;
		mass_ = rho_0_ *base_particle_data.Vol_;
	}
//=============================================================================================//
	MuscleParticleData::MuscleParticleData()
		: voltage_n_(0.0), grad_voltage_(0.0),
		dvoltage_dt_(0), gate_var_(0), T_a_(0.0), active_stress_(Matd(0.0))
	{

	}
//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody *body)
		: Particles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i) 
		{
			Point pnt = base_particle_data_[i].pos_n_;
			solid_body_data_.push_back(SolidParticleData(pnt));
		}
	}
	//=============================================================================================//
	void SolidParticles::OffsetInitialParticlePosition(Vecd offset)
	{
		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			base_particle_data_[i].pos_n_ += offset;
			solid_body_data_[i].pos_0_ += offset;
		}
	}
//=============================================================================================//
	void SolidParticles::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = position;
			solid_body_data_[number_of_particles].pos_0_ = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
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
//=============================================================================================//
	SolidParticles* SolidParticles::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void SolidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			restart_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void SolidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		/** Nothing should be done for non-moving BCs. */
	}
	//=============================================================================================//
	ElasticSolidParticles::ElasticSolidParticles(SPHBody *body)
		: SolidParticles(body)
	{
		elastic_solid_ = dynamic_cast<ElasticSolid*>(body->base_material_);
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			elastic_body_data_.push_back(ElasticSolidParticleData(base_particle_data_[i], elastic_solid_));
	}
//=============================================================================================//
	ElasticSolidParticles* ElasticSolidParticles::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			restart_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			restart_xml->AddAttributeToElement<Real>("Density", elastic_body_data_[i].rho_n_);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", base_particle_data_[i].vel_n_);
			restart_xml->AddAttributeToElement<Vecd>("Displacement", elastic_body_data_[i].pos_temp_);
			restart_xml->AddAttributeToElement("DefTensor", elastic_body_data_[i].F_);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ElasticSolidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			base_particle_data_[number_of_particles].pos_n_ = pos_;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			base_particle_data_[number_of_particles].Vol_ = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			base_particle_data_[number_of_particles].vel_n_ = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			elastic_body_data_[number_of_particles].rho_n_ = rst_rho_n_;
			Vecd rst_dis_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Displacement");
			elastic_body_data_[number_of_particles].pos_temp_ = rst_dis_;
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			elastic_body_data_[number_of_particles].F_ = rst_F_;
			number_of_particles++;
		}
	}
	//=============================================================================================//
	MuscleParticles::MuscleParticles(SPHBody *body)
		: ElasticSolidParticles(body)
	{
		muscle_ = dynamic_cast<Muscle*>(body->base_material_);
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
		{
			muscle_body_data_.push_back(MuscleParticleData());
		}
	}
//===============================================================//
	MuscleParticles* MuscleParticles::PointToThisObject()
	{
		return this;
	}
//=============================================================================================//
}
