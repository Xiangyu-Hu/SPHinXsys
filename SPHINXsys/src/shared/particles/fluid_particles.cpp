/**
 * @file fluid_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"
#include "base_body.h"

namespace SPH
{
	//===============================================================//
	FluidParticleData::FluidParticleData()
		: p_(0.0), drho_dt_(0.0), div_correction_(1.0), vel_trans_(0),
		dvel_dt_trans_(0), vorticity_(0), vort_2d_(0.0),
		rho_0_(1.0), rho_n_(1.0), mass_(1.0), dvel_dt_inner_(0)
	{

	}
	//===============================================================//
	FluidParticleData::FluidParticleData(BaseParticleData &base_particle_data, 
		WeaklyCompressibleFluid *weakly_compressible_fluid)
		: p_(0.0), drho_dt_(0.0), div_correction_(1.0), vel_trans_(0),
		dvel_dt_trans_(0), vorticity_(0), vort_2d_(0.0), dvel_dt_inner_(0)
	{
		rho_0_ = weakly_compressible_fluid->ReinitializeRho(p_);
		rho_n_ = rho_0_;
		mass_ = base_particle_data.Vol_ * rho_0_;
	}
	//===============================================================//
	FluidParticles::FluidParticles(SPHBody *body) : Particles(body)
	{
		weakly_compressible_fluid_ = dynamic_cast<WeaklyCompressibleFluid*>(body->base_material_);
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			fluid_particle_data_.push_back(FluidParticleData(base_particle_data_[i], weakly_compressible_fluid_));
		}
	}
	//===============================================================//
	FluidParticles* FluidParticles::PointToThisObject() 
	{
		return this;
	}
	//===============================================================//
	void FluidParticles::AddABufferParticle()
	{
		Particles::AddABufferParticle();
		fluid_particle_data_.push_back(FluidParticleData());

	}
	//===============================================================//
	void FluidParticles::RealizeABufferParticle(size_t buffer_particle_index, size_t real_particle_index)
	{
		Particles::RealizeABufferParticle(buffer_particle_index, real_particle_index);
		fluid_particle_data_[buffer_particle_index] = fluid_particle_data_[real_particle_index];
		FluidParticleData& fluid_particle_data_i
			= fluid_particle_data_[buffer_particle_index];
	}
	//===============================================================//
	void FluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		XmlEngine* restart_xml = new XmlEngine(xml_name, ele_name);

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			restart_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", base_particle_data_[i].vel_n_);
			restart_xml->AddAttributeToElement<Real>("Density", fluid_particle_data_[i].rho_n_);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//===============================================================//
	void FluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
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
			fluid_particle_data_[number_of_particles].rho_n_ = rst_rho_n_;
			fluid_particle_data_[number_of_particles].vel_trans_(0);
			base_particle_data_[number_of_particles].dvel_dt_(0);
			number_of_particles++;
		}
	}
	//===============================================================//
	ViscoelasticFluidParticleData::ViscoelasticFluidParticleData()
		: tau_(0), dtau_dt_(0)
	{

	}
	//===============================================================//
	ViscoelasticFluidParticles::ViscoelasticFluidParticles(SPHBody *body)
		: FluidParticles(body)
	{
		oldroyd_b_fluid_ = dynamic_cast<Oldroyd_B_Fluid*>(body->base_material_);
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			viscoelastic_particle_data_.push_back(ViscoelasticFluidParticleData());
		}
	}
	//===============================================================//
	ViscoelasticFluidParticles* ViscoelasticFluidParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
	void ViscoelasticFluidParticles::AddABufferParticle()
	{
		FluidParticles::AddABufferParticle();
		viscoelastic_particle_data_.push_back(ViscoelasticFluidParticleData());
	}
	//===============================================================//
	void ViscoelasticFluidParticles
		::RealizeABufferParticle(size_t buffer_particle_index, size_t real_particle_index)
	{
		FluidParticles::RealizeABufferParticle(buffer_particle_index, real_particle_index);
		viscoelastic_particle_data_[buffer_particle_index] = viscoelastic_particle_data_[real_particle_index];
	}
	//===============================================================//
	void ViscoelasticFluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		XmlEngine* restart_xml = new XmlEngine("particles_xml", "particles");
		size_t number_of_particles = base_particle_data_.size();
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
	//===============================================================//
	void ViscoelasticFluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function ViscoelasticFluidParticles::ReadParticleFromXmlForRestart is not done. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
}