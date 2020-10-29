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
	//=================================================================================================//
	FluidParticles::FluidParticles(SPHBody *body, Fluid* fluid)
		: BaseParticles(body, fluid)
	{
		fluid->assignFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(p_, registered_scalars_, scalars_map_, scalars_to_write_, "Pressure", false);
		registerAVariable(drho_dt_, registered_scalars_, scalars_map_, scalars_to_write_, "DensityChangeRate", false);
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		sortable_vectors_.push_back(&pos_n_);
		sortable_vectors_.push_back(&vel_n_);
		sortable_scalars_.push_back(&mass_);
		sortable_scalars_.push_back(&rho_n_);
		sortable_scalars_.push_back(&p_);
	}
	//=================================================================================================//
	FluidParticles* FluidParticles::pointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		FluidParticles::writeParticlesToVtuFile(output_file);
	}
	//=================================================================================================//
	void FluidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		unique_ptr<XmlEngine> restart_xml(new XmlEngine(xml_name, ele_name));

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", vel_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Density", rho_n_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void FluidParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[number_of_particles] = pos_;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[number_of_particles] = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			vel_n_[number_of_particles] = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			rho_n_[number_of_particles] = rst_rho_n_;
			dvel_dt_[number_of_particles] = Vecd(0);
			number_of_particles++;
		}
	}
	//=================================================================================================//
	ViscoelasticFluidParticles
		::ViscoelasticFluidParticles(SPHBody *body, Oldroyd_B_Fluid* oldroyd_b_fluid)
		: FluidParticles(body, oldroyd_b_fluid)
	{
		oldroyd_b_fluid->assignViscoelasticFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(tau_, registered_matrices_, matrices_map_, matrices_to_write_, "ElasticStress", true);
		registerAVariable(dtau_dt_, registered_matrices_, matrices_map_, matrices_to_write_, "ElasticStressChangeRate", false);
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		sortable_matrices_.push_back(&tau_);
	}
	//=================================================================================================//
	ViscoelasticFluidParticles* ViscoelasticFluidParticles::pointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));
		size_t number_of_particles = pos_n_.size();
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function ViscoelasticFluidParticles::ReadParticleFromXmlForRestart is not done. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
}