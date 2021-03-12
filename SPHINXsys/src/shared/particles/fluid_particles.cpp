/** 
 * @file fluid_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "fluid_particles.h"

#include "base_body.h"
#include "weakly_compressible_fluid.h"
#include "xml_engine.h"

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
		registerAVariable<indexScalar, Real>(p_, "Pressure");
		registerAVariable<indexScalar, Real>(drho_dt_, "DensityChangeRate");
		registerAVariable<indexScalar, Real>(rho_sum_, "DensitySummation");
		registerAVariable<indexBoolean, bool>(is_free_surface_, "FreeSurfaceIndication");
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
	void FluidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		unique_ptr<XmlEngine> restart_xml(new XmlEngine(xml_name, ele_name));

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
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
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[total_real_particles] = pos_;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[total_real_particles] = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			vel_n_[total_real_particles] = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			rho_n_[total_real_particles] = rst_rho_n_;
			dvel_dt_[total_real_particles] = Vecd(0);
			total_real_particles++;
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
		registerAVariable<indexMatrix, Matd>(tau_, "ElasticStress");
		registerAVariable<indexMatrix, Matd>(dtau_dt_, "ElasticStressChangeRate");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		sortable_matrices_.push_back(&tau_);
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));
		size_t total_real_particles = pos_n_.size();
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
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
		cout << "\n This function ViscoelasticFluidParticles::ReadParticleFromXmlForRestart is not done! \n";
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
	}
	//=================================================================================================//
}
