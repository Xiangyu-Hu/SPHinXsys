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
	FluidParticleData::FluidParticleData()
		: p_(0.0), rho_0_(1.0), rho_n_(1.0), mass_(1.0),
		drho_dt_(0.0), dvel_dt_trans_(0), vel_trans_(0),
		vorticity_(0) 
		
	{

	}
	//=================================================================================================//
	FluidParticleData::FluidParticleData(BaseParticleData &base_particle_data, Fluid *fluid)
		: p_(0.0), rho_0_(fluid->ReinitializeRho(p_)), rho_n_(rho_0_),
		mass_(base_particle_data.Vol_* rho_0_),
		drho_dt_(0.0), dvel_dt_trans_(0), vel_trans_(0),
		vorticity_(0)
	{
	
	}
	//=================================================================================================//
	FluidParticles::FluidParticles(SPHBody *body, Fluid* fluid)
		: BaseParticles(body, fluid), signal_speed_max_(0.0)
	{
		fluid->assignFluidParticles(this);
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			fluid_particle_data_.push_back(FluidParticleData(base_particle_data_[i], fluid));
		}
	}
	//=================================================================================================//
	FluidParticles* FluidParticles::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void FluidParticles::AddABufferParticle()
	{
		BaseParticles::AddABufferParticle();
		fluid_particle_data_.push_back(FluidParticleData());

	}
	//=================================================================================================//
	void FluidParticles::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		BaseParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		fluid_particle_data_[this_particle_index] = fluid_particle_data_[another_particle_index];
	}
	//=================================================================================================//
	void FluidParticles::UpdateFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		BaseParticles::UpdateFromAnotherParticle(this_particle_index, another_particle_index);
		fluid_particle_data_[this_particle_index].rho_n_ = fluid_particle_data_[another_particle_index].rho_n_;
		fluid_particle_data_[this_particle_index].p_ = fluid_particle_data_[another_particle_index].p_;
	}
	//=================================================================================================//
	void FluidParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		BaseParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(fluid_particle_data_[this_particle_index], fluid_particle_data_[that_particle_index]);
	}
	//=================================================================================================//
	void  FluidParticles
		::mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction)
	{
		BaseParticles::mirrorInAxisDirection(particle_index_i, body_bound, axis_direction);
		FluidParticleData& fluid_particle_data_i = fluid_particle_data_[particle_index_i];
		fluid_particle_data_i.vel_trans_[axis_direction] *= -1.0;
	}
	//=================================================================================================//
	void FluidParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		BaseParticles::WriteParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"Vorticity\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fluid_particle_data_[i].vorticity_[0] << " " << fluid_particle_data_[i].vorticity_[1] << " " << fluid_particle_data_[i].vorticity_[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		FluidParticles::WriteParticlesToVtuFile(output_file);
	}
	//=================================================================================================//
	void FluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		const SimTK::String xml_name("particles_xml"), ele_name("particles");
		unique_ptr<XmlEngine> restart_xml(new XmlEngine(xml_name, ele_name));

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<Real>("Sigma0", base_particle_data_[i].sigma_0_);
			restart_xml->AddAttributeToElement<Vecd>("Position", base_particle_data_[i].pos_n_);
			restart_xml->AddAttributeToElement<Real>("Volume", base_particle_data_[i].Vol_);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", base_particle_data_[i].vel_n_);
			restart_xml->AddAttributeToElement<Real>("Density", fluid_particle_data_[i].rho_n_);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void FluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Real sigma_0 = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Sigma0");
			base_particle_data_[number_of_particles].sigma_0_ = sigma_0;
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
	//=================================================================================================//
	ViscoelasticFluidParticleData::ViscoelasticFluidParticleData()
		: tau_(0), dtau_dt_(0)
	{

	}
	//=================================================================================================//
	ViscoelasticFluidParticles
		::ViscoelasticFluidParticles(SPHBody *body, Oldroyd_B_Fluid* oldroyd_b_fluid)
		: FluidParticles(body, oldroyd_b_fluid)
	{
		oldroyd_b_fluid->assignViscoelasticFluidParticles(this);
		for (size_t i = 0; i < base_particle_data_.size(); ++i) {
			viscoelastic_particle_data_.push_back(ViscoelasticFluidParticleData());
		}
	}
	//=================================================================================================//
	ViscoelasticFluidParticles* ViscoelasticFluidParticles::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::AddABufferParticle()
	{
		FluidParticles::AddABufferParticle();
		viscoelastic_particle_data_.push_back(ViscoelasticFluidParticleData());
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles
		::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		FluidParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		viscoelastic_particle_data_[this_particle_index] = viscoelastic_particle_data_[another_particle_index];
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		FluidParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(viscoelastic_particle_data_[this_particle_index], viscoelastic_particle_data_[that_particle_index]);
	}
	//=================================================================================================//
	void ViscoelasticFluidParticles::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));
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
	//=================================================================================================//
	void ViscoelasticFluidParticles::ReadParticleFromXmlForRestart(std::string &filefullpath)
	{
		cout << "\n This function ViscoelasticFluidParticles::ReadParticleFromXmlForRestart is not done. Exit the program! \n";
		exit(0);
	}
	//=================================================================================================//
}