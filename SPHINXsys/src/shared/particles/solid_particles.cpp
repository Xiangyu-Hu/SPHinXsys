/**
 * @file solid_body_particles.cpp
 * @brief Definition of funcitons declared in solid_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 * @version 0.2.1
 * 			add muscle particles and muscle data.
 */
#include "base_body.h"
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
		rho_0_ = elastic_solid->getReferenceDensity();
		rho_n_ = rho_0_;
		mass_ = rho_0_ *base_particle_data.Vol_;
	}
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body)
		: BaseParticles(body, new BaseMaterial())
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
		{
			Point pnt = base_particle_data_[i].pos_n_;
			solid_body_data_.push_back(SolidParticleData(pnt));
		}
	}
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body, BaseMaterial* base_material)
		: BaseParticles(body, base_material)
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
	//=================================================================================================//
	void SolidParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		BaseParticles::WriteParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"NormalDirection\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << solid_body_data_[i].n_[0] << " " << solid_body_data_[i].n_[1] << " " << 0.0 << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
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
	//===============================================================//
	void SolidParticles
		::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		BaseParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		solid_body_data_[this_particle_index] = solid_body_data_[another_particle_index];
	}
	//===============================================================//
	void SolidParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		BaseParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(solid_body_data_[this_particle_index], solid_body_data_[that_particle_index]);
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
	//=================================================================================================//	
	Vecd SolidParticles::normalizeGradient(size_t particle_index_i, Vecd& gradient) 
	{
		Matd&   B_i = solid_body_data_[particle_index_i].B_;
		return  B_i * gradient;
	}
	//=================================================================================================//
	Vecd SolidParticles::getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) 
	{
		Matd& B_i = solid_body_data_[particle_index_i].B_;
		Matd& B_j = solid_body_data_[particle_index_j].B_;
		return 0.5 * dW_ij * (B_i + B_j) * e_ij;
	}
	//=============================================================================================//
	ElasticSolidParticles::ElasticSolidParticles(SPHBody* body, BaseMaterial* base_material)
		: SolidParticles(body, base_material)
	{
		ElasticSolid* elastic_solid = dynamic_cast<ElasticSolid*>(base_material_);
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			elastic_body_data_.push_back(ElasticSolidParticleData(base_particle_data_[i], elastic_solid));
	}
	//===============================================================//
	void ElasticSolidParticles
		::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		SolidParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		elastic_body_data_[this_particle_index] = elastic_body_data_[another_particle_index];
	}
	//===============================================================//
	void ElasticSolidParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		SolidParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(elastic_body_data_[this_particle_index], elastic_body_data_[that_particle_index]);
	}
	//=============================================================================================//
	ElasticSolidParticles* ElasticSolidParticles::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ElasticSolidParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		SolidParticles::WriteParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"Density\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << elastic_body_data_[i].rho_n_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
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
	ActiveMuscleParticles::ActiveMuscleParticles(SPHBody* body, BaseMaterial* base_material)
		: ElasticSolidParticles(body, base_material)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			active_muscle_data_.push_back(ActiveMuscleData());
	}
	//=============================================================================================//
	void ActiveMuscleParticles
		::CopyFromAnotherParticle(size_t this_particle_index, size_t another_particle_index)
	{
		ElasticSolidParticles::CopyFromAnotherParticle(this_particle_index, another_particle_index);
		active_muscle_data_[this_particle_index] = active_muscle_data_[another_particle_index];
	}
	//=============================================================================================//
	void ActiveMuscleParticles::swapParticles(size_t this_particle_index, size_t that_particle_index)
	{
		ElasticSolidParticles::swapParticles(this_particle_index, that_particle_index);
		std::swap(active_muscle_data_[this_particle_index], active_muscle_data_[that_particle_index]);
	}
	//=============================================================================================//
	ActiveMuscleParticles* ActiveMuscleParticles::PointToThisObject()
	{
		return this;
	}
	//=============================================================================================//
	void ActiveMuscleParticles::WriteParticlesToVtuFile(ofstream& output_file)
	{
		ElasticSolidParticles::WriteParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"Active Stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << active_muscle_data_[i].active_contraction_stress_ << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ActiveMuscleParticles::WriteParticlesToXmlForRestart(std::string& filefullpath)
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
			restart_xml->AddAttributeToElement<Real>("ActiveStress", active_muscle_data_[i].active_contraction_stress_);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ActiveMuscleParticles::ReadParticleFromXmlForRestart(std::string& filefullpath)
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
			Real rst_active_stress_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "ActiveStress");
			active_muscle_data_[number_of_particles].active_contraction_stress_ = rst_active_stress_;
			number_of_particles++;
		}
	}
	//=================================================================================================//
}
