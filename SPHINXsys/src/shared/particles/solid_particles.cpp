/**
 * @file solid_particles.cpp
 * @brief Definition of funcitons declared in solid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 * @version 0.2.1
 * 			add muscle particles and muscle data.
 */
#include "solid_particles.h"
#include "base_body.h"
#include "elastic_solid.h"
//=============================================================================================//
namespace SPH {
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body)
		: SolidParticles(body, new Solid())
	{
	}
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body, Solid* solid)
		: BaseParticles(body, solid)
	{
		solid->assignSolidParticles(this);
	
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(pos_0_, registered_vectors_, vectors_map_, vectors_to_write_, "InitialPosition", false);
		registerAVariable(n_, registered_vectors_, vectors_map_, vectors_to_write_, "NormalDirection", true);
		registerAVariable(n_0_, registered_vectors_, vectors_map_, vectors_to_write_, "InitialNormalDirection", false); //seems to be moved to method
		registerAVariable(B_, registered_matrices_, matrices_map_, matrices_to_write_, "CorrectionMatrix", false, Matd(1.0));
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable(vel_ave_, registered_vectors_, vectors_map_, vectors_to_write_, "AverageVelocity", false);
		registerAVariable(dvel_dt_ave_, registered_vectors_, vectors_map_, vectors_to_write_, "AverageAcceleration", false);
		registerAVariable(force_from_fluid_, registered_vectors_, vectors_map_, vectors_to_write_, "ForceFromFluid", false);
		registerAVariable(viscous_force_from_fluid_, registered_vectors_, vectors_map_, vectors_to_write_, "ViscousForceFromFluid", false);
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable(contact_density_, registered_scalars_, scalars_map_, scalars_to_write_, "ContactDensity", true);
		registerAVariable(contact_force_, registered_vectors_, vectors_map_, vectors_to_write_, "ViscousForceFromFluid", false);

		//set the initial value
		for (size_t i = 0; i != pos_n_.size(); ++i) pos_0_[i] =  pos_n_[i];

	}
	//=============================================================================================//
	void SolidParticles::OffsetInitialParticlePosition(Vecd offset)
	{
		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			pos_n_[i] += offset;
			pos_0_[i] += offset;
		}
	}
	//=================================================================================================//
	void SolidParticles::initializeNormalDirectionFromGeometry()
	{
		for (size_t i = 0; i != body_->number_of_particles_; ++i)
		{
			Vecd normal_direction = body_->body_shape_->findNormalDirection(pos_n_[i]);
			n_[i] = normal_direction;
			n_0_[i] = normal_direction;
		}
	}
	//=============================================================================================//
	void SolidParticles::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[number_of_particles] = position;
			pos_0_[number_of_particles] = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[number_of_particles] = volume;
			number_of_particles++;
		}

		if (number_of_particles != pos_n_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=============================================================================================//
	SolidParticles* SolidParticles::pointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void SolidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("InitialPosition", pos_0_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void SolidParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		/** Nothing should be done for non-moving BCs. */
	}
	//=================================================================================================//	
	Vecd SolidParticles::normalizeKernelGradient(size_t particle_index_i, Vecd& kernel_gradient) 
	{
		return  B_[particle_index_i] * kernel_gradient;
	}
	//=================================================================================================//
	Vecd SolidParticles::getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) 
	{
		return 0.5 * dW_ij * (B_[particle_index_i] + B_[particle_index_j]) * e_ij;
	}
	//=============================================================================================//
	ElasticSolidParticles::ElasticSolidParticles(SPHBody* body, ElasticSolid* elastic_solid)
		: SolidParticles(body, elastic_solid)
	{
		elastic_solid->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(F_, registered_matrices_, matrices_map_, matrices_to_write_, "DeformationGradient", false, Matd(1.0));
		registerAVariable(dF_dt_, registered_matrices_, matrices_map_, matrices_to_write_, "DeformationRate", false);
		registerAVariable(stress_, registered_matrices_, matrices_map_, matrices_to_write_, "Stress", false);
	}
	//=============================================================================================//
	ElasticSolidParticles* ElasticSolidParticles::pointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		SolidParticles::writeParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("InitialPosition", pos_0_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", vel_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Density", rho_n_[i]);
			restart_xml->AddAttributeToElement("DefTensor", F_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ElasticSolidParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t number_of_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[number_of_particles] = pos_;
			Vecd pos_0 = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "InitialPosition");
			pos_0_[number_of_particles] = pos_0;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[number_of_particles] = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			vel_n_[number_of_particles] = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			rho_n_[number_of_particles] = rst_rho_n_;
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			F_[number_of_particles] = rst_F_;
			number_of_particles++;
		}
	}
	//=============================================================================================//
	ActiveMuscleParticles::ActiveMuscleParticles(SPHBody* body, ActiveMuscle* active_muscle)
		: ElasticSolidParticles(body, active_muscle)
	{
		active_muscle->assignActiveMuscleParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(active_stress_, registered_matrices_, matrices_map_, matrices_to_write_, "ActiveStress", false);
		registerAVariable(active_contraction_stress_, registered_scalars_, scalars_map_, scalars_to_write_, "ActiveContractionStress", true);
	}
	//=============================================================================================//
	ActiveMuscleParticles* ActiveMuscleParticles::pointToThisObject()
	{
		return this;
	}
	//=============================================================================================//
	void ActiveMuscleParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		ElasticSolidParticles::writeParticlesToVtuFile(output_file);

		size_t number_of_particles = body_->number_of_particles_;

		output_file << "    <DataArray Name=\"Active Stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != number_of_particles; ++i) {
			output_file << fixed << setprecision(9) << active_contraction_stress_[i] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ActiveMuscleParticles::writeParticlesToXmlForRestart(std::string& filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t number_of_particles = body_->number_of_particles_;
		for (size_t i = 0; i != number_of_particles; ++i)
		{
			restart_xml->CreatXmlElement("particle");
			restart_xml->AddAttributeToElement<size_t>("ID", i);
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", vel_n_[i]);
			restart_xml->AddAttributeToElement<Real>("Density", rho_n_[i]);
			restart_xml->AddAttributeToElement("DefTensor", F_[i]);
			restart_xml->AddAttributeToElement<Real>("ActiveStress", active_contraction_stress_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ActiveMuscleParticles::readParticleFromXmlForRestart(std::string& filefullpath)
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
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			F_[number_of_particles] = rst_F_;
			Real rst_active_stress_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "ActiveStress");
			active_contraction_stress_[number_of_particles] = rst_active_stress_;
			number_of_particles++;
		}
	}
	//=================================================================================================//
}
