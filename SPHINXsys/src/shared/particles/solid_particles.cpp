/**
 * @file solid_particles.cpp
 * @brief Definition of functions declared in solid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "solid_particles.h"

#include "geometry.h"
#include "base_body.h"
#include "elastic_solid.h"
#include "xml_engine.h"

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
		registerAVariable<indexVector, Vecd>(pos_0_, "InitialPosition");
		registerAVariable<indexVector, Vecd>(n_, "NormalDirection", true);
		registerAVariable<indexVector, Vecd>(n_0_, "InitialNormalDirection");
		registerAVariable<indexMatrix, Matd>(B_, "CorrectionMatrix", false, Matd(1.0));
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(vel_ave_, "AverageVelocity");
		registerAVariable<indexVector, Vecd>(dvel_dt_ave_, "AverageAcceleration");
		registerAVariable<indexVector, Vecd>(force_from_fluid_, "ForceFromFluid");
		registerAVariable<indexVector, Vecd>(viscous_force_from_fluid_, "ViscousForceFromFluid");
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable<indexScalar, Real>(contact_density_, "ContactDensity");
		registerAVariable<indexVector, Vecd>(contact_force_, "ContactForce");
		//-----------------------------------------------------------------------------------------
		//		register sortable particle data before building up particle configuration
		//-----------------------------------------------------------------------------------------
		sortable_vectors_.push_back(&pos_n_);
		sortable_vectors_.push_back(&pos_0_);
		sortable_scalars_.push_back(&Vol_);

		//set the initial value for intial particle position
		for (size_t i = 0; i != pos_n_.size(); ++i) pos_0_[i] =  pos_n_[i];

		//sorting particle once
		dynamic_cast<RealBody*>(body)->sortParticleWithMeshCellLinkedList();
	}
	//=============================================================================================//
	void SolidParticles::offsetInitialParticlePosition(Vecd offset)
	{
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			pos_n_[i] += offset;
			pos_0_[i] += offset;
		}
	}
	//=================================================================================================//
	void SolidParticles::initializeNormalDirectionFromGeometry()
	{
		ComplexShape* body_shape = body_->body_shape_;
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			Vecd normal_direction = body_shape->findNormalDirection(pos_n_[i]);
			n_[i] = normal_direction;
			n_0_[i] = normal_direction;
		}
	}
	//=============================================================================================//
	void SolidParticles::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd position = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[total_real_particles] = position;
			pos_0_[total_real_particles] = position;
			Real volume = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[total_real_particles] = volume;
			total_real_particles++;
		}

		if (total_real_particles != pos_n_.size())
		{
			std::cout << "\n Error: reload particle number does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	void SolidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
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
		registerAVariable<indexMatrix, Matd>(F_, "DeformationGradient", false, Matd(1.0));
		registerAVariable<indexMatrix, Matd>(dF_dt_, "DeformationRate");
		registerAVariable<indexMatrix, Matd>(stress_PK1_, "FirstPiolaKirchhoffStress");
		registerAVariable<indexMatrix, Matd>(corrected_stress_, "CorrectedStress");
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToVtuFile(ofstream& output_file)
	{
		SolidParticles::writeParticlesToVtuFile(output_file);

		size_t total_real_particles = total_real_particles_;

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << fixed << setprecision(9) << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
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
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[total_real_particles] = pos_;
			Vecd pos_0 = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "InitialPosition");
			pos_0_[total_real_particles] = pos_0;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[total_real_particles] = rst_Vol_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			vel_n_[total_real_particles] = rst_vel_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			rho_n_[total_real_particles] = rst_rho_n_;
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			F_[total_real_particles] = rst_F_;
			total_real_particles++;
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
		registerAVariable<indexMatrix, Matd>(active_stress_, "ActiveStress");
		registerAVariable<indexScalar, Real>(active_contraction_stress_, "ActiveContractionStress", true);
	}
	//=================================================================================================//
	void ActiveMuscleParticles::writeParticlesToXmlForRestart(std::string& filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
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
			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DefTensor");
			F_[total_real_particles] = rst_F_;
			Real rst_active_stress_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "ActiveStress");
			active_contraction_stress_[total_real_particles] = rst_active_stress_;
			total_real_particles++;
		}
	}
	//=============================================================================================//
	ShellParticles::ShellParticles(SPHBody* body, ElasticSolid* elastic_solid, Real thickness)
		: ElasticSolidParticles(body, elastic_solid)
	{
		elastic_solid->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexScalar, Real>(shell_thickness_, "Thickness", false, thickness);
		registerAVariable<indexVector, Vecd>(pseudo_n_, "InitialPseudoNormal");
		registerAVariable<indexVector, Vecd>(dpseudo_n_dt_, "PseudoNormalChangeRate");
		registerAVariable<indexVector, Vecd>(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
		registerAVariable<indexVector, Vecd>(rotation_, "Rotation", true);
		registerAVariable<indexVector, Vecd>(angular_vel_, "AngularVelocity");
		registerAVariable<indexVector, Vecd>(dangular_vel_dt_, "AngularAcceleration");
		registerAVariable<indexMatrix, Matd>(F_bending_, "BendingDeformationGradient");
		registerAVariable<indexMatrix, Matd>(dF_bending_dt_, "BendingDeformationGradientChangeRate");
		registerAVariable<indexMatrix, Matd>(resultant_stress_, "ResultantStress");
		registerAVariable<indexMatrix, Matd>(resultant_moment_, "ResultantMoment");
		registerAVariable<indexMatrix, Matd>(transformation_matrix_, "TransformationMatrix");
	}
	//=================================================================================================//
	void ShellParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		unique_ptr<XmlEngine> restart_xml(new XmlEngine("particles_xml", "particles"));

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			restart_xml->creatXmlElement("particle");
			restart_xml->AddAttributeToElement<Vecd>("Position", pos_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("InitialPosition", pos_0_[i]);
			restart_xml->AddAttributeToElement<Real>("Volume", Vol_[i]);
			restart_xml->AddAttributeToElement<Real>("Density", rho_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("Velocity", vel_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("PseudoNormal", pseudo_n_[i]);
			restart_xml->AddAttributeToElement<Vecd>("PseudoNormalChangeRate", dpseudo_n_dt_[i]);
			restart_xml->AddAttributeToElement<Real>("Thickness", shell_thickness_[i]);
			restart_xml->AddAttributeToElement("DeformationGradient", F_[i]);
			restart_xml->AddElementToXmlDoc();
		}
		restart_xml->WriteToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void ShellParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		size_t total_real_particles = 0;
		unique_ptr<XmlEngine> read_xml(new XmlEngine());
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd pos_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Position");
			pos_n_[total_real_particles] = pos_;
			Vecd pos_0 = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "InitialPosition");
			pos_0_[total_real_particles] = pos_0;
			Real rst_Vol_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Volume");
			Vol_[total_real_particles] = rst_Vol_;
			Real rst_rho_n_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Density");
			rho_n_[total_real_particles] = rst_rho_n_;
			Vecd rst_vel_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Velocity");
			vel_n_[total_real_particles] = rst_vel_;

			Vecd rst_pseudo_normal_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "PseudoNormal");
			pseudo_n_[total_real_particles] = rst_pseudo_normal_;
			Vecd rst_pseudo_normal_change_rate_ = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "PseudoNormalChangeRate");
			dpseudo_n_dt_[total_real_particles] = rst_pseudo_normal_change_rate_;
			Real rst_thickness_ = read_xml->GetRequiredAttributeValue<Real>(ele_ite_, "Thickness");
			shell_thickness_[total_real_particles] = rst_thickness_;

			Matd rst_F_ = read_xml->GetRequiredAttributeMatrixValue(ele_ite_, "DeformationGradient");
			F_[total_real_particles] = rst_F_;
			total_real_particles++;
		}
	}
	//=================================================================================================//
}
