/** 
 * @file elastic_solid.cpp
 * @author Chi Zhang and Xiangyu Hu
 * @version  0.1
 * @version  0.2.1
 */

#include "elastic_solid.h"

namespace SPH {
	//=================================================================================================//
	Real ElasticSolid::GetSoundSpeed(size_t particle_index_i)
	{
		return c_0_;
	}
	//=================================================================================================//
	Real ElasticSolid::getViscousTimeStepSize(Real smoothing_length)
	{
		return 0.5 * smoothing_length * smoothing_length * rho_0_ / (eta_0_ + 1.0e-15);
	}
	//=================================================================================================//
	Matd ElasticSolid::DampingStress(Matd& F, Matd& dF_dt, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (~dF_dt * F + ~F * dF_dt);
		Matd sigmaPK2 = eta_0_ * strain_rate;
		return F * sigmaPK2;
	}
	//=================================================================================================//
	Matd ElasticSolid::NumericalDampingStress(Matd& F, Matd& dF_dt, Real numerical_viscoisty)
	{
		Matd strain_rate = 0.5 * (~dF_dt * F + ~F * dF_dt);
		Matd sigmaPK2 = numerical_viscoisty * strain_rate;
		return F * sigmaPK2;
	}
	//=================================================================================================//
	Real ElasticSolid::getNumericalViscosity(Real smoothing_length)
	{
		return 0.5 * rho_0_ * c_0_ * smoothing_length;
	}
	//=================================================================================================//
	void LinearElasticSolid::assignDerivedMaterialParameters() 
	{
		Solid::assignDerivedMaterialParameters();
		lambda_0_ = SetLambda();
		G_0_ = SetShearModulus();
		c_0_ = SetSoundSpeed();
		std::cout << "The speed of sound: " << c_0_ << std::endl;
		std::cout << "The Lambda: " << lambda_0_ << std::endl;
	};
	//=================================================================================================//
	Matd LinearElasticSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = lambda_0_ * strain.trace() * Matd(1.0)
			+ 2.0 * G_0_ * strain;
		return F * sigmaPK2;
	}
	//=================================================================================================//
	Real LinearElasticSolid::SetSoundSpeed()
	{
		return  sqrt(E_0_ / 3.0 / (1.0 - 2.0 * nu_) / rho_0_);
	}
	//=================================================================================================//
	Real LinearElasticSolid::SetShearModulus()
	{
		return 0.5 * E_0_ / (1.0 + nu_);
	}
	//=================================================================================================//
	Real LinearElasticSolid::SetLambda()
	{
		return nu_ * E_0_ / (1.0 + nu_) / (1.0 - 2.0 * nu_);
	}
	//=================================================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = G_0_ * Matd(1.0)
			+ (lambda_0_ * log(det(F)) - G_0_) * inverse(right_cauchy);
		return F * sigmaPK2;
	}
	//=================================================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = G_0_ / (1.0 - 2.0 * strain.trace() / j1_m_) * Matd(1.0)
			+ (lambda_0_ * log(det(F)) - G_0_) * inverse(right_cauchy);
		return F * sigmaPK2;
	}
	//=================================================================================================//
	Real Muscle::SetSoundSpeed()
	{
		Real young_modulus = 2.0 * a_0_[0] * b_0_[0] * (1.0 + nu_);
		return  sqrt( young_modulus / 3.0 / (1.0 - 2.0 * nu_) / rho_0_);
	}
	//=================================================================================================//
	Real Muscle::SetLambda()
	{
		Real young_modulus = 2.0 * a_0_[0] * b_0_[0] * (1.0 + nu_);
		return nu_ * young_modulus / (1.0 + nu_) / (1.0 - 2.0 * nu_);
	}
	//=================================================================================================//
	void Muscle::assignDerivedMaterialParameters()
	{
		f0f0_ = SimTK::outer(f0_, f0_);
		f0s0_ = SimTK::outer(f0_, s0_);
		s0s0_ = SimTK::outer(s0_, s0_);
		c_0_ = SetSoundSpeed();
		lambda_0_ = SetLambda();
		std::cout << "The speed of sound: " << c_0_ << std::endl;
		std::cout << "The Lambda: " << lambda_0_ << std::endl;
	}
	//=================================================================================================//
	Matd Muscle::ConstitutiveRelation(Matd& F, size_t i)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy * f0_, f0_) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy * s0_, s0_) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy * f0_, s0_);
		Real ln_J = log(det(F));
		Real I_1_1 = right_cauchy.trace() - Real(f0_.size());
		Matd sigmaPK2 = a_0_[0] * exp(b_0_[0] * I_1_1) * Matd(1.0)
			+ (lambda_0_ * ln_J - a_0_[0]) * inverse(right_cauchy)
			+ 2.0 * a_0_[1] * I_ff_1 * exp(b_0_[1] * I_ff_1 * I_ff_1) * f0f0_
			+ 2.0 * a_0_[2] * I_ss_1 * exp(b_0_[2] * I_ss_1 * I_ss_1) * s0s0_
			+ a_0_[3] * I_fs * exp(b_0_[3] * I_fs * I_fs) * f0s0_;

		return F * sigmaPK2;
	}
	//=================================================================================================//
	Matd LocallyOrthotropicMuscle::ConstitutiveRelation(Matd& F, size_t i)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy * local_f0_[i], local_f0_[i]) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy * local_s0_[i], local_s0_[i]) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy * local_f0_[i], local_s0_[i]);
		Real ln_J = log(det(F));
		Real I_1_1 = right_cauchy.trace() - Real(Vecd(0).size());
		Matd sigmaPK2 = a_0_[0] * exp(b_0_[0] * I_1_1) * Matd(1.0)
			+ (lambda_0_ * ln_J - a_0_[0]) * inverse(right_cauchy)
			+ 2.0 * a_0_[1] * I_ff_1 * exp(b_0_[1] * I_ff_1 * I_ff_1) * local_f0f0_[i]
			+ 2.0 * a_0_[2] * I_ss_1 * exp(b_0_[2] * I_ss_1 * I_ss_1) * local_s0s0_[i]
			+ a_0_[3] * I_fs * exp(b_0_[3] * I_fs * I_fs) * local_f0s0_[i];

		return F * sigmaPK2;
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::writeToXmlForReloadMaterialProperty(std::string &filefullpath)
	{
		std::cout << "\n material properties writing " << std::endl;
		const SimTK::String xml_name("material_xml"), ele_name("material");
		XmlEngine* reload_xml = new XmlEngine(xml_name, ele_name);

		for (size_t i = 0; i != local_f0_.size(); ++i)
		{
			reload_xml->CreatXmlElement("muscle");
			reload_xml->AddAttributeToElement<Vecd>("Fibre", local_f0_[i]);
			reload_xml->AddAttributeToElement<Vecd>("Sheet", local_s0_[i]);
			reload_xml->AddElementToXmlDoc();
		}
		reload_xml->WriteToXmlFile(filefullpath);
		std::cout << "\n material properties writing finished " << std::endl;
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::readFromXmlForMaterialProperty(std::string &filefullpath)
	{
		size_t number_of_element = 0;
		XmlEngine* read_xml = new XmlEngine();
		read_xml->LoadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite_ = read_xml->root_element_.element_begin();
		for (; ele_ite_ != read_xml->root_element_.element_end(); ++ele_ite_)
		{
			Vecd fibre = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Fibre");
			local_f0_.push_back(fibre);
			Vecd sheet = read_xml->GetRequiredAttributeValue<Vecd>(ele_ite_, "Sheet");
			local_s0_.push_back(sheet);
			number_of_element++;
		}

		if(number_of_element != local_f0_.size())
		{
			std::cout << "\n Error: reload material properties does not matrch" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}else
		{
			std::cout << "\n Material properties writing finished and " << number_of_element << " number of propertyies have been readed !" << std::endl;
		}
		

		for(size_t i = 0; i < number_of_element; i++)
 		{
 			local_f0f0_.push_back(SimTK::outer(local_f0_[i], local_f0_[i]));
 			local_s0s0_.push_back(SimTK::outer(local_s0_[i], local_s0_[i]));
 			local_f0s0_.push_back(SimTK::outer(local_f0_[i], local_s0_[i]));
 		}
	}
	//=================================================================================================//
	Matd ActiveMuscle::ConstitutiveRelation(Matd& F, size_t i)
	{
		Matd passive_stress = Muscle::ConstitutiveRelation(F, i);
		return passive_stress + base_particles_->accessAParticleDataTypeReal(i) * F * f0f0_;
	}
	//=================================================================================================//
	Matd ActiveLocallyOrthotropicMuscle::ConstitutiveRelation(Matd& F, size_t i)
	{
		Matd passive_stress = LocallyOrthotropicMuscle::ConstitutiveRelation(F, i);
		return passive_stress + base_particles_->accessAParticleDataTypeReal(i) * F * local_f0f0_[i];
	}
	//=================================================================================================//
}
