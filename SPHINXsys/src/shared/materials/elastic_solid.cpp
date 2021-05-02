/** 
 * @file elastic_solid.cpp
 * @author Chi Zhang and Xiangyu Hu
 */

#include "elastic_solid.h"

#include "base_body.h"
#include "solid_particles.h"

namespace SPH {
	//=================================================================================================//
	void ElasticSolid::assignElasticSolidParticles(ElasticSolidParticles* elastic_particles) 
	{
		elastic_particles_ = elastic_particles;
	}
	//=================================================================================================//
	Real ElasticSolid::ViscousTimeStepSize(Real smoothing_length)
	{
		Real total_viscosity = eta_0_;
		return 0.5 * smoothing_length * smoothing_length * rho_0_ / (total_viscosity + TinyReal);
	}
	//=================================================================================================//
	Matd ElasticSolid::NumericalDampingStress(Matd& F, Matd& dF_dt, Real numerical_viscosity, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (~dF_dt * F + ~F * dF_dt);
		Matd sigmaPK2 = numerical_viscosity * strain_rate;
		return sigmaPK2;
	}
	//=================================================================================================//
	Real ElasticSolid::NumericalViscosity(Real smoothing_length)
	{
		return 0.5 * rho_0_ * c_0_ * smoothing_length;
	}
	//=================================================================================================//
	void LinearElasticSolid::assignDerivedMaterialParameters() 
	{
		Solid::assignDerivedMaterialParameters();
		setLambda();
		setShearModulus();
		setSoundSpeed();
		setContactStiffness();
		std::cout << "The speed of sound: " << c_0_ << std::endl;
		std::cout << "The Lambda: " << lambda_0_ << std::endl;
		std::cout << "Contact stiffness: " << contact_stiffness_ << std::endl;
	};
	//=================================================================================================//
	Matd LinearElasticSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = lambda_0_ * strain.trace() * Matd(1.0)
			+ 2.0 * G_0_ * strain;
		return sigmaPK2;
	}
	//=================================================================================================//
	void LinearElasticSolid::setSoundSpeed()
	{
		c_0_ = sqrt(E_0_ / 3.0 / (1.0 - 2.0 * nu_) / rho_0_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setShearModulus()
	{
		G_0_ =  0.5 * E_0_ / (1.0 + nu_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setLambda()
	{
		lambda_0_ = nu_ * E_0_ / (1.0 + nu_) / (1.0 - 2.0 * nu_);
	}
	//=================================================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = G_0_ * Matd(1.0)
			+ (lambda_0_ * log(det(F)) - G_0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = G_0_ / (1.0 - 2.0 * strain.trace() / j1_m_) * Matd(1.0)
			+ (lambda_0_ * log(det(F)) - G_0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	void Muscle::setSoundSpeed()
	{
		c_0_ = sqrt(bulk_modulus_ / rho_0_);
	}
	//=================================================================================================//
	void Muscle::setLambda()
	{
		Real shear_modulus_ref = a_0_[0]* b_0_[0] + 2.0 * a_0_[1] * b_0_[1] 
						 + 2.0 * a_0_[2] * b_0_[2] + a_0_[3] * b_0_[3];
		lambda_0_ = bulk_modulus_ - 2.0 * shear_modulus_ref / 3.0;
	}
	//=================================================================================================//
	void Muscle::assignDerivedMaterialParameters()
	{
		f0f0_ = SimTK::outer(f0_, f0_);
		f0s0_ = SimTK::outer(f0_, s0_);
		s0s0_ = SimTK::outer(s0_, s0_);
		setSoundSpeed();
		setLambda();
		setContactStiffness();
		std::cout << "The speed of sound: " << c_0_ << std::endl;
		std::cout << "The Lambda: " << lambda_0_ << std::endl;
		std::cout << "Contact stiffness: " << contact_stiffness_ << std::endl;
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

		return sigmaPK2;
	}
	//=================================================================================================//
	Real Muscle::YoungsModulus()
	{ 
		return 9.0 * bulk_modulus_ * a_0_[0] / (3.0 * bulk_modulus_ + a_0_[0]); 
	};
	//=================================================================================================//
	Real Muscle::PoissonRatio()
	{ 
		return (3.0 * bulk_modulus_ - 2.0 * a_0_[0]) / (3.0 * bulk_modulus_ + a_0_[0]) / 2.0; 
	};
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

		return sigmaPK2;
	}		
	//=================================================================================================//
	void LocallyOrthotropicMuscle::assignElasticSolidParticles(ElasticSolidParticles* elastic_particles) 
	{
		Muscle::assignElasticSolidParticles(elastic_particles);
		initializeFiberAndSheet();
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::initializeFiberAndSheet()
	{
		size_t total_real_particles = base_particles_->total_real_particles_;
		registerAVariableToParticleData<indexVector, Vecd>(parameter_data_, parameter_maps_, local_f0_, "Fiber");
		registerAVariableToParticleData<indexVector, Vecd>(parameter_data_, parameter_maps_, local_s0_, "Sheet");
		local_f0_.resize(total_real_particles, Vecd(0));
		local_s0_.resize(total_real_particles, Vecd(0));
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::readFromXmlForLocalParameters(std::string &filefullpath)
	{
 		BaseMaterial::readFromXmlForLocalParameters(filefullpath);
		size_t total_real_particles = base_particles_->total_real_particles_;
		for(size_t i = 0; i != total_real_particles; i++)
 		{
 			local_f0f0_.push_back(SimTK::outer(local_f0_[i], local_f0_[i]));
 			local_s0s0_.push_back(SimTK::outer(local_s0_[i], local_s0_[i]));
 			local_f0s0_.push_back(SimTK::outer(local_f0_[i], local_s0_[i]));
 		}
	}
	//=================================================================================================//
}
