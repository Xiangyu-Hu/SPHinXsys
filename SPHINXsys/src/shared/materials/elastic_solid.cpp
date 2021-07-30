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
	void ElasticSolid::assignDerivedMaterialParameters()
	{
		Solid::assignDerivedMaterialParameters();
		setReferenceSoundSpeed();
		setTensileWaveSpeed();
		setShearWaveSpeed();
		setYoungsModulus();
		setShearModulus();
		setBulkModulus();
		setPoissonRatio();
		setContactStiffness();
	}
	//=================================================================================================//
	Matd ElasticSolid::NumericalDampingRightCauchy(Matd& F, Matd& dF_dt, Real smoothing_length, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (~dF_dt * F + ~F * dF_dt);
		Matd normal_rate = getDiagonal(strain_rate);
		return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * smoothing_length;
	}
	//=================================================================================================//
	Matd ElasticSolid::NumericalDampingLeftCauchy(Matd& F, Matd& dF_dt, Real smoothing_length, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (dF_dt * ~F + F * ~dF_dt);
		Matd normal_rate = getDiagonal(strain_rate);
		return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * smoothing_length;
	}
	//=================================================================================================//
	Real ElasticSolid::NumericalDamping(Real dE_dt_ij, Real smoothing_length)
	{
		return 0.5 * rho0_ * c0_ * dE_dt_ij * smoothing_length;
	}
	//=================================================================================================//
	Real ElasticSolid::NumericalViscosity(Real smoothing_length)
	{
		return 0.5 * rho0_ * c0_ * smoothing_length;
	}
	//=================================================================================================//
	Matd ElasticSolid::DeviatoricKirchhoff(const Matd& deviatoric_be)
	{
		return G0_ * deviatoric_be;
	}
	//=================================================================================================//
	Real LinearElasticSolid::getBulkModulus()
	{
		return youngs_modulus_ / 3.0 / (1.0 - 2.0 * poisson_ratio_);
	}
	//=================================================================================================//
	Real LinearElasticSolid::getShearModulus()
	{
		return 0.5 * youngs_modulus_ / (1.0 + poisson_ratio_);
	}
	//=================================================================================================//
	Real LinearElasticSolid::getLambda()
	{
		return nu_ * youngs_modulus_ / (1.0 + poisson_ratio_) / (1.0 - 2.0 * poisson_ratio_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setReferenceSoundSpeed()
	{
		c0_ = sqrt(getBulkModulus() / rho0_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setTensileWaveSpeed()
	{
		ct0_ = sqrt(youngs_modulus_ / rho0_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setShearWaveSpeed()
	{
		cs0_ = sqrt(getShearModulus() / rho0_);
	}
	//=================================================================================================//
	void LinearElasticSolid::setShearModulus()
	{
		G0_ = getShearModulus();
	}
	//=================================================================================================//
	void LinearElasticSolid::setBulkModulus()
	{
		K0_ = getBulkModulus();
	}
	//=================================================================================================//
	void LinearElasticSolid::assignDerivedMaterialParameters()
	{
		ElasticSolid::assignDerivedMaterialParameters();
		lambda0_ = getLambda();
		std::cout << "The speed of sound: " << c0_ << std::endl;
		std::cout << "The Lambda: " << lambda0_ << std::endl;
		std::cout << "Contact stiffness: " << contact_stiffness_ << std::endl;
	};
	//=================================================================================================//
	Matd LinearElasticSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = lambda0_ * strain.trace() * Matd(1.0) + 2.0 * G0_ * strain;
		return sigmaPK2;
	}
	//=================================================================================================//
	Real  LinearElasticSolid::VolumetricKirchhoff(Real J)
	{
		return  K0_ * J * (J - 1);
	}
	//=================================================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = G0_ * Matd(1.0) + (lambda0_ * log(det(F)) - G0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	Real NeoHookeanSolid::VolumetricKirchhoff(Real J)
	{
		return  0.5 * K0_ * (J * J - 1);
	}
	//=================================================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd& F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = G0_ / (1.0 - 2.0 * strain.trace() / j1_m_) * Matd(1.0)
			+ (lambda0_ * log(det(F)) - G0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	Real Muscle::getShearModulus()
	{
		return a0_[0] * b0_[0] + 2.0 * a0_[1] * b0_[1] + 2.0 * a0_[2] * b0_[2] + a0_[3] * b0_[3];
	}
	//=================================================================================================//
	Real Muscle::getPoissonRatio()
	{
		return 0.5 * (3.0 * bulk_modulus_ - 2.0 * getShearModulus()) / (3.0 * bulk_modulus_ + getShearModulus());
	}
	//=================================================================================================//
	Real Muscle::getLambda()
	{
		return bulk_modulus_ - 2.0 * getShearModulus() / 3.0;
	}
	//=================================================================================================//
	Real Muscle::getYoungsModulus()
	{
		return 3.0 * bulk_modulus_ * (1.0 - 2.0 * getPoissonRatio());
	}
	//=================================================================================================//
	void Muscle::setReferenceSoundSpeed()
	{
		c0_ = sqrt(bulk_modulus_ / rho0_);
	}
	//=================================================================================================//
	void Muscle::setTensileWaveSpeed()
	{
		ct0_ = sqrt(getYoungsModulus() / rho0_);
	}
	//=================================================================================================//
	void Muscle::setShearWaveSpeed()
	{
		cs0_ = sqrt(getShearModulus() / rho0_);
	}
	//=================================================================================================//
	void Muscle::setYoungsModulus()
	{
		E0_ = getYoungsModulus();
	}
	//=================================================================================================//
	void Muscle::setShearModulus()
	{
		G0_ = getShearModulus();
	}
	//=================================================================================================//
	void Muscle::setPoissonRatio()
	{
		nu_ = getPoissonRatio();
	}
	//=================================================================================================//
	void Muscle::assignDerivedMaterialParameters()
	{
		ElasticSolid::assignDerivedMaterialParameters();
		lambda0_ = getLambda();
		f0f0_ = SimTK::outer(f0_, f0_);
		f0s0_ = SimTK::outer(f0_, s0_);
		s0s0_ = SimTK::outer(s0_, s0_);
		std::cout << "The speed of sound: " << c0_ << std::endl;
		std::cout << "The Lambda: " << lambda0_ << std::endl;
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
		Matd sigmaPK2 = a0_[0] * exp(b0_[0] * I_1_1) * Matd(1.0)
			+ (lambda0_ * ln_J - a0_[0]) * inverse(right_cauchy)
			+ 2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * f0f0_
			+ 2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * s0s0_
			+ a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * f0s0_;

		return sigmaPK2;
	}
	//=================================================================================================//
	Real  Muscle::VolumetricKirchhoff(Real J)
	{
		return  K0_ * J * (J - 1);
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
		Matd sigmaPK2 = a0_[0] * exp(b0_[0] * I_1_1) * Matd(1.0)
			+ (lambda0_ * ln_J - a0_[0]) * inverse(right_cauchy)
			+ 2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * local_f0f0_[i]
			+ 2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * local_s0s0_[i]
			+ a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * local_f0s0_[i];

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
		base_particles_->registerAVariable<indexVector, Vecd>(local_f0_, "Fiber");
		base_particles_->registerAVariable<indexVector, Vecd>(local_s0_, "Sheet");
		base_particles_->addAVariableNameToList<indexVector, Vecd>(reload_local_parameters_, "Fiber");
		base_particles_->addAVariableNameToList<indexVector, Vecd>(reload_local_parameters_, "Sheet");
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
