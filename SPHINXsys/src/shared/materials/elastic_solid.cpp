/** 
 * @file elastic_solid.cpp
 * @author Chi Zhang and Xiangyu Hu
 */

#include "elastic_solid.h"

#include "base_body.h"
#include "solid_particles.h"

namespace SPH
{
	//=================================================================================================//
	void ElasticSolid::setSoundSpeeds()
	{
		c0_ = sqrt(K0_ / rho0_);
		ct0_ = sqrt(E0_ / rho0_);
		cs0_ = sqrt(G0_ / rho0_);
	};
	//=================================================================================================//
	void ElasticSolid::assignElasticSolidParticles(ElasticSolidParticles *elastic_particles)
	{
		elastic_particles_ = elastic_particles;
	}
	//=================================================================================================//
	Matd ElasticSolid::
		NumericalDampingRightCauchy(Matd &F, Matd &dF_dt, Real smoothing_length, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (~dF_dt * F + ~F * dF_dt);
		Matd normal_rate = getDiagonal(strain_rate);
		return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * smoothing_length;
	}
	//=================================================================================================//
	Matd ElasticSolid::
		NumericalDampingLeftCauchy(Matd &F, Matd &dF_dt, Real smoothing_length, size_t particle_index_i)
	{
		Matd strain_rate = 0.5 * (dF_dt * ~F + F * ~dF_dt);
		Matd normal_rate = getDiagonal(strain_rate);
		return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * smoothing_length;
	}
	//=================================================================================================//
	Real ElasticSolid::PairNumericalDamping(Real dE_dt_ij, Real smoothing_length)
	{
		return 0.5 * rho0_ * c0_ * dE_dt_ij * smoothing_length;
	}
	//=================================================================================================//
	Matd ElasticSolid::DeviatoricKirchhoff(const Matd &deviatoric_be)
	{
		return G0_ * deviatoric_be;
	}
	//=================================================================================================//
	LinearElasticSolid::
		LinearElasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio) : ElasticSolid(rho0)
	{
		material_type_ = "LinearElasticSolid";
		E0_ = youngs_modulus;
		nu_ = poisson_ratio;
		G0_ = getShearModulus(youngs_modulus, poisson_ratio);
		K0_ = getBulkModulus(youngs_modulus, poisson_ratio);
		lambda0_ = getLambda(youngs_modulus, poisson_ratio);
		setSoundSpeeds();
		setContactStiffness(c0_);
	}
	//=================================================================================================//
	Real LinearElasticSolid::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
	{
		return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
	}
	//=================================================================================================//
	Real LinearElasticSolid::getShearModulus(Real youngs_modulus, Real poisson_ratio)
	{
		return 0.5 * youngs_modulus/ (1.0 + poisson_ratio);
	}
	//=================================================================================================//
	Real LinearElasticSolid::getLambda(Real youngs_modulus, Real poisson_ratio)
	{
		return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
	}
	//=================================================================================================//
	Matd LinearElasticSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd strain = 0.5 * (~F * F - Matd(1.0));
		Matd sigmaPK2 = lambda0_ * strain.trace() * Matd(1.0) + 2.0 * G0_ * strain;
		return sigmaPK2;
	}
	//=================================================================================================//
	Matd LinearElasticSolid::EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i)
	{
		return lambda0_ * almansi_strain.trace() * Matd(1.0) + 2.0 * G0_ * almansi_strain;
	}
	//=================================================================================================//
	Real LinearElasticSolid::VolumetricKirchhoff(Real J)
	{
		return K0_ * J * (J - 1);
	}
	//=================================================================================================//
	Matd NeoHookeanSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd sigmaPK2 = G0_ * Matd(1.0) + (lambda0_ * log(det(F)) - G0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	Matd NeoHookeanSolid::EulerianConstitutiveRelation(Matd &almansi_strain, Matd &F, size_t particle_index_i)
	{
		Real J = det(F);
		Matd B = inverse(-2.0 * almansi_strain + Matd(1.0));
		Matd cauchy_stress = 0.5 * K0_ * (J - 1.0 / J) * Matd(1.0)
			+ G0_ * pow(J, -2.0 / (Real)Dimensions - 1.0) * (B - B.trace() / (Real)Dimensions * Matd(1.0));
		return cauchy_stress;
	}
	//=================================================================================================//
	Real NeoHookeanSolid::VolumetricKirchhoff(Real J)
	{
		return 0.5 * K0_ * (J * J - 1);
	}
	//=================================================================================================//
	Matd FeneNeoHookeanSolid::ConstitutiveRelation(Matd &F, size_t particle_index_i)
	{
		Matd right_cauchy = ~F * F;
		Matd strain = 0.5 * (right_cauchy - Matd(1.0));
		Matd sigmaPK2 = G0_ / (1.0 - 2.0 * strain.trace() / j1_m_) * Matd(1.0) +
						(lambda0_ * log(det(F)) - G0_) * inverse(right_cauchy);
		return sigmaPK2;
	}
	//=================================================================================================//
	Real Muscle::getShearModulus(const Real (&a0)[4], const Real (&b0)[4])
	{
		// This is only the background material property.
		// The previous version seems not correct because it leads to 
		// that shear modulus is even bigger than bulk modulus.
		return a0[0];
	}
	//=================================================================================================//
	Real Muscle::getPoissonRatio(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4])
	{
		Real shear_modulus = getShearModulus(a0, b0);
		return 0.5 * (3.0 * bulk_modulus - 2.0 * shear_modulus) /
			   (3.0 * bulk_modulus + shear_modulus);
	}
	//=================================================================================================//
	Real Muscle::getYoungsModulus(Real bulk_modulus, const Real (&a0)[4], const Real (&b0)[4])
	{
		return 3.0 * bulk_modulus * (1.0 - 2.0 * getPoissonRatio(bulk_modulus, a0, b0));
	}
	//=================================================================================================//
	Matd Muscle::ConstitutiveRelation(Matd &F, size_t i)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy * f0_, f0_) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy * s0_, s0_) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy * f0_, s0_);
		Real I_1_1 = right_cauchy.trace() - Real(Dimensions);
		Matd sigmaPK2 = a0_[0] * exp(b0_[0] * I_1_1) * Matd(1.0) +
						(lambda0_ * log(det(F)) - a0_[0]) * inverse(right_cauchy) +
						2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * f0f0_ +
						2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * s0s0_ +
						a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * f0s0_;

		return sigmaPK2;
	}
	//=================================================================================================//
	Real Muscle::VolumetricKirchhoff(Real J)
	{
		return K0_ * J * (J - 1);
	}
	//=================================================================================================//
	Matd LocallyOrthotropicMuscle::ConstitutiveRelation(Matd &F, size_t i)
	{
		Matd right_cauchy = ~F * F;
		Real I_ff_1 = SimTK::dot(right_cauchy * local_f0_[i], local_f0_[i]) - 1.0;
		Real I_ss_1 = SimTK::dot(right_cauchy * local_s0_[i], local_s0_[i]) - 1.0;
		Real I_fs = SimTK::dot(right_cauchy * local_f0_[i], local_s0_[i]);
		Real I_1_1 = right_cauchy.trace() - Real(Dimensions);
		Matd sigmaPK2 = a0_[0] * exp(b0_[0] * I_1_1) * Matd(1.0) +
						(lambda0_ * log(det(F)) - a0_[0]) * inverse(right_cauchy) +
						2.0 * a0_[1] * I_ff_1 * exp(b0_[1] * I_ff_1 * I_ff_1) * local_f0f0_[i] +
						2.0 * a0_[2] * I_ss_1 * exp(b0_[2] * I_ss_1 * I_ss_1) * local_s0s0_[i] +
						a0_[3] * I_fs * exp(b0_[3] * I_fs * I_fs) * local_f0s0_[i];

		return sigmaPK2;
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::assignElasticSolidParticles(ElasticSolidParticles *elastic_particles)
	{
		Muscle::assignElasticSolidParticles(elastic_particles);
		initializeFiberAndSheet();
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::initializeFiberAndSheet()
	{
		base_particles_->registerAVariable<Vecd>(local_f0_, "Fiber");
		base_particles_->registerAVariable<Vecd>(local_s0_, "Sheet");
		base_particles_->addAVariableNameToList<Vecd>(reload_local_parameters_, "Fiber");
		base_particles_->addAVariableNameToList<Vecd>(reload_local_parameters_, "Sheet");
	}
	//=================================================================================================//
	void LocallyOrthotropicMuscle::readFromXmlForLocalParameters(const std::string &filefullpath)
	{
		BaseMaterial::readFromXmlForLocalParameters(filefullpath);
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t i = 0; i != total_real_particles; i++)
		{
			local_f0f0_.push_back(SimTK::outer(local_f0_[i], local_f0_[i]));
			local_s0s0_.push_back(SimTK::outer(local_s0_[i], local_s0_[i]));
			local_f0s0_.push_back(SimTK::outer(local_f0_[i], local_s0_[i]));
		}
	}
	//=================================================================================================//
}
