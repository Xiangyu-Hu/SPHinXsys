/**
 * @file elastic_solid.cpp
 * @author Chi Zhang and Xiangyu Hu
 */

#include "inelastic_solid.h"

namespace SPH {
	//=================================================================================================//
	void HardeningPlasticSolid::initializePlasticParameters()
	{
		base_particles_->registerAVariable<indexMatrix, Matd>(inverse_plastic_strain_, "InversePlasticRightCauchyStrain", Matd(1.0));
		base_particles_->registerAVariable<indexScalar, Real>(hardening_parameter_, "HardeningParameter");
		base_particles_->addAVariableToRestart<indexMatrix, Matd>("InversePlasticRightCauchyStrain");
		base_particles_->addAVariableToRestart<indexScalar, Real>("HardeningParameter");
	}
	//=================================================================================================//
	void HardeningPlasticSolid::assignElasticSolidParticles(ElasticSolidParticles* elastic_particles)
	{
		ElasticSolid::assignElasticSolidParticles(elastic_particles);
		initializePlasticParameters();
	}
	//=================================================================================================//
	Matd HardeningPlasticSolid::PlasticConstitutiveRelation(const Matd& F, size_t index_i, Real dt)
	{
		Matd be = F * inverse_plastic_strain_[index_i] * (~F);
		Matd normalized_be = be * pow(SimTK::det(be), -one_over_dimensions_);
		Real normalized_be_isentropic = normalized_be.trace() * one_over_dimensions_;
		Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd(1.0));
		Real deviatoric_PK_norm = deviatoric_PK.norm();
		Real trial_function =
			deviatoric_PK_norm - sqrt_2_over_3_ * (hardening_modulus_ * hardening_parameter_[index_i] + yield_stress_);
		if (trial_function > 0.0)
		{
			Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
			Real relax_increment = 0.5 * trial_function / (renormalized_shear_modulus + hardening_modulus_ / 3.0);
			hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
			deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
			Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic;
			normalized_be = relaxed_be * pow(det(relaxed_be), -one_over_dimensions_);
		}
		Matd inverse_F = SimTK::inverse(F);
		Matd inverse_F_T = ~inverse_F;
		inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

		return (deviatoric_PK + VolumetricKirchhoff(SimTK::det(F)) * Matd(1.0)) * inverse_F_T;
	}
	//=================================================================================================//
}
