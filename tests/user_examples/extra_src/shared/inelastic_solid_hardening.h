/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	inelastic_solid_hardening.h 
 * @brief 	These are classes for define properties of elastic solid materials.
 *			These classes are based on isotropic linear elastic solid.
 * 			Several more complex materials, including neo-hookean, FENE noe-hookean
 *			and anisotropic muscle, are derived from the basic elastic solid class.
 * @author	Xiangyu Hu and Chi Zhang
 */
#pragma once

#include "inelastic_solid.h"

namespace SPH
{

/**
 * @class NonlinearHardeningPlasticSolid
 * @brief Class for plastic solid with nonlinear hardening
 */
class NonLinearHardeningPlasticSolid : public HardeningPlasticSolid
{
protected:
	Real saturation_flow_stress_, saturation_exponent_;
public:
  /** Constructor */
  explicit NonLinearHardeningPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
    Real hardening_modulus, Real saturation_flow_stress, Real saturation_exponent)
    : HardeningPlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress, hardening_modulus),
    saturation_flow_stress_(saturation_flow_stress), saturation_exponent_(saturation_exponent)
  {
    material_type_name_ = "NonLinearHardeningPlasticSolid";

  };
  virtual ~NonLinearHardeningPlasticSolid() {};

 Real NonlinearHardening(Real hardening_parameter_pre)
	{
	  return  (hardening_modulus_ * hardening_parameter_pre + yield_stress_
			       + (saturation_flow_stress_ - yield_stress_)* (1 - exp(-saturation_exponent_ * hardening_parameter_pre)));
	};

	Real NonlinearHardeningDerivative(Real hardening_parameter_pre)
	{
		return  (hardening_modulus_+ saturation_exponent_ * (saturation_flow_stress_ - yield_stress_)
			       * exp(-saturation_exponent_ * hardening_parameter_pre));
	};

	/** compute the stress through defoemation, and plastic relaxation. */
	virtual Matd PlasticConstitutiveRelation(const Matd  &F,  size_t index_i, Real dt = 0.0)
	{
 
  		Matd normalized_F =  F * pow(F.determinant(), -OneOverDimensions);
		Matd normalized_be = normalized_F * inverse_plastic_strain_[index_i] * normalized_F.transpose();
		Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
		Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
		Real deviatoric_PK_norm = deviatoric_PK.norm();

		Real relax_increment = 0.0;
		Real trial_function =  deviatoric_PK_norm - sqrt_2_over_3_ * NonlinearHardening(hardening_parameter_[index_i]);
		if (trial_function > 0.0)
		{
			Real renormalized_shear_modulus = normalized_be_isentropic * G0_;        
			while (trial_function > 0.0)
			{		
				Real function_relax_increment_derivative = -2.0 * renormalized_shear_modulus
					* (1.0 + NonlinearHardeningDerivative(hardening_parameter_[index_i] + sqrt_2_over_3_ * relax_increment) / 3.0 / renormalized_shear_modulus);
				relax_increment -=  trial_function / function_relax_increment_derivative;
        
				trial_function = deviatoric_PK_norm
					- sqrt_2_over_3_ * NonlinearHardening(hardening_parameter_[index_i] + sqrt_2_over_3_ * relax_increment)
					- 2.0 * renormalized_shear_modulus * relax_increment;	
			}      
			hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
			deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
			normalized_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();	
    
		}
 
		Matd inverse_normalized_F = normalized_F.inverse();
		Matd inverse_normalized_F_T = inverse_normalized_F.transpose();;
		inverse_plastic_strain_[index_i] = inverse_normalized_F * normalized_be * inverse_normalized_F_T;
		
		return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_normalized_F_T;
	
	};

	virtual NonLinearHardeningPlasticSolid *ThisObjectPtr() override { return this; };
};

} // namespace SPH
