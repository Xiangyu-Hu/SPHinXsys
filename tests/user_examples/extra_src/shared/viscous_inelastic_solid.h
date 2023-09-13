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
 * @file 	viscous_inelastic_solid.h 
 * @brief 	These is a class for defining properties of viscous solid materials,
 * 			 based on isotropic linear elastic solid.
 * @author	Xiaojing Tang and  Xiangyu Hu 
 */
#pragma once
#include "inelastic_solid.h"

namespace SPH
{
/**
 * @class ViscousPlasticSolid
 * @brief Class for plastic solid with viscosity
 */
class ViscousPlasticSolid : public PlasticSolid
{
  protected:
    Real viscous_modulus_;
    Real Herschel_Bulkley_power_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    StdLargeVec<Matd> inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */

  public:
    /** Constructor */
    explicit ViscousPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real viscous_modulus, Real herschel_bulkley_power)    
        : PlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress), viscous_modulus_(viscous_modulus), Herschel_Bulkley_power_(herschel_bulkley_power) 
    {
        material_type_name_ = "ViscousPlasticSolid";
    };
    virtual ~ViscousPlasticSolid(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override
	{
		PlasticSolid::initializeLocalParameters(base_particles);
		base_particles->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain",
										[&](size_t i) -> Matd
										{ return Matd::Identity(); });
		base_particles->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
		
	};
    Real ViscousModulus() { return viscous_modulus_; };
    /** compute the stress through deformation, and plastic relaxation. */
    virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) 
	{
		Matd F = deformation;
		Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
		Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
		Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
		Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
		Real deviatoric_PK_norm = deviatoric_PK.norm();
		Real trial_function = deviatoric_PK_norm - sqrt_2_over_3_ * yield_stress_; 
		if (trial_function > 0.0)
   		{
			Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
			Real deviatoric_PK_norm_Mid = 0.0;
			Real deviatoric_PK_norm_Max = deviatoric_PK_norm;
			Real deviatoric_PK_norm_Min = sqrt_2_over_3_ * yield_stress_;
			Real predicted_func = 0.0;
			Real Precision = 1.0e-6;
			Real Relative_Error;
			do
			{
				deviatoric_PK_norm_Mid = (deviatoric_PK_norm_Max + deviatoric_PK_norm_Min) / 2.0;
				predicted_func = pow(viscous_modulus_, 1.0 / Herschel_Bulkley_power_) * (deviatoric_PK_norm_Mid - deviatoric_PK_norm) + 
						2.0 * renormalized_shear_modulus * dt * pow((deviatoric_PK_norm_Mid - sqrt_2_over_3_ * yield_stress_), 1.0 / Herschel_Bulkley_power_);
				if (predicted_func < 0.0)
				{
					deviatoric_PK_norm_Min = deviatoric_PK_norm_Mid;
				}
				else
				{
					deviatoric_PK_norm_Max = deviatoric_PK_norm_Mid;
				}
				Relative_Error = predicted_func / deviatoric_PK_norm;
			} while (fabs(Relative_Error) >= Precision); 

			deviatoric_PK = deviatoric_PK_norm_Mid * deviatoric_PK / deviatoric_PK_norm;
			Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();
			normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
		}

		Matd inverse_F = F.inverse();
		Matd inverse_F_T = inverse_F.transpose();
		be = pow(F.determinant(), (2.0 / 3.0)) * normalized_be;
		inverse_plastic_strain_[index_i] = inverse_F * be * inverse_F_T;

		return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
	 
	};

    virtual ViscousPlasticSolid *ThisObjectPtr() override { return this; };
};
 
} // namespace SPH
