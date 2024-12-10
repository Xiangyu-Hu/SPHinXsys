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
 * @file 	general_continuum.h
 * @brief 	Describe the linear elastic, J2 plasticity, and Drucker-Prager's plastic model
 * @author	Shuaihao Zhang and Xiangyu Hu
 */

#ifndef GENERAL_CONTINUUM_H
#define GENERAL_CONTINUUM_H

#include "weakly_compressible_fluid.h"

namespace SPH
{
class GeneralContinuum : public WeaklyCompressibleFluid
{
  protected:
    Real E_;                 /* Youngs or tensile modules  */
    Real G_;                 /* shear modules  */
    Real K_;                 /* bulk modules  */
    Real nu_;                /* Poisson ratio  */
    Real contact_stiffness_; /* contact-force stiffness related to bulk modulus*/
  public:
    explicit GeneralContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio)
        : WeaklyCompressibleFluid(rho0, c0), E_(0.0), G_(0.0), K_(0.0), nu_(0.0), contact_stiffness_(rho0_ * c0 * c0)
    {
        material_type_name_ = "GeneralContinuum";
        E_ = youngs_modulus;
        nu_ = poisson_ratio;
        G_ = getShearModulus(youngs_modulus, poisson_ratio);
        K_ = getBulkModulus(youngs_modulus, poisson_ratio);
        lambda0_ = getLambda(youngs_modulus, poisson_ratio);
    };
    virtual ~GeneralContinuum(){};

    Real lambda0_; /* first Lame parameter */
    Real getYoungsModulus() { return E_; };
    Real getPoissonRatio() { return nu_; };
    Real getDensity() { return rho0_; };
    Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
    Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
    Real getLambda(Real youngs_modulus, Real poisson_ratio);

    Real ContactStiffness() { return contact_stiffness_; };

    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress);

    virtual GeneralContinuum *ThisObjectPtr() override { return this; };

     //add
    class GeneralContinuumKernel:WeaklyCompressibleFluid::EosKernel
    {
      public:
        GeneralContinuumKernel(GeneralContinuum &encloser): WeaklyCompressibleFluid::EosKernel(encloser),
        E_(encloser.E_), G_(encloser.G_),K_(encloser.K_),nu_(encloser.nu_),contact_stiffness_(encloser.contact_stiffness_){};

        Real lambda0_; /* first Lame parameter */
        Real getYoungsModulus() { return E_; };
        Real getPoissonRatio() { return nu_; };
        //Real getDensity() { return rho0_; };

        Real getBulkModulus(Real youngs_modulus, Real poisson_ratio)
        {
            return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
        };

        Real getShearModulus(Real youngs_modulus, Real poisson_ratio)
        {
            return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
        };

        Real getLambda(Real youngs_modulus, Real poisson_ratio)
        {
            return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
        };

      protected:
        Real E_;                 /* Youngs or tensile modules  */
        Real G_;                 /* shear modules  */
        Real K_;                 /* bulk modules  */
        Real nu_;                /* Poisson ratio  */
        Real contact_stiffness_; /* contact-force stiffness related to bulk modulus*/
    };
};

class PlasticContinuum : public GeneralContinuum
{
  protected:
    Real c_;                            /* cohesion  */
    Real phi_;                          /* friction angle  */
    Real psi_;                          /* dilatancy angle  */
    Real alpha_phi_;                    /* Drucker-Prager's constants */
    Real k_c_;                          /* Drucker-Prager's constants */
    const Real stress_dimension_ = 3.0; /* plain strain condition */
  public:
    explicit PlasticContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real friction_angle, Real cohesion = 0, Real dilatancy = 0)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
          c_(cohesion), phi_(friction_angle), psi_(dilatancy), alpha_phi_(0.0), k_c_(0.0)
    {
        material_type_name_ = "PlasticContinuum";
        alpha_phi_ = getDPConstantsA(friction_angle);
        k_c_ = getDPConstantsK(cohesion, friction_angle);
    };
    virtual ~PlasticContinuum(){};

    Real getDPConstantsA(Real friction_angle);
    Real getDPConstantsK(Real cohesion, Real friction_angle);
    Real getFrictionAngle() { return phi_; };

    virtual Mat3d ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor);
    virtual Mat3d ReturnMapping(Mat3d &stress_tensor);

    virtual GeneralContinuum *ThisObjectPtr() override { return this; };

    //add
    class PlasticKernel:GeneralContinuum::GeneralContinuumKernel
    {
      public:

        PlasticKernel(PlasticContinuum &encloser) : GeneralContinuum::GeneralContinuumKernel(encloser),
        c_(encloser.c_), phi_(encloser.phi_),psi_(encloser.psi_),alpha_phi_(encloser.alpha_phi_),k_c_(encloser.k_c_)
        {};

        Real getDPConstantsA(Real friction_angle)
        {
          return tan(friction_angle) / sqrt(9.0 + 12.0 * tan(friction_angle) * tan(friction_angle));
        };
        
        Mat3d ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor)
        {
          Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
          Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
          Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
          Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() + stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
          Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
          Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
          Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor.trace() - k_c_;
          Real lambda_dot_ = 0;
          Mat3d g = Mat3d::Zero();
          if (f >= TinyReal)
          {
              Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
              // non-associate flow rule
              lambda_dot_ = (3.0 * alpha_phi_ * K_ * strain_rate.trace() + (G_ / sqrt(stress_tensor_J2)) * deviatoric_stress_times_strain_rate) / (9.0 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
              g = lambda_dot_ * (3.0 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() + G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
          }
          Mat3d stress_rate_temp = stress_rate_elastic - g;
          return stress_rate_temp;
        };
  
        Mat3d ReturnMapping(Mat3d &stress_tensor)
        {
          Real stress_tensor_I1 = stress_tensor.trace();
          if (-alpha_phi_ * stress_tensor_I1 + k_c_ < 0)
            stress_tensor -= (1.0 / stress_dimension_) * (stress_tensor_I1 - k_c_ / alpha_phi_) * Mat3d::Identity();
          stress_tensor_I1 = stress_tensor.trace();
          Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
          Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
          if (-alpha_phi_ * stress_tensor_I1 + k_c_ < sqrt(stress_tensor_J2))
          {
              Real r = (-alpha_phi_ * stress_tensor_I1 + k_c_) / (sqrt(stress_tensor_J2) + TinyReal);
              stress_tensor = r * deviatoric_stress_tensor + (1.0 / stress_dimension_) * stress_tensor_I1 * Mat3d::Identity();
          }
          return stress_tensor;
        };


      protected:
          Real c_;                            /* cohesion  */
          Real phi_;                          /* friction angle  */
          Real psi_;                          /* dilatancy angle  */
          Real alpha_phi_;                    /* Drucker-Prager's constants */
          Real k_c_;                          /* Drucker-Prager's constants */
          const Real stress_dimension_ = 3.0; /* plain strain condition */


    };
};

class J2Plasticity : public GeneralContinuum
{
  protected:
    Real yield_stress_;
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);

  public:
    explicit J2Plasticity(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real hardening_modulus = 0.0)
        : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
          yield_stress_(yield_stress), hardening_modulus_(hardening_modulus)
    {
        material_type_name_ = "J2Plasticity";
    };
    virtual ~J2Plasticity(){};

    Real YieldStress() { return yield_stress_; };
    Real HardeningModulus() { return hardening_modulus_; };

    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress, Real &hardening_factor);
    virtual Matd ReturnMappingShearStress(Matd &shear_stress,  Real &hardening_factor);
    virtual Real ScalePenaltyForce(Matd &shear_stress, Real &hardening_factor);
    virtual Real HardeningFactorRate(const Matd &shear_stress, Real &hardening_factor);
    virtual J2Plasticity *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // GENERAL_CONTINUUM_H
