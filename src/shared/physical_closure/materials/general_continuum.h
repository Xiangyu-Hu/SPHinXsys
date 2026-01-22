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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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
    explicit GeneralContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio);
    virtual ~GeneralContinuum() {};
    Real getYoungsModulus() { return E_; };
    Real getPoissonRatio() { return nu_; };
    Real getDensity() { return rho0_; };
    Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
    Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
    Real getLambda(Real youngs_modulus, Real poisson_ratio);
    Real ContactStiffness() { return contact_stiffness_; };
    virtual Matd ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress);

    class ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, GeneralContinuum &encloser);
        inline Real getYoungsModulus() { return E_; };
        inline Real getPoissonRatio() { return nu_; };
        inline Real getDensity() { return rho0_; };
        inline Real getBulkModulus(Real youngs_modulus, Real poisson_ratio);
        inline Real getShearModulus(Real youngs_modulus, Real poisson_ratio);
        inline Real getLambda(Real youngs_modulus, Real poisson_ratio);

      protected:
        Real E_;                 /* Youngs or tensile modules  */
        Real G_;                 /* shear modules  */
        Real K_;                 /* bulk modules  */
        Real nu_;                /* Poisson ratio  */
        Real contact_stiffness_; /* contact-force stiffness related to bulk modulus*/
        Real rho0_;              /* contact-force stiffness related to bulk modulus*/
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
    explicit PlasticContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio,
                              Real friction_angle, Real cohesion = 0, Real dilatancy = 0);
    virtual ~PlasticContinuum() {};
    Real getDPConstantsA(Real friction_angle);
    Real getDPConstantsK(Real cohesion, Real friction_angle);
    Real getFrictionAngle() { return phi_; };
    virtual Mat3d ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor);
    virtual Mat3d ReturnMapping(Mat3d &stress_tensor);

    class ConstituteKernel : public GeneralContinuum::ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, PlasticContinuum &encloser);
        inline Real getDPConstantsA(Real friction_angle);
        inline Real getFrictionAngle() { return phi_; };
        inline Mat3d StressTensorRate(UnsignedInt index_i, const Mat3d &velocity_gradient, const Mat3d &stress_tensor);
        inline Mat3d updateStressTensor(UnsignedInt index_i, const Mat3d &prev_stress_tensor, const Mat3d &stress_tensor_increment);

      protected:
        Real c_;                                                   /* cohesion  */
        Real phi_;                                                 /* friction angle  */
        Real psi_;                                                 /* dilatancy angle  */
        Real alpha_phi_;                                           /* Drucker-Prager's constants */
        Real k_c_;                                                 /* Drucker-Prager's constants */
        Real stress_dimension_ = 3.0; /* plain strain condition */ // Temporarily cancel const --need to check

        inline Mat3d ReturnMapping(Mat3d try_stress_tensor);
    };
};

class J2Plasticity : public GeneralContinuum
{
  protected:
    Real yield_stress_;
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    DiscreteVariable<Real> *dv_hardening_factor_;

  public:
    explicit J2Plasticity(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio,
                          Real yield_stress, Real hardening_modulus = 0.0);
    virtual ~J2Plasticity() {};
    Real YieldStress() { return yield_stress_; };
    Real HardeningModulus() { return hardening_modulus_; };
    Matd ConstitutiveRelationShearStressWithHardening(Matd &velocity_gradient, Matd &shear_stress, Real &hardening_factor);
    virtual Matd ReturnMappingShearStress(Matd &shear_stress, Real &hardening_factor);
    virtual Real ScalePenaltyForce(Matd &shear_stress, Real &hardening_factor);
    virtual Real HardeningFactorRate(const Matd &shear_stress, Real &hardening_factor);
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
};
} // namespace SPH
#endif // GENERAL_CONTINUUM_H