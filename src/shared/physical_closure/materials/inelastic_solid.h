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
 * @file 	inelastic_solid.h
 * @brief 	These are classes for define properties of plastic solid materials.
 * 			Several plastic materials including linear hardening, non-linear hardening and viscous plastic, are presented.
 * @author	Xiaojing Tang, Dong Wu, and Xiangyu Hu
 */
#ifndef INELASTIC_SOLID_H
#define INELASTIC_SOLID_H

#include "elastic_solid.h"

namespace SPH
{
/**
 * @class PlasticSolid
 * @brief Abstract class for a generalized plastic solid
 */
class PlasticSolid : public NeoHookeanSolid
{
  protected:
    Real yield_stress_;

  public:
    /** Constructor */
    explicit PlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress)
        : NeoHookeanSolid(rho0, youngs_modulus, poisson_ratio), yield_stress_(yield_stress)
    {
        material_type_name_ = "PlasticSolid";
    };
    virtual ~PlasticSolid() {};

    Real YieldStress() { return yield_stress_; };
    /** compute the elastic part of normalized left Cauchy-Green deformation gradient tensor. */
    virtual Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0) = 0;
};

/**
 * @class HardeningPlasticSolid
 * @brief Class for plastic solid with hardening
 */
class HardeningPlasticSolid : public PlasticSolid
{
    DiscreteVariable<Matd> *dv_inverse_plastic_strain_;
    DiscreteVariable<Real> *dv_hardening_parameter_;

  protected:
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    Matd *inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */
    Real *hardening_parameter_;    /**< hardening parameter */

  public:
    explicit HardeningPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                   Real yield_stress, Real hardening_modulus);
    virtual ~HardeningPlasticSolid() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    Real HardeningModulus() { return hardening_modulus_; };
    /** compute the elastic part of normalized left Cauchy-Green deformation gradient tensor. */
    virtual Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0) override;

    class ConstituteKernel : public PlasticSolid::ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, HardeningPlasticSolid &encloser);
        Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0);

      protected:
        Real hardening_modulus_, yield_stress_;
        Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
        Matd *inverse_plastic_strain_;
        Real *hardening_parameter_;
    };
};

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
        return (hardening_modulus_ * hardening_parameter_pre + yield_stress_ + (saturation_flow_stress_ - yield_stress_) * (1 - exp(-saturation_exponent_ * hardening_parameter_pre)));
    };

    Real NonlinearHardeningDerivative(Real hardening_parameter_pre)
    {
        return (hardening_modulus_ + saturation_exponent_ * (saturation_flow_stress_ - yield_stress_) * exp(-saturation_exponent_ * hardening_parameter_pre));
    };
    /** compute the elastic part of normalized left Cauchy-Green deformation gradient tensor. */
    virtual Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0) override;
};

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
    Matd *inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */

  public:
    /** Constructor */
    explicit ViscousPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real viscous_modulus, Real herschel_bulkley_power)
        : PlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress), viscous_modulus_(viscous_modulus), Herschel_Bulkley_power_(herschel_bulkley_power)
    {
        material_type_name_ = "ViscousPlasticSolid";
    };
    virtual ~ViscousPlasticSolid() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    Real ViscousModulus() { return viscous_modulus_; };
    /** compute the elastic part of normalized left Cauchy-Green deformation gradient tensor. */
    virtual Matd ElasticLeftCauchy(const Matd &deformation, size_t index_i, Real dt = 0.0) override;
};
} // namespace SPH
#endif // INELASTIC_SOLID_H
