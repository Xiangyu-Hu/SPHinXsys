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
 * @file 	weakly_compressible_fluid.h
 * @brief 	Describe the weakly compressible fluid which is used
 * 			model incompressible fluids. Here, we have included several equation of states.
 * 			Futhermore, A typical non-newtonian fluid model is included.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef WEAKLY_COMPRESSIBLE_FLUID_H
#define WEAKLY_COMPRESSIBLE_FLUID_H

#include "base_material.h"

namespace SPH
{
/**
 * @class WeaklyCompressibleFluid
 * @brief Linear equation of state (EOS).
 */
class WeaklyCompressibleFluid : public Fluid
{
  protected:
    Real p0_; /**< reference pressure */
  public:
    explicit WeaklyCompressibleFluid(Real rho0, Real c0, Real mu = 0.0)
        : Fluid(rho0, c0, mu), p0_(rho0 * c0 * c0)
    {
        material_type_name_ = "WeaklyCompressibleFluid";
    };
    explicit WeaklyCompressibleFluid(ConstructArgs<Real, Real, Real> args)
        : WeaklyCompressibleFluid(std::get<0>(args), std::get<1>(args), std::get<2>(args)) {};
    virtual ~WeaklyCompressibleFluid() {};

    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
    virtual WeaklyCompressibleFluid *ThisObjectPtr() override { return this; };

    class EosKernel : Fluid::EosKernel
    {
      public:
        EosKernel(WeaklyCompressibleFluid &encloser);

        Real getPressure(Real rho)
        {
            return p0_ * (rho / rho0_ - 1.0);
        };

        Real DensityFromPressure(Real p)
        {
            return rho0_ * (p / p0_ + 1.0);
        };

        Real getSoundSpeed(Real p = 0.0, Real rho = 1.0)
        {
            return c0_;
        };

      protected:
        Real p0_;
    };
};

/**
 * @class WeaklyCompressibleFluidFreeSurface
 * @brief Equation of state (EOS) with cut-off pressure.
 */
template <class WeaklyCompressibleFluidType>
class WeaklyCompressibleFluidFreeSurface : public WeaklyCompressibleFluidType
{
  protected:
    Real cutoff_pressure_, cutoff_density_;

  public:
    template <typename... Args>
    explicit WeaklyCompressibleFluidFreeSurface(Real cutoff_pressure, Args &&...args)
        : WeaklyCompressibleFluidType(std::forward<Args>(args)...),
          cutoff_pressure_(cutoff_pressure),
          cutoff_density_(WeaklyCompressibleFluidType::DensityFromPressure(cutoff_pressure))
    {
        WeaklyCompressibleFluidType::material_type_ += "FreeSurface";
    };
    virtual ~WeaklyCompressibleFluidFreeSurface() {};

    virtual Real getPressure(Real rho) override
    {
        return rho < cutoff_density_ ? cutoff_pressure_ : WeaklyCompressibleFluid::getPressure(rho);
    };
};

/**
 * @class SymmetricTaitFluid
 * @brief Tait EOS for positive and negative pressure symmetrically.
 */
class SymmetricTaitFluid : public WeaklyCompressibleFluid
{
  protected:
    int gamma_; /**< determine the stiffness of the fluid */
  public:
    explicit SymmetricTaitFluid(Real rho0, Real c0, int gamma)
        : WeaklyCompressibleFluid(rho0, c0), gamma_(gamma)
    {
        material_type_name_ = "SymmetricTaitFluid";
    };
    virtual ~SymmetricTaitFluid() {};

    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
};

/**
 * @class Oldroyd_B_Fluid
 * @brief linear EOS with relaxation time and polymeric viscosity.
 */
class Oldroyd_B_Fluid : public WeaklyCompressibleFluid
{
  protected:
    Real lambda_; /**< relaxation time */
    Real mu_p_;   /**< polymeric viscosity */

  public:
    explicit Oldroyd_B_Fluid(Real rho0, Real c0, Real mu, Real lambda, Real mu_p)
        : WeaklyCompressibleFluid(rho0, c0, mu), lambda_(lambda), mu_p_(mu_p)
    {
        material_type_name_ = "Oldroyd_B_Fluid";
    };
    virtual ~Oldroyd_B_Fluid() {};

    Real getReferenceRelaxationTime() { return lambda_; };
    Real ReferencePolymericViscosity() { return mu_p_; };
    virtual Oldroyd_B_Fluid *ThisObjectPtr() override { return this; };
};

/**
 * @class GeneralizedNewtonianFluid
 * @brief Herschel Bulkley Model [for more information see: https://en.wikipedia.org/wiki/Generalized_Newtonian_fluid]
 */
class GeneralizedNewtonianFluid : public WeaklyCompressibleFluid
{
  protected:
    Real min_shear_rate_;
    Real max_shear_rate_;

  public:
    explicit GeneralizedNewtonianFluid(Real rho0, Real c0, Real min_shear_rate, Real max_shear_rate)
        : WeaklyCompressibleFluid(rho0, c0),
          min_shear_rate_(min_shear_rate), max_shear_rate_(max_shear_rate) {}

    virtual ~GeneralizedNewtonianFluid() {};

    virtual Real getViscosity(Real shear_rate) = 0;

    Real getMinShearRate() { return min_shear_rate_; };
    Real getMaxShearRate() { return max_shear_rate_; };
};

/**
 * @class HerschelBulkleyFluid
 * @brief https://en.wikipedia.org/wiki/Herschel%E2%80%93Bulkley_fluid
 */
class HerschelBulkleyFluid : public GeneralizedNewtonianFluid
{
  protected:
    Real consistency_index_;
    Real power_index_;
    Real yield_stress_;

  public:
    explicit HerschelBulkleyFluid(Real rho0, Real c0, Real min_shear_rate, Real max_shear_rate,
                                  Real consistency_index, Real power_index, Real yield_stress)
        : GeneralizedNewtonianFluid(rho0, c0, min_shear_rate, max_shear_rate),
          consistency_index_(consistency_index), power_index_(power_index), yield_stress_(yield_stress)
    {
        material_type_name_ = "HerschelBulkleyFluid";
    };
    virtual ~HerschelBulkleyFluid() {};

    Real getConsistencyIndex() { return consistency_index_; };
    Real getPowerIndex() { return power_index_; };
    Real getYieldStress() { return yield_stress_; };

    Real getViscosity(Real shear_rate) override;
    virtual HerschelBulkleyFluid *ThisObjectPtr() override { return this; };
};

/**
 * @class CarreauFluid
 * @brief https://en.wikipedia.org/wiki/Carreau_fluid
 */
class CarreauFluid : public GeneralizedNewtonianFluid
{
  protected:
    Real characteristic_time_;
    Real mu_infty_;
    Real mu0_;
    Real power_index_;

  public:
    explicit CarreauFluid(Real rho0, Real c0, Real min_shear_rate_, Real max_shear_rate_,
                          Real characteristic_time, Real mu_infty, Real mu0, Real power_index)
        : GeneralizedNewtonianFluid(rho0, c0, min_shear_rate_, max_shear_rate_),
          characteristic_time_(characteristic_time), mu_infty_(mu_infty),
          mu0_(mu0), power_index_(power_index)
    {
        material_type_name_ = "CarreauFluid";
    };
    virtual ~CarreauFluid() {};

    Real getCharacteristicTime() { return characteristic_time_; };
    Real getMuInfty() { return mu_infty_; };
    Real getMu0() { return mu0_; };
    Real getPowerIndex() { return power_index_; };

    Real getViscosity(Real shear_rate) override;
    virtual CarreauFluid *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_H