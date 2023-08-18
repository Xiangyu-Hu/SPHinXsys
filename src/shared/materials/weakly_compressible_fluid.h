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
    virtual ~WeaklyCompressibleFluid(){};

    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;
    virtual WeaklyCompressibleFluid *ThisObjectPtr() override { return this; };
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
    template <typename... ConstructorArgs>
    explicit WeaklyCompressibleFluidFreeSurface(Real cutoff_pressure, ConstructorArgs &&...args)
        : WeaklyCompressibleFluidType(std::forward<ConstructorArgs>(args)...),
          cutoff_pressure_(cutoff_pressure),
          cutoff_density_(WeaklyCompressibleFluidType::DensityFromPressure(cutoff_pressure))
    {
        WeaklyCompressibleFluidType::material_type_ += "FreeSurface";
    };
    virtual ~WeaklyCompressibleFluidFreeSurface(){};

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
    virtual ~SymmetricTaitFluid(){};

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
    virtual ~Oldroyd_B_Fluid(){};

    Real getReferenceRelaxationTime() { return lambda_; };
    Real ReferencePolymericViscosity() { return mu_p_; };
    virtual Oldroyd_B_Fluid *ThisObjectPtr() override { return this; };
};
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_H