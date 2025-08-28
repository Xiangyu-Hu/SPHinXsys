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
 * 			Furthermore, A typical non-newtonian fluid model is included.
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
    explicit WeaklyCompressibleFluid(Real rho0, Real c0);
    explicit WeaklyCompressibleFluid(ConstructArgs<Real, Real> args);
    virtual ~WeaklyCompressibleFluid() {};

    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;

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
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_H