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
 * @file 	porous_media_solid.h
 * @brief 	These are classes for define properties of elastic solid materials.
 *			These classes are based on isotropic linear elastic solid.
 * 			Several more complex materials, including neo-Hookean, FENE neo-Hookean
 *			and anisotropic muscle, are derived from the basic elastic solid class.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef POROUS_SOLID_H
#define POROUS_SOLID_H

#include "elastic_solid.h"
#include <fstream>

namespace SPH
{
namespace multi_species_continuum
{
/**
 * @class PorousMediaSolid
 * @brief Abstract class for a porous media solid
 * Note that porous media can be plastic solid if necessary
 */
class PorousMediaSolid : public LinearElasticSolid
{
  public:
    /*< material parameter */
    Real fluid_initial_density_;
    Real diffusivity_constant_;
    Real water_pressure_constant_;

    Real getFluidInitialDensity() { return fluid_initial_density_; };
    Real getDiffusivityConstant() { return diffusivity_constant_; };
    Real getWaterPressureConstant() { return water_pressure_constant_; };

    /** Constructor */
    explicit PorousMediaSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real diffusivity_constant,
                              Real fluid_initial_density, Real water_pressure_constant)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio),
          fluid_initial_density_(fluid_initial_density),
          diffusivity_constant_(diffusivity_constant),
          water_pressure_constant_(water_pressure_constant)
    {
        material_type_name_ = "PorousMediaSolid";
    };

    virtual ~PorousMediaSolid() {};
};

} // namespace multi_species_continuum
} // namespace SPH
#endif // POROUS_SOLID_H
