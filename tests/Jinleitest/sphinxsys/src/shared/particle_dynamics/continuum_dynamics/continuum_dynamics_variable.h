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
 * @file continuum_dynamics_variable.h
 * @brief Here, we define the algorithm classes for computing derived solid dynamics variables.
 * @details These variable can be added into variable list for state output.
 * @author Shuaihao Zhang and Xiangyu Hu
 */

#ifndef CONTINUUM_DYNAMICS_VARIABLE_H
#define CONTINUUM_DYNAMICS_VARIABLE_H

#include "base_general_dynamics.h"
#include "general_continuum.h"
#include "general_continuum.hpp"

namespace SPH
{
namespace continuum_dynamics
{
/**
 * @class VonMisesStress
 * @brief computing von_Mises_stress
 */
class VonMisesStress : public BaseDerivedVariable<Real>
{
  public:
    explicit VonMisesStress(SPHBody &sph_body);
    virtual ~VonMisesStress(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *p_;
    Matd *shear_stress_;
};
/**
 * @class VonMisesStrain
 * @brief computing von_Mises_strain
 */
class VonMisesStrain : public BaseDerivedVariable<Real>
{
  public:
    explicit VonMisesStrain(SPHBody &sph_body);
    virtual ~VonMisesStrain(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *strain_tensor_;
};
/**
 * @class VerticalStress
 */
class VerticalStress : public BaseDerivedVariable<Real>
{
  public:
    explicit VerticalStress(SPHBody &sph_body);
    virtual ~VerticalStress(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Mat3d *stress_tensor_3D_;
};
/**
 * @class AccumulatedDeviatoricPlasticStrain
 */
class AccDeviatoricPlasticStrain : public BaseDerivedVariable<Real>
{
  public:
    explicit AccDeviatoricPlasticStrain(SPHBody &sph_body);
    virtual ~AccDeviatoricPlasticStrain(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    Mat3d *stress_tensor_3D_, *strain_tensor_3D_;
    Real E_, nu_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_DYNAMICS_VARIABLE_H