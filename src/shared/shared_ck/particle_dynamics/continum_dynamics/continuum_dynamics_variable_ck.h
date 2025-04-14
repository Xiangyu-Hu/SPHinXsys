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
 * @file continuum_dynamics_variable_ck.h
 * @brief Here, we define the algorithm classes for computing derived solid dynamics variables.
 * @details These variable can be added into variable list for state output.
 * @author Shuang Li, Xiangyu Hu and Xiangyu Hu
 */

#ifndef CONTINUUM_DYNAMICS_VARIABLE_CK_H
#define CONTINUUM_DYNAMICS_VARIABLE_CK_H

#include "base_general_dynamics.h"
#include "general_continuum.h"
#include "general_continuum.hpp"

namespace SPH
{
namespace continuum_dynamics
{
/**
 * @class VerticalStress
 */
class VerticalStressCK : public BaseDerivedVariable<Real>
{
  public:
    explicit VerticalStressCK(SPHBody &sph_body);
    virtual ~VerticalStressCK(){};
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Mat3d *stress_tensor_3D_;
        Real *derived_variable_;
    };
  protected:
    DiscreteVariable<Mat3d> *dv_stress_tensor_3D_;
    DiscreteVariable<Real> *dv_derived_variable_;
};

/**
 * @class AccumulatedDeviatoricPlasticStrain
 */
class AccDeviatoricPlasticStrainCK : public BaseDerivedVariable<Real>
{
  using PlasticKernel = typename PlasticContinuum::PlasticKernel;
  public:
    explicit AccDeviatoricPlasticStrainCK(SPHBody &sph_body);
    virtual ~AccDeviatoricPlasticStrainCK(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        PlasticKernel plastic_kernel_;
        Mat3d *stress_tensor_3D_, *strain_tensor_3D_;
        Real *derived_variable_;
        Real E_, nu_;
    };

  protected:
    PlasticContinuum &plastic_continuum_;
    DiscreteVariable<Mat3d> *dv_stress_tensor_3D_, *dv_strain_tensor_3D_;
    DiscreteVariable<Real> *dv_derived_variable_;
    Real E_, nu_;
};

} // namespace continuum_dynamics
} // namespace SPH

#endif // CONTINUUM_DYNAMICS_VARIABLE_CK_H