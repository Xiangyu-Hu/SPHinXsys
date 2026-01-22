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
 * @file 	shear_integration.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details CK and SYCL version.
 * @author	Shuang Li,Xiangyu Hu and Shuaihao Zhang
 */
#ifndef SHEAR_INTEGRATION_H
#define SHEAR_INTEGRATION_H

#include "acoustic_step_1st_half.h"
#include "general_continuum.h"
#include "general_continuum.hpp"
#include "general_gradient.h"
namespace SPH
{
namespace continuum_dynamics
{

template <typename...>
class ShearForce;

template <class MaterialType, typename... Parameters>
class ShearForce<Inner<WithInitialization, MaterialType, Parameters...>>
    : public Interaction<Inner<Parameters...>>
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;
    using MaterialKernel = typename MaterialType::MaterialKernel;

  public:
    explicit ShearForce(Inner<Parameters...> &inner_relation);
    virtual ~ShearForce() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        MaterialKernel material_;
        Matd *velocity_gradient_, *shear_stress_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Matd *shear_stress_;
        Real *scale_penalty_force_, xi_;
    };

  protected:
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // SHEAR_INTEGRATION_H