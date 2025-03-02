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
 * @file 	initilization_dynamics_ck.h
 * @brief 	Here, we define the ck_version for stress diffusion. 
 * @details Refer to Feng et al(2021).
 * @author	Shuang Li, Xiangyu Hu and Shuaihao Zhang
 */
#ifndef INITILIZATION_DYNAMICS_CK_H
#define INITILIZATION_DYNAMICS_CK_H

#include "base_continuum_dynamics.h"
#include "constraint_dynamics.h"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "general_continuum.hpp"

namespace SPH
{
namespace continuum_dynamics
{
class ContinuumInitialConditionCK : public LocalDynamics
{
  public:
    explicit ContinuumInitialConditionCK(SPHBody &sph_body);
    virtual ~ContinuumInitialConditionCK(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0){};
      protected:
        Vecd *pos_, *vel_;
        Mat3d *stress_tensor_3D_;
    };

  protected:

    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
    DiscreteVariable<Mat3d> *dv_stress_tensor_3D_;
};

} // namespace continuum_dynamics
} // namespace SPH
#endif // INITILIZATION_DYNAMICS_CK_H