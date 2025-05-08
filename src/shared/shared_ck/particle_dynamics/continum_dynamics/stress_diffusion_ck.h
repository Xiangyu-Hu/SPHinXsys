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
 * @file 	stress_diffusion.h
 * @brief 	Here, we define the ck_version for stress diffusion.
 * @details Refer to Feng et al(2021).
 * @author	Shuang Li, Xiangyu Hu and Shuaihao Zhang
 */

#ifndef STRESS_DIFFUSION_CK_H
#define STRESS_DIFFUSION_CK_H

#include "base_continuum_dynamics.h"
#include "constraint_dynamics.h"
#include "continuum_integration_1st_ck.h"
#include "continuum_integration_1st_ck.hpp"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "general_continuum.hpp"
namespace SPH
{
namespace continuum_dynamics
{
template <typename...>
class StressDiffusionCK;

template <typename... Parameters>
class StressDiffusionCK<Inner<Parameters...>> : public PlasticAcousticStep<Interaction<Inner<Parameters...>>>
{
    using PlasticKernel = typename PlasticContinuum::PlasticKernel;
    using BaseInteraction = PlasticAcousticStep<Interaction<Inner<Parameters...>>>;

  public:
    explicit StressDiffusionCK(Inner<Parameters...> &inner_relation);
    virtual ~StressDiffusionCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        PlasticKernel plastic_kernel_;
        Real zeta_, phi_;
        Real smoothing_length_, sound_speed_;
        Real *mass_, *Vol_;
        Vecd *pos_, *force_prior_;
        Mat3d *stress_tensor_3D_, *stress_rate_3D_;
    };

  protected:
    Real dv_zeta_ = 0.1, dv_phi_; /*diffusion coefficient*/
    Real dv_smoothing_length_, dv_sound_speed_;
    DiscreteVariable<Vecd> *dv_pos_;
};

using StressDiffusionInnerCK = StressDiffusionCK<Inner<>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // STRESS_DIFFUSION_CK_H