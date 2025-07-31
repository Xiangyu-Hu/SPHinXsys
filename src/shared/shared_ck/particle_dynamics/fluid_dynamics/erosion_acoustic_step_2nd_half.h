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
 * @file 	acoustic_step_2nd_half.h
 * @brief 	Here, we define the algorithm classes for erosion fluid dynamics.
 * @details Erosion dynamics.
 * @author	Shuang Li & Xiangyu Hu
 */

#ifndef EROSION_ACOUSTIC_STEP_2ND_HALF_H
#define EROSION_ACOUSTIC_STEP_2ND_HALF_H

#include "erosion_acoustic_step_1st_half.h"

namespace SPH
{
namespace fluid_dynamics
{
    
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep2ndHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>, public Interaction<Soil>
{
    using FluidType = typename RiemannSolverType::SourceFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep2ndHalf(Contact<Parameters...> &soil_contact_relation);
    virtual ~AcousticStep2ndHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
        Real *soil_Vol_, *soil_p_;
        Vecd *soil_n_,*soil_vel_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // EROSION_ACOUSTIC_STEP_2ND_HALF_H
