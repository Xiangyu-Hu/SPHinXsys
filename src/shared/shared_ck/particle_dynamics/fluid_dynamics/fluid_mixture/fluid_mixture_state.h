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
 * @file 	fluid_mixture_state.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_MIXTURE_STATE_H
#define FLUID_MIXTURE_STATE_H

#include "base_fluid_dynamics.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename FluidMixtureType>
class ConstantMixtureFraction : public ReturnFunction<Real>
{
    using EosKernel = typename FluidMixtureType::EosKernel;

  public:
    ConstantMixtureFraction(
        BaseParticles *base_particles,
        FluidMixtureType &fluid_mixture, StdVec<Real> mixture_fractions)
        : fluid_mixture_(fluid_mixture),
          ca_mixture_fractions_("ConstantMixtureFractions", mixture_fractions)
    {
        if (fluid_mixture.NumberOfMixtures() != mixture_fractions.size())
        {
            std::cout << "\n Error: the mixture fraction list and the number of "
                      << "mixtures should be the same !" << std::endl;
            exit(1);
        }
    };
    ~ConstantMixtureFraction() {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : eos_(ex_policy, encloser.fluid_mixture_),
              mixture_fractions_(
                  encloser.ca_mixture_fractions_.DelegatedArrayView(ex_policy)){};

        void operator()(size_t index_i)
        {
            eos_.setMixtureFractions(index_i, mixture_fractions_);
        };

      protected:
        EosKernel eos_;
        ArrayView<Real> mixture_fractions_;
    };

  protected:
    FluidMixtureType &fluid_mixture_;
    ConstantArray<Real> ca_mixture_fractions_;
};

template <typename FluidMixtureType>
class UpdateReferenceDensity : public ReturnFunction<Real>
{
    using EosKernel = typename FluidMixtureType::EosKernel;

  public:
    UpdateReferenceDensity(
        BaseParticles *base_particles, FluidMixtureType &fluid_mixture)
        : fluid_mixture_(fluid_mixture),
          dv_rho0_(fluid_mixture.dvReferenceDensity()),
          dv_Vol_ref_(base_particles->getVariableByName<Real>(
              "VolumetricMeasureRef")),
          dv_mass_(base_particles->getVariableByName<Real>("Mass")) {};
    ~UpdateReferenceDensity() {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : eos_(ex_policy, encloser.fluid_mixture_),
              rho0_(encloser.dv_rho0_->DelegatedDataView(ex_policy)),
              Vol_ref_(encloser.dv_Vol_ref_->DelegatedDataView(ex_policy)),
              mass_(encloser.dv_mass_->DelegatedDataView(ex_policy)){};

        void operator()(size_t index_i)
        {
            rho0_[index_i] = eos_.computeReferenceDensity(index_i);
            mass_[index_i] = rho0_[index_i] * Vol_ref_[index_i];
        };

      protected:
        EosKernel eos_;
        DataView<Real> rho0_, Vol_ref_, mass_;
    };

  protected:
    FluidMixtureType &fluid_mixture_;
    DiscreteVariable<Real> *dv_rho0_, *dv_Vol_ref_, *dv_mass_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_MIXTURE_STATE_H
