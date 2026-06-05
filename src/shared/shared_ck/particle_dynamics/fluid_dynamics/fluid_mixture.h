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
 * @file 	fluid_mixture.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_MIXTURE_H
#define FLUID_MIXTURE_H

#include "base_fluid_dynamics.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
class ReferenceDesnitySetup : public LocalDynamics
{
  public:
    explicit ReferenceDesnitySetup(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          fluid_mixture_(DynamicCast<WeaklyCompressibleMixture>(
              this, this->sph_body_->getMatterMaterial())),
          ca_inv_rho0_list_(fluid_mixture_.caInvReferenceDensity()),
          dv_Y_list_(fluid_mixture_.dvMassFraction()),
          dv_rho0_(fluid_mixture_.dvReferenceDensity()),
          dv_mass_(fluid_mixture_.dvMass()){};
    virtual ~ReferenceDesnitySetup(){};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, ReferenceDesnitySetup &encloser)
            : inv_rho0_list_(encloser.ca_inv_rho0_list_->DelegatedDataView(ex_policy)),
              Y_list_(encloser.dv_Y_list_->DelegatedMultiEntryView(ex_policy)),
              rho0_(encloser.dv_rho0_->DelegatedDataView(ex_policy)),
              mass_(encloser.dv_mass_->DelegatedDataView(ex_policy)){};

        void update(size_t index_i, Real dt = 0.0)
        {
            Real sum = 0.0;
            Real old_rho0_inv = 1.0 / rho0_[index_i];
            Real old_mass = mass_[index_i];
            for (size_t k = 0; k != Y_list_.Width(); ++k)
            {
                sum += Y_list_[index_i][k] * inv_rho0_list_[k];
            }
            rho0_[index_i] = 1.0 / sum;
            mass_[index_i] = rho0_[index_i] * old_mass * old_rho0_inv;
        };

      protected:
        DataView<Real> inv_rho0_list_;
        MultiEntryView<Real> Y_list_;
        DataView<Real> rho0_, mass_;
    };

  protected:
    WeaklyCompressibleMixture &fluid_mixture_;
    ConstantArray<Real> *ca_inv_rho0_list_;
    DiscreteVariable<Real> *dv_Y_list_;
    DiscreteVariable<Real> *dv_rho0_, *dv_mass_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_MIXTURE_H
