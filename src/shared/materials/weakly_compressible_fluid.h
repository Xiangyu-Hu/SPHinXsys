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
 * @file 	weakly_compressible_fluid.h
 * @brief 	Describe the weakly compressible fluid which is used
 * 			model incompressible fluids. Here, we have included several equation of states.
 * 			Furthermore, A typical non-newtonian fluid model is included.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef WEAKLY_COMPRESSIBLE_FLUID_H
#define WEAKLY_COMPRESSIBLE_FLUID_H

#include "base_material.h"
#include "common_functors.h"
#include "sphinxsys_constant.h"

namespace SPH
{
/**
 * @class WeaklyCompressibleFluid
 * @brief Linear equation of state (EOS).
 */
class WeaklyCompressibleFluid : public Fluid
{
  protected:
    Real rho0_; /**< reference density. */
    Real c0_;   /**< reference sound speed. */
    Real p0_;   /**< reference pressure */

  public:
    explicit WeaklyCompressibleFluid(Real rho0, Real c0);
    explicit WeaklyCompressibleFluid(ConstructArgs<Real, Real> args);
    virtual ~WeaklyCompressibleFluid(){};
    virtual Real ReferenceDensity() const override { return rho0_; };
    virtual Real ReferenceSoundSpeed() const override { return c0_; };
    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;

    class EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser)
            : rho0_(encloser.rho0_), p0_(encloser.p0_), c0_(encloser.c0_){};

        Real PressureFromDensity(UnsignedInt, Real rho)
        {
            return p0_ * (rho / rho0_ - 1.0);
        };

        Real DensityFromPressure(UnsignedInt, Real p)
        {
            return rho0_ * (p / p0_ + 1.0);
        };

        Real getSoundSpeed(UnsignedInt, Real, Real)
        {
            return c0_;
        };

        Real getReferenceDensity(UnsignedInt)
        {
            return rho0_;
        };

      protected:
        Real rho0_, p0_, c0_;
    };
};

class WeaklyCompressibleMixture : public Fluid
{
  protected:
    StdVec<std::string> species_name_list_; /**< species name list. */
    StdVec<Real> rho0_list_;                /**< reference density list. */
    ConstantArray<Real> *ca_inv_rho0_list_; /**< inverse reference density list. */
    DiscreteVariable<Real> *dv_Y_list_;     /**< species mass fraction list. */
    DiscreteVariable<Real> *dv_rho0_;       /**< local reference density. */
    Real c0_;                               /**< reference sound speed. */

  public:
    WeaklyCompressibleMixture(StdVec<std::string> species_name_list, StdVec<Real> rho0_list, Real c0);
    virtual ~WeaklyCompressibleMixture();
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Real ReferenceDensity() const override { return rho0_list_[0]; };
    virtual Real ReferenceSoundSpeed() const override { return c0_; };

    class EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser)
            : c0_(encloser.c0_), c0_sq_(c0_ * c0_),
              rho0_(encloser.dv_rho0_.DelegateDataView(ex_policy)){};

        Real PressureFromDensity(UnsignedInt index_i, Real rho)
        {
            return c0_sq_ * (rho - rho0_[index_i]);
        };

        Real DensityFromPressure(UnsignedInt index_i, Real p)
        {
            return rho0_[index_i] + p / c0_sq_;
        };

        Real getSoundSpeed(UnsignedInt, Real, Real)
        {
            return c0_;
        };

        Real getReferenceDensity(UnsignedInt index_i)
        {
            return rho0_[index_i];
        };

      protected:
        Real c0_, c0_sq_;
        DataView<Real> rho0_;
    };

    class MixtureKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        MixtureKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser)
            : inv_rho0_list_(encloser.ca_inv_rho0_list_.DelegatedData(ex_policy)),
              Y_list_(encloser.dv_Y_list_.DelegateMultiEntryView(ex_policy)),
              rho0_(encloser.dv_rho0_.DelegateDataView(ex_policy)){};

        void update(UnsignedInt index_i)
        {
            Real sum = 0.0;
            for (size_t k = 0; k != Y_list_.Width(); ++k)
            {
                sum += Y_list_[index_i][k] * inv_rho0_list_[k];
            }
            rho0_[index_i] = 1.0 / sum;
        };

      protected:
        Real *inv_rho0_list_;
        MultiEntryView<Real> Y_list_;
        DataView<Real> rho0_;
    };
};
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_H