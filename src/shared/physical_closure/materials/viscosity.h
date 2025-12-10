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
 * @file 	viscosity.h
 * @brief 	tbd.
 * @author	Xiangyu Hu
 */

#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "base_data_type_package.h"
#include "particle_functors.h"

namespace SPH
{
class BaseParticles;

class Viscosity
{
  protected:
    Real mu_; /**< reference viscosity. */
  public:
    explicit Viscosity(Real mu) : mu_(mu) {};
    virtual ~Viscosity() {};
    Real ReferenceViscosity() { return mu_; };
    virtual void registerLocalParameters(BaseParticles *base_particles) {};
    virtual void registerLocalParametersFromReload(BaseParticles *base_particles) {};
    virtual void initializeLocalParameters(BaseParticles *base_particles) {};

    typedef ParameterFixed<Real> OneSideViscosity;

    template <typename ExecutionPolicy>
    ParameterFixed<Real> getOneSideViscosity(const ExecutionPolicy &ex_policy)
    {
        return ParameterFixed<Real>(mu_);
    };

    typedef PairGeomAverageFixed<Real> InterParticleViscosity;

    template <typename ExecutionPolicy>
    PairGeomAverageFixed<Real> getInterParticleViscosity(const ExecutionPolicy &ex_policy, Viscosity &other_viscosity)
    {
        return PairGeomAverageFixed<Real>(mu_, other_viscosity.ReferenceViscosity());
    };
};

class OldroydBViscosity : public Viscosity
{
  protected:
    Real lambda_; /**< relaxation time */
    Real mu_p_;   /**< polymeric viscosity */

  public:
    explicit OldroydBViscosity(Real mu, Real lambda, Real mu_p);
    explicit OldroydBViscosity(ConstructArgs<Real, Real, Real> args);
    virtual ~OldroydBViscosity() {};
    Real ReferenceRelaxationTime() { return lambda_; };
    Real ReferencePolymericViscosity() { return mu_p_; };
};

class GeneralizedNewtonianViscosity : public Viscosity
{
  protected:
    Real min_shear_rate_;
    Real max_shear_rate_;

  public:
    GeneralizedNewtonianViscosity(Real min_shear_rate, Real max_shear_rate);
    explicit GeneralizedNewtonianViscosity(ConstructArgs<Real, Real> args);
    virtual ~GeneralizedNewtonianViscosity() {};
    Real getMinShearRate() { return min_shear_rate_; };
    Real getMaxShearRate() { return max_shear_rate_; };
    virtual Real getViscosity(Real shear_rate) = 0;
};

/**
 * @class HerschelBulkleyViscosity
 * @brief https://en.wikipedia.org/wiki/Herschel%E2%80%93Bulkley_fluid
 */
class HerschelBulkleyViscosity : public GeneralizedNewtonianViscosity
{
  protected:
    Real consistency_index_;
    Real power_index_;
    Real yield_stress_;

  public:
    HerschelBulkleyViscosity(Real min_shear_rate, Real max_shear_rate,
                             Real consistency_index, Real power_index, Real yield_stress);
    explicit HerschelBulkleyViscosity(ConstructArgs<Real, Real, Real, Real, Real> args);
    virtual ~HerschelBulkleyViscosity() {};
    Real getConsistencyIndex() { return consistency_index_; };
    Real getPowerIndex() { return power_index_; };
    Real getYieldStress() { return yield_stress_; };
    Real getViscosity(Real shear_rate) override;
};

/**
 * @class CarreauViscosity
 * @brief https://en.wikipedia.org/wiki/Carreau_fluid
 */
class CarreauViscosity : public GeneralizedNewtonianViscosity
{
  protected:
    Real characteristic_time_;
    Real mu_infty_;
    Real mu0_;
    Real power_index_;

  public:
    CarreauViscosity(Real min_shear_rate_, Real max_shear_rate_,
                     Real characteristic_time, Real mu_infty, Real mu0, Real power_index);
    virtual ~CarreauViscosity() {};
    Real getCharacteristicTime() { return characteristic_time_; };
    Real getMuInfty() { return mu_infty_; };
    Real getMu0() { return mu0_; };
    Real getPowerIndex() { return power_index_; };
    Real getViscosity(Real shear_rate) override;
};
} // namespace SPH
#endif // VISCOSITY_H