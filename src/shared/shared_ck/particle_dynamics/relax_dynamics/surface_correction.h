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
 * @file surface_correction.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef SURFACE_CORRECTION_H
#define SURFACE_CORRECTION_H

#include "base_body_part.h"
#include "base_general_dynamics.h"

namespace SPH
{
class LevelsetBounding : public BaseLocalDynamics<BodyPartByCell>
{
    using ProbeSignedDistance = LevelSet::ProbeLevelSet<Real>;
    using ProbeLevelsetGradient = LevelSet::ProbeLevelSet<Vecd>;

  public:
    LevelsetBounding(NearShapeSurface &body_part);
    virtual ~LevelsetBounding() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0)
        {
            Real phi = signed_distance_(pos_[index_i]);

            if (phi > -constrained_distance_)
            {
                pos_[index_i] -= (phi + constrained_distance_) *
                                 level_set_gradient_(pos_[index_i]).normalized();
            }
        };

      protected:
        Vecd *pos_;
        ProbeSignedDistance signed_distance_;
        ProbeLevelsetGradient level_set_gradient_;
        Real constrained_distance_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    LevelSet &level_set_;
    Real constrained_distance_;
};

class LevelsetKernelGradientIntegral : public LocalDynamics
{
    using ProbeKernelGradientIntegral = LevelSet::ProbeLevelSet<Vecd>;

  public:
    LevelsetKernelGradientIntegral(SPHBody &sph_body, LevelSetShape &level_set_shape);
    virtual ~LevelsetKernelGradientIntegral() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0)
        {
            residual_[index_i] -= 2.0 * kernel_gradient_integral_(pos_[index_i], 1.0);
        };

      protected:
        Vecd *pos_, *residual_;
        ProbeKernelGradientIntegral kernel_gradient_integral_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Vecd> *dv_residual_;
    LevelSet &level_set_;
};

} // namespace SPH
#endif // LEVEL_SET_CORRECTION_H
