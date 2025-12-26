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
 * @file relaxation_stepping_ck.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef RELAXATION_STEPPING_CK_H
#define RELAXATION_STEPPING_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
class RelaxationScalingCK : public LocalDynamicsReduce<ReduceMax>
{
  public:
    RelaxationScalingCK(SPHBody &sph_body);
    virtual ~RelaxationScalingCK() {};

    class FinishDynamics
    {
        Real h_ref_;

      public:
        using OutputType = Real;
        FinishDynamics(RelaxationScalingCK &encloser);
        Real Result(Real reduced_value);
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy>
        ReduceKernel(const ExecutionPolicy &ex_policy, RelaxationScalingCK &encloser)
            : residual_(encloser.dv_residual_->DelegatedData(ex_policy)),
              h_ref_(encloser.h_ref_){};

        Real reduce(size_t index_i, Real dt)
        {
            return residual_[index_i].norm();
        };

      protected:
        Vecd *residual_;
        Real h_ref_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_residual_;
    Real h_ref_;
};

class PositionRelaxationCK : public LocalDynamics
{
  public:
    explicit PositionRelaxationCK(SPHBody &sph_body);
    virtual ~PositionRelaxationCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy, PositionRelaxationCK &encloser)
            : pos_(encloser.pos_->DelegatedData(ex_policy)),
              residual_(encloser.residual_->DelegatedData(ex_policy)){};

        void update(size_t index_i, Real dt_square)
        {
            pos_[index_i] += residual_[index_i] * dt_square * 0.5;
        };

      protected:
        Vecd *pos_, *residual_;
    };

  protected:
    DiscreteVariable<Vecd> *pos_, *residual_;
};

template <class DynamicIdentifier>
class UpdateSmoothingLengthRatio : public BaseLocalDynamics<DynamicIdentifier>
{
    using SpacingAdaptation = typename DynamicIdentifier::SpacingAdaptation;
    using AdaptationKerenl = typename SpacingAdaptation::ComptutingKernel;

  public:
    template <typename... Args>
    UpdateSmoothingLengthRatio(DynamicIdentifier &identfier, Args &&...args)
        : BaseLocalDynamics<DynamicIdentifier>(identfier),
          dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
          dv_h_ratio_(this->particles_->template getVariableByName<Real>("SmoothingLengthRatio")),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          adaptation_method_(identfier.getAdaptation(), std::forward<Args>(args)...),
          reference_spacing_(identfier.getAdaptation().ReferenceSpacing()){};
    virtual ~UpdateSmoothingLengthRatio() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              h_ratio_(encloser.dv_h_ratio_->DelegatedData(ex_policy)),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              spacing_adaptation_(encloser.adaptation_method_),
              reference_spacing_(encloser.reference_spacing_){};

        void update(size_t index_i, Real dt = 0.0)
        {
            Real local_spacing = spacing_adaptation_(pos_[index_i]);
            h_ratio_[index_i] = reference_spacing_ / local_spacing;
            Vol_[index_i] = math::pow(local_spacing, Dimensions);
        };

      protected:
        Vecd *pos_;
        Real *h_ratio_, *Vol_;
        AdaptationKerenl spacing_adaptation_;
        Real reference_spacing_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_h_ratio_, *dv_Vol_;
    SpacingAdaptation adaptation_method_;
    Real reference_spacing_;
};
} // namespace SPH
#endif // RELAXATION_STEPPING_CK_H
