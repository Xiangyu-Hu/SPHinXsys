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
 * @file neighbor_method.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef NEIGHBOR_METHOD_H
#define NEIGHBOR_METHOD_H

#include "base_body.h"
#include "base_particles.hpp"

namespace SPH
{
class BodyPartitionSpatial;

template <typename...>
class SmoothingLength;

template <>
class SmoothingLength<Base>
{
  public:
    template <typename... Args>
    SmoothingLength(Args &&...args){};

  protected:
    template <class DynamicsIdentifier>
    Real getSmoothingLength(const SingleValued &single_valued, DynamicsIdentifier &identifier)
    {
        return identifier.getSPHAdaptation().ReferenceSmoothingLength();
    };
    Real getSmoothingLength(const SingleValued &single_valued, BodyPartitionSpatial &body_partition_spatial);

    template <class DynamicsIdentifier>
    DiscreteVariable<Real> *getSmoothingLength(const Continuous &continuous, DynamicsIdentifier &identifier)
    {
        Real h_spacing_ratio = identifier.getSPHAdaptation().SmoothingLengthSpacingRatio();
        return identifier.getBaseParticles()
            .template registerStateVariableOnly<Real>(
                "SmoothingLength",
                [&](size_t i)
                { return 0.99 * h_spacing_ratio * identifier.getBaseParticles().ParticleSpacing(i); });
    };
};

template <>
class SmoothingLength<SingleValued> : public SmoothingLength<Base>
{
  public:
    template <class DynamicsIdentifier>
    SmoothingLength(DynamicsIdentifier &identifier)
        : SmoothingLength<Base>(identifier),
          inv_h_(1.0 / getSmoothingLength(SingleValued{}, identifier)){};

    template <class SourceIdentifier, class TargetIdentifier>
    SmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : SmoothingLength<Base>(source_identifier, contact_identifier),
          inv_h_(1.0 / SMAX(getSmoothingLength(SingleValued{}, source_identifier),
                            getSmoothingLength(SingleValued{}, contact_identifier))){};

    class ComputingKernel
    {
        Real inv_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, SmoothingLength<SingleValued> &smoothing_length)
            : inv_h_(smoothing_length.inv_h_){};
        Real operator()(UnsignedInt i, UnsignedInt j) const { return inv_h_; };
    };

  protected:
    Real inv_h_;
};

template <>
class SmoothingLength<Continuous> : public SmoothingLength<Base>
{
  public:
    template <class DynamicsIdentifier>
    SmoothingLength(DynamicsIdentifier &identifier)
        : SmoothingLength<Base>(identifier),
          dv_source_h_(getSmoothingLength(Continuous{}, identifier)),
          dv_target_h_(dv_source_h_){};

    template <class SourceIdentifier, class TargetIdentifier>
    SmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : SmoothingLength<Base>(source_identifier, contact_identifier),
          dv_source_h_(getSmoothingLength(Continuous{}, source_identifier)),
          dv_target_h_(getSmoothingLength(Continuous{}, contact_identifier)){};

    class ComputingKernel
    {
        Real *source_h_;
        Real *target_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, SmoothingLength<Continuous> &smoothing_length)
            : source_h_(smoothing_length.dv_source_h_->DelegatedData(ex_policy)),
              target_h_(smoothing_length.dv_target_h_->DelegatedData(ex_policy)){};
        Real operator()(UnsignedInt i, UnsignedInt j) const { return 1.0 / SMAX(source_h_[i], target_h_[j]); };
    };

  protected:
    DiscreteVariable<Real> *dv_source_h_;
    DiscreteVariable<Real> *dv_target_h_;
};

template <>
class SmoothingLength<SingleValued, Continuous> : public SmoothingLength<Base>
{
  public:
    template <class SourceIdentifier, class TargetIdentifier>
    SmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : SmoothingLength<Base>(source_identifier, contact_identifier),
          source_h_(getSmoothingLength(SingleValued{}, source_identifier)),
          dv_target_h_(getSmoothingLength(Continuous{}, contact_identifier)){};

    class ComputingKernel
    {
        Real source_h_;
        Real *target_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, SmoothingLength<SingleValued, Continuous> &smoothing_length)
            : source_h_(smoothing_length.source_h_),
              target_h_(smoothing_length.dv_target_h_->DelegatedData(ex_policy)){};
        Real operator()(UnsignedInt i, UnsignedInt j) const { return 1.0 / SMAX(source_h_, target_h_[j]); };
    };

  protected:
    Real source_h_;
    DiscreteVariable<Real> *dv_target_h_;
};

template <>
class SmoothingLength<Continuous, SingleValued> : public SmoothingLength<Base>
{
  public:
    template <class SourceIdentifier, class TargetIdentifier>
    SmoothingLength(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : SmoothingLength<Base>(source_identifier, contact_identifier),
          dv_source_h_(getSmoothingLength(Continuous{}, source_identifier)),
          target_h_(getSmoothingLength(SingleValued{}, contact_identifier)){};

    class ComputingKernel
    {
        Real *source_h_;
        Real target_h_;

      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy, SmoothingLength<Continuous, SingleValued> &smoothing_length)
            : source_h_(smoothing_length.dv_source_h_->DelegatedData(ex_policy)),
              target_h_(smoothing_length.target_h_){};
        Real operator()(UnsignedInt i, UnsignedInt j) const { return 1.0 / SMAX(source_h_[i], target_h_); };
    };

  protected:
    DiscreteVariable<Real> *dv_source_h_;
    Real target_h_;
};
} // namespace SPH
#endif // NEIGHBOR_METHOD_H