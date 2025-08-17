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
 * @file neighbor_method.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef NEIGHBOR_METHOD_H
#define NEIGHBOR_METHOD_H

#include "kernel_tabulated_ck.h"
#include "sphinxsys_containers.h"

namespace SPH
{
class BodyPartitionSpatial;

template <typename...>
class NeighborMethod;

template <>
class NeighborMethod<Base>
{
  public:
    NeighborMethod(Kernel &base_kernel) : base_kernel_(base_kernel) {};

  protected:
    Kernel &base_kernel_;

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
            .template registerStateVariable<Real>(
                "SmoothingLength",
                [&](size_t i)
                { return 0.99 * h_spacing_ratio * identifier.getBaseParticles().ParticleSpacing(i); });
    };
};

template <>
class NeighborMethod<SingleValued> : public NeighborMethod<Base>
{
  public:
    template <class SourceIdentifier, class TargetIdentifier>
    NeighborMethod(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier)
        : NeighborMethod<Base>(*source_identifier.getSPHAdaptation().getKernel())
    {
        Real source_h = getSmoothingLength(SingleValued{}, source_identifier);
        Real target_h = getSmoothingLength(SingleValued{}, contact_identifier);
        inv_h_ = 1.0 / SMAX(source_h, target_h);
        search_depth_ = static_cast<int>(std::ceil((source_h - Eps) / target_h));
    }

    class SmoothingKernel : public KernelTabulatedCK
    {
        Real inv_h_, inv_h_squared_, inv_h_cubed_, inv_h_fourth_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : KernelTabulatedCK(encloser.base_kernel_),
              inv_h_(encloser.inv_h_), inv_h_squared_(inv_h_ * inv_h_),
              inv_h_cubed_(inv_h_squared_ * inv_h_), inv_h_fourth_(inv_h_cubed_ * inv_h_){};

        inline Real W(const Vecd &displacement) const
        {
            return Factor(displacement) * normalized_W((displacement * inv_h_).norm());
        };

        inline Real dW(const Vecd &displacement) const
        {
            return GradientFactor(displacement) * normalized_dW((displacement * inv_h_).norm());
        };

      protected:
        inline Real Factor(const Vec2d &) const { return inv_h_squared_ * dimension_factor_2D_; };
        inline Real Factor(const Vec3d &) const { return inv_h_cubed_ * dimension_factor_3D_; };
        inline Real GradientFactor(const Vec2d &) const { return inv_h_cubed_ * dimension_factor_2D_; };
        inline Real GradientFactor(const Vec3d &) const { return inv_h_fourth_ * dimension_factor_3D_; };
    };

    class CriterionKernel
    {
        Vecd *source_pos_;
        Vecd *target_pos_;
        Real kernel_size_squared_, inv_h_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        CriterionKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                        DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos)
            : source_pos_(dv_source_pos->DelegatedData(ex_policy)),
              target_pos_(dv_target_pos->DelegatedData(ex_policy)),
              kernel_size_squared_(math::pow(encloser.base_kernel_.KernelSize(), 2)),
              inv_h_(encloser.inv_h_) {}

        inline bool operator()(UnsignedInt i, UnsignedInt j) const
        {
            return (inv_h_ * (source_pos_[i] - target_pos_[j])).squaredNorm() < kernel_size_squared_;
        };
    };

    class SmoothingRatioKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingRatioKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser){};

        Real operator()(UnsignedInt i) const { return 1.0; };
    };

    class SearchDepthKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        SearchDepthKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : search_depth_(encloser.search_depth_){};

        int operator()(UnsignedInt i) const { return search_depth_; };

      private:
        int search_depth_;
    };

  protected:
    Real inv_h_;
    int search_depth_; /**< Search depth for neighbor search. */
};
} // namespace SPH
#endif // NEIGHBOR_METHOD_H
