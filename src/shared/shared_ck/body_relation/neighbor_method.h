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

#include "adaptation.h"
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
    NeighborMethod(SharedPtr<Kernel> base_kernel) : base_kernel_(base_kernel) {};

    class SmoothingKernel : public KernelTabulatedCK
    {
      public:
        SmoothingKernel(NeighborMethod<Base> &encloser) : KernelTabulatedCK(*encloser.base_kernel_) {};

      protected:
        inline Real W(const Real &inv_h_squared, const Vec2d &displacement, const Real &inv_h) const;
        inline Real W(const Real &inv_h_cubed, const Vec3d &displacement, const Real &inv_h) const;
        inline Real dW(const Real &inv_h_cubed, const Vec2d &displacement, const Real &inv_h) const;
        inline Real dW(const Real &inv_h_fourth, const Vec3d &displacement, const Real &inv_h) const;
        inline Real d2W(const Real &inv_h_fourth, const Vec2d &displacement, const Real &inv_h) const;
        inline Real d2W(const Real &inv_h_fifth, const Vec3d &displacement, const Real &inv_h) const;
    };

  protected:
    SharedPtr<Kernel> base_kernel_;
};

template <>
class NeighborMethod<SPHAdaptation, SPHAdaptation> : public NeighborMethod<Base>
{
    using BaseKernel = NeighborMethod<Base>::SmoothingKernel;

  public:
    template <typename SourceIdentifier, typename TargetIdentifier>
    NeighborMethod(SourceIdentifier &source_identifier, TargetIdentifier &target_identifier);
    NeighborMethod(SharedPtr<Kernel> base_kernel, Real h, Real search_increment);

    class SmoothingKernel : public BaseKernel
    {
        Real inv_h_, inv_h_squared_, inv_h_cubed_, inv_h_fourth_, inv_h_fifth_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        inline Real W(const Vec2d &displacement, UnsignedInt, UnsignedInt) const;
        inline Real W(const Vec3d &displacement, UnsignedInt, UnsignedInt) const;
        inline Real dW(const Vec2d &displacement, UnsignedInt, UnsignedInt) const;
        inline Real dW(const Vec3d &displacement, UnsignedInt, UnsignedInt) const;
        inline Real d2W(const Vec2d &displacement, UnsignedInt, UnsignedInt) const;
        inline Real d2W(const Vec3d &displacement, UnsignedInt, UnsignedInt) const;
       
        Real CutOffRadius() const { return kernel_size_ / inv_h_; }
    };

    class NeighborCriterion
    {
        Vecd *source_pos_;
        Vecd *target_pos_;
        Real kernel_size_squared_, inv_h_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                          DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos);
        inline bool operator()(UnsignedInt i, UnsignedInt j) const;
    };

    class ReverseNeighborCriterion
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReverseNeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser){};

        inline bool operator()(UnsignedInt i, UnsignedInt j) const { return true; };
    };

    class SmoothingRatio
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingRatio(const ExecutionPolicy &ex_policy, EncloserType &encloser){};

        inline Real operator()(UnsignedInt i) const { return 1.0; };
    };

    class SearchBox
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        SearchBox(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        inline BoundingBoxi operator()(UnsignedInt i) const { return search_box_; };

      private:
        BoundingBoxi search_box_;
    };

  protected:
    Real inv_h_;
    int search_depth_;
    BoundingBoxi search_box_; /**< Search depth for neighbor search. */
};

template <>
class NeighborMethod<AdaptiveSmoothingLength, AdaptiveSmoothingLength> : public NeighborMethod<Base>
{
  public:
    template <typename SourceIdentifier, typename TargetIdentifier>
    NeighborMethod(SourceIdentifier &source_identifier, TargetIdentifier &target_identifier);

    class SmoothingKernel : public NeighborMethod<Base>::SmoothingKernel
    {
        using BaseKernel = NeighborMethod<Base>::SmoothingKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        inline Real W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const;
        inline Real W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const;
        inline Real dW(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const;
        inline Real dW(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const;
        inline Real d2W(const Vec2d &displacement, UnsignedInt i, UnsignedInt j) const;
        inline Real d2W(const Vec3d &displacement, UnsignedInt i, UnsignedInt j) const;

      protected:
        Real src_inv_h_ref_, tar_inv_h_ref_;
        Real *src_h_ratio_, *tar_h_ratio_;

        inline Real invH(UnsignedInt i, UnsignedInt j) const;
    };

  protected:
    Real h_ref_;
    DiscreteVariable<Vecd> *dv_src_h_ratio_, *dv_tar_h_ratio_;
};
} // namespace SPH
#endif // NEIGHBOR_METHOD_H
