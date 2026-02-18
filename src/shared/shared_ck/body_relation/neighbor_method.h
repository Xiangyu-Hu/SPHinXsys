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
class SPHAdaptation;

template <typename...>
class Neighbor;

template <>
class Neighbor<Base>
{
  public:
    Neighbor(SharedPtr<Kernel> base_kernel,
             DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos);
    Neighbor(SharedPtr<Kernel> base_kernel);
    ~Neighbor() {};

    class SmoothingKernel : public KernelTabulatedCK
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        SmoothingKernel(Neighbor<Base> &encloser);

        inline Vecd vec_r_ij(UnsignedInt i, UnsignedInt j) const { return src_pos_[i] - tar_pos_[j]; };
        inline Vecd e_ij(UnsignedInt i, UnsignedInt j) const { return vec_r_ij(i, j).normalized(); };

      protected:
        Vecd *src_pos_, *tar_pos_;

        Real W2D(const Real &inv_h_squared, const Real &scaled_r) const;
        Real W3D(const Real &inv_h_cubed, const Real &scaled_r) const;
        Real W02D(const Real &inv_h_squared) const;
        Real W03D(const Real &inv_h_cubed) const;
        Real dW2D(const Real &inv_h_cubed, const Real &scaled_r) const;
        Real dW3D(const Real &inv_h_fourth, const Real &scaled_r) const;
        Real d2W2D(const Real &inv_h_fourth, const Real &scaled_r) const;
        Real d2W3D(const Real &inv_h_fifth, const Real &scaled_r) const;
    };

  protected:
    SharedPtr<Kernel> base_kernel_;
    DiscreteVariable<Vecd> *dv_src_pos_;
    DiscreteVariable<Vecd> *dv_tar_pos_;
};

template <>
class Neighbor<SPHAdaptation, SPHAdaptation> : public Neighbor<Base>
{
    using BaseKernel = Neighbor<Base>::SmoothingKernel;

  public:
    template <typename SourceIdentifier, typename TargetIdentifier>
    Neighbor(SourceIdentifier &source_identifier, TargetIdentifier &target_identifier,
             DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos);
    Neighbor(SharedPtr<Kernel> base_kernel, Real h, Real search_increment);

    class SmoothingKernel : public BaseKernel
    {
        Real src_inv_h_, src_inv_h_squared_, src_inv_h_cubed_;
        Real inv_h_, inv_h_squared_, inv_h_cubed_, inv_h_fourth_, inv_h_fifth_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        SmoothingKernel(Neighbor<SPHAdaptation, SPHAdaptation> &encloser);
        inline Vecd nablaW_ij(UnsignedInt i, UnsignedInt j) const;
        inline Real W_ij(UnsignedInt i, UnsignedInt j) const { return W(vec_r_ij(i, j)); };
        inline Real dW_ij(UnsignedInt i, UnsignedInt j) const { return dW(vec_r_ij(i, j)); };
        Real W0(UnsignedInt i, const Vec2d &zero) const;
        Real W0(UnsignedInt i, const Vec3d &zero) const;

        Real W(const Vec2d &displacement) const;
        Real W(const Vec3d &displacement) const;
        Real dW(const Vec2d &displacement) const;
        Real dW(const Vec3d &displacement) const;
        Real d2W(const Vec2d &displacement) const;
        Real d2W(const Vec3d &displacement) const;
        inline Real CutOffRadius() const { return kernel_size_ / inv_h_; }
    };
    typedef SmoothingKernel NeighborKernel;

    class NeighborCriterion
    {
        Vecd *src_pos_, *tar_pos_;
        Real kernel_size_squared_, inv_h_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        bool operator()(UnsignedInt j, UnsignedInt i) const; // Note the reverse order of indices
    };

    class ReverseNeighborCriterion
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReverseNeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser){};

        inline bool operator()(UnsignedInt i, UnsignedInt j) const { return true; };
    };

    class CutOff
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        CutOff(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        inline Vecd operator()(UnsignedInt i) const { return cut_off_; };

      private:
        Vecd cut_off_;
    };

  protected:
    Real src_inv_h_, inv_h_;
};

template <class SourceAdaptationType, class TargetAdaptationType>
class Neighbor<SourceAdaptationType, TargetAdaptationType> : public Neighbor<Base>
{
    using SourceSmoothingLengthRatio = typename SourceAdaptationType::SmoothingLengthRatioType;
    using TargetSmoothingLengthRatio = typename TargetAdaptationType::SmoothingLengthRatioType;

  public:
    template <typename SourceIdentifier, typename TargetIdentifier>
    Neighbor(SourceIdentifier &source_identifier, TargetIdentifier &target_identifier,
             DiscreteVariable<Vecd> *dv_src_pos, DiscreteVariable<Vecd> *dv_tar_pos);

    class SmoothingKernel : public Neighbor<Base>::SmoothingKernel
    {
        using BaseKernel = Neighbor<Base>::SmoothingKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        inline Vecd nablaW_ij(UnsignedInt i, UnsignedInt j) const;
        inline Real W_ij(UnsignedInt i, UnsignedInt j) const;
        inline Real dW_ij(UnsignedInt i, UnsignedInt j) const;
        Real W0(UnsignedInt i, const Vec2d &zero) const;
        Real W0(UnsignedInt i, const Vec3d &zero) const;

      protected:
        Real src_inv_h_ref_, tar_inv_h_ref_;
        SourceSmoothingLengthRatio src_h_ratio_;
        TargetSmoothingLengthRatio tar_h_ratio_;

        std::tuple<Vecd, Real, bool> getTransformedMeasure(UnsignedInt i, UnsignedInt j) const;
        Real W(const Vec2d &disp_transform, Real inv_h) const;
        Real W(const Vec3d &disp_transform, Real inv_h) const;
        Real dW(const Vec2d &disp_transform, Real inv_h) const;
        Real dW(const Vec3d &disp_transform, Real inv_h) const;
        Real d2W(const Vec2d &disp_transform, Real inv_h) const;
        Real d2W(const Vec3d &disp_transform, Real inv_h) const;
    };
    typedef SmoothingKernel NeighborKernel;

    class NeighborCriterion
    {
        Vecd *src_pos_;
        Vecd *tar_pos_;
        Real kernel_size_squared_, src_inv_h_ref_;
        SourceSmoothingLengthRatio src_h_ratio_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        bool operator()(UnsignedInt j, UnsignedInt i) const; // Note the reverse order of indices
    };

    class ReverseNeighborCriterion
    {
        Vecd *src_pos_;
        Vecd *tar_pos_;
        Real kernel_size_squared_, tar_inv_h_ref_;
        TargetSmoothingLengthRatio tar_h_ratio_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        ReverseNeighborCriterion(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        bool operator()(UnsignedInt i, UnsignedInt j) const;
    };

    class CutOff
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        CutOff(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Vecd operator()(UnsignedInt i) const;

      private:
        Real kernel_size_, src_inv_h_ref_;
        SourceSmoothingLengthRatio src_h_ratio_;
    };

  protected:
    Real src_inv_h_ref_, tar_inv_h_ref_;
    Real src_inv_h_min_, tar_inv_h_min_;
    SourceAdaptationType &src_adaptation_;
    TargetAdaptationType &tar_adaptation_;
};
} // namespace SPH
#endif // NEIGHBOR_METHOD_H
