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
 * @file 	neighborhood_ck.h
 * @brief 	There are the classes for particle neighborhood.
 * It saves the information for carrying out inter-particle (or pair) interaction,
 * and also considered as the local configuration of the particles.
 * @author	Xiangyu Hu
 */

#ifndef NEIGHBORHOOD_CK_H
#define NEIGHBORHOOD_CK_H

#include "neighbor_method.h"

namespace SPH
{
template <class NeighborMethodType>
class Neighbor
{
  public:
    template <class SourceIdentifier, class TargetIdentifier>
    Neighbor(SourceIdentifier &source_identifier, TargetIdentifier &contact_identifier,
             DiscreteVariable<Vecd> *dv_source_pos, DiscreteVariable<Vecd> *dv_target_pos)
        : neighbor_method_(source_identifier, contact_identifier),
          dv_source_pos_(dv_source_pos), dv_target_pos_(dv_target_pos){};

    class NeighborKernel
    {
        using SmoothingKernel = typename NeighborMethodType::SmoothingKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        NeighborKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : smoothing_kernel_(ex_policy, encloser.neighbor_method_),
              source_pos_(encloser.dv_source_pos_->DelegatedData(ex_policy)),
              target_pos_(encloser.dv_target_pos_->DelegatedData(ex_policy)){};

        inline Vecd vec_r_ij(UnsignedInt i, UnsignedInt j) const { return source_pos_[i] - target_pos_[j]; };

        inline Vecd e_ij(UnsignedInt i, UnsignedInt j) const
        {
            Vecd displacement = vec_r_ij(i, j);
            return displacement / (displacement.norm() + TinyReal);
        };

        inline Real W_ij(UnsignedInt i, UnsignedInt j) const { return smoothing_kernel_.W(vec_r_ij(i, j)); };
        inline Real dW_ij(UnsignedInt i, UnsignedInt j) const { return smoothing_kernel_.dW(vec_r_ij(i, j)); };
        inline Real W(const Vecd &displacement) const { return smoothing_kernel_.W(displacement); };

      protected:
        SmoothingKernel smoothing_kernel_;
        Vecd *source_pos_;
        Vecd *target_pos_;
    };

    typedef typename NeighborMethodType::CriterionKernel CriterionKernelType;

    template <class ExecutionPolicy>
    CriterionKernelType createCriterionKernel(const ExecutionPolicy &ex_policy)
    {
        return CriterionKernelType(ex_policy, neighbor_method_, dv_source_pos_, dv_target_pos_);
    }

    class SearchDepth
    {
        using SearchDepthKernel = typename NeighborMethodType::SearchDepthKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SearchDepth(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : search_depth_kernel_(ex_policy, encloser.neighbor_method_){};

        inline int operator()(UnsignedInt i) const { return search_depth_kernel_(i); };

      protected:
        SearchDepthKernel search_depth_kernel_;
    };

    class SmoothingRatio
    {
        using SmoothingRatioKernel = typename NeighborMethodType::SmoothingRatioKernel;

      public:
        template <class ExecutionPolicy, class EncloserType>
        SmoothingRatio(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : smoothing_ratio_kernel_(ex_policy, encloser.neighbor_method_){};

        inline Real operator()(UnsignedInt i) const { return smoothing_ratio_kernel_(i); };

      protected:
        SmoothingRatioKernel smoothing_ratio_kernel_;
    };

  protected:
    NeighborMethodType neighbor_method_; /**< The neighbor method for the neighborhood. */
    DiscreteVariable<Vecd> *dv_source_pos_;
    DiscreteVariable<Vecd> *dv_target_pos_;
};
} // namespace SPH
#endif // NEIGHBORHOOD_CK_H