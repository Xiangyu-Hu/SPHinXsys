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
template <class NeighborMethod>
class Neighbor
{
    using SmoothingKernel = typename NeighborMethod::SmoothingKernel;

  public:
    template <class ExecutionPolicy>
    Neighbor(const ExecutionPolicy &ex_policy, NeighborMethod &neighbor_method)
        : smoothing_kernel_(ex_policy, neighbor_method){};

    inline Vecd vec_r_ij(size_t i, size_t j) const { return smoothing_kernel_.vec_r_ij(i, j); };
    inline Vecd vec_r_i(UnsignedInt i) const { return smoothing_kernel_.vec_r_i[i]; };
    inline Vecd vec_r_j(UnsignedInt j) const { return smoothing_kernel_.vec_r_j[j]; };
    inline Real W_ij(size_t i, size_t j) const { return smoothing_kernel_.W_ij(i, j); };
    inline Real dW_ij(size_t i, size_t j) const { return smoothing_kernel_.dW_ij(i, j); };
    inline Vecd e_ij(size_t i, size_t j) const { return smoothing_kernel_.e_ij(i, j); };
    inline Real W(const Vecd &displacement) const { return smoothing_kernel_.W(displacement); };

  protected:
    SmoothingKernel smoothing_kernel_;
};

template <class NeighborMethod>
class NeighborCriterion
{
    using CriterionKernel = typename NeighborMethod::CriterionKernel;

  public:
    template <class ExecutionPolicy>
    NeighborCriterion(const ExecutionPolicy &ex_policy, NeighborMethod &neighbor_method)
        : criterion_kernel_(ex_policy, neighbor_method){};

    inline bool operator()(UnsignedInt target_index, UnsignedInt source_index) const
    {
        return criterion_kernel_(source_index, target_index);
    };

  protected:
    CriterionKernel criterion_kernel_;
};

} // namespace SPH
#endif // NEIGHBORHOOD_CK_H