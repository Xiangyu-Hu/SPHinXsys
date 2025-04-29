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

#include "kernel_tabulated_ck.h"
#include "neighborhood.h"

namespace SPH
{
template <typename... T>
class Neighbor;

template <>
class Neighbor<Base>
{
  public:
    template <class ExecutionPolicy>
    Neighbor(const ExecutionPolicy &ex_policy,
             SPHAdaptation *sph_adaptation, SPHAdaptation *target_adaptation,
             DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_target_pos);

    KernelTabulatedCK &getKernel() { return kernel_; }
    inline Vecd vec_r_ij(UnsignedInt i, UnsignedInt j) const { return source_pos_[i] - target_pos_[j]; };

    inline Vecd e_ij(UnsignedInt i, UnsignedInt j) const
    {
        Vecd displacement = vec_r_ij(i, j);
        return displacement / (displacement.norm() + TinyReal);
    }

  protected:
    KernelTabulatedCK kernel_;
    Real kernel_size_square_;
    Vecd *source_pos_;
    Vecd *target_pos_;
};

template <class NeighborMethod>
class Neighbor<NeighborMethod> : public Neighbor<Base>
{
    using ScalingFactor = typename NeighborMethod::ComputingKernel;

  public:
    template <class ExecutionPolicy>
    Neighbor(const ExecutionPolicy &ex_policy,
             SPHAdaptation *sph_adaptation, SPHAdaptation *target_adaptation,
             DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_target_pos, NeighborMethod &smoothing_length);
    inline Vecd normalizedVecRij(UnsignedInt i, UnsignedInt j) const { return inv_h_(i, j) * (source_pos_[i] - target_pos_[j]); }
    inline Real W_ij(UnsignedInt i, UnsignedInt j) const { return kernel_.W(normalizedVecRij(i, j)); }
    inline Real dW_ij(UnsignedInt i, UnsignedInt j) const { return inv_h_(i, j) * kernel_.dW(normalizedVecRij(i, j)); }

    class NeighborCriterion
    {
      public:
        NeighborCriterion(Neighbor<NeighborMethod> &neighbor);
        bool operator()(UnsignedInt target_index, UnsignedInt source_index) const
        {
            Vecd normalized_displacement = inv_h_(source_index, target_index) *
                                           (source_pos_[source_index] - target_pos_[target_index]);
            return normalized_displacement.squaredNorm() < kernel_size_square_;
        };

      protected:
        Vecd *source_pos_;
        Vecd *target_pos_;
        ScalingFactor inv_h_;
        Real kernel_size_square_;
    };

  protected:
    ScalingFactor inv_h_;
};
} // namespace SPH
#endif // NEIGHBORHOOD_CK_H