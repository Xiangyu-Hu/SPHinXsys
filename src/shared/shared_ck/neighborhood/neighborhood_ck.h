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

#include "kernel_wenland_c2_ck.h"
#include "neighborhood.h"

namespace SPH
{
template <typename... T>
class Neighbor;

template <>
class Neighbor<>
{
  public:
    template <class ExecutionPolicy>
    Neighbor(const ExecutionPolicy &ex_policy, SPHAdaptation *sph_adaptation, DiscreteVariable<Vecd> *dv_pos);

    template <class ExecutionPolicy>
    Neighbor(const ExecutionPolicy &ex_policy, SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
             DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_target_pos);

    inline Vecd vec_r_ij(size_t i, size_t j) const { return source_pos_[i] - target_pos_[j]; };
    inline Real W_ij(size_t i, size_t j) const { return kernel_.W(vec_r_ij(i, j)); }
    inline Real dW_ij(size_t i, size_t j) const { return kernel_.dW(vec_r_ij(i, j)); }

    inline Vecd e_ij(size_t i, size_t j) const
    {
        Vecd displacement = vec_r_ij(i, j);
        return displacement / (displacement.norm() + TinyReal);
    }

  protected:
    KernelWendlandC2CK kernel_;
    Vecd *source_pos_;
    Vecd *target_pos_;
};

class NeighborList
{
  public:
    template <class ExecutionPolicy>
    NeighborList(const ExecutionPolicy &ex_policy,
                 DiscreteVariable<UnsignedInt> *dv_neighbor_index,
                 DiscreteVariable<UnsignedInt> *dv_particle_offset);

  protected:
    UnsignedInt *neighbor_index_;
    UnsignedInt *particle_offset_;
    inline UnsignedInt FirstNeighbor(UnsignedInt i) { return particle_offset_[i]; };
    inline UnsignedInt LastNeighbor(UnsignedInt i) { return particle_offset_[i + 1]; };
};
} // namespace SPH
#endif // NEIGHBORHOOD_CK_H