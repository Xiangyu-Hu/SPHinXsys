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
    template <typename DataType>
    inline Real W_ij(const DataType &r_ij) const
    {
        return kernel_.W(r_ij);
    }
    template <typename DataType>
    inline Real dW_ij(const DataType &r_ij) const
    {
        return kernel_.dW(r_ij);
    }
    inline Vecd r_ij(size_t index_i, size_t index_j) const { return source_pos_[index_i] - target_pos_[index_j]; };
    inline Vecd e_ij(size_t index_i, size_t index_j) const
    {
        Vecd displacement = r_ij(index_i, index_j);
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
    inline UnsignedInt FirstNeighbor(UnsignedInt index_i) { return particle_offset_[index_i]; };
    inline UnsignedInt LastNeighbor(UnsignedInt index_i) { return particle_offset_[index_i] + 1; };

  private:
    UnsignedInt *particle_offset_;
};
} // namespace SPH
#endif // NEIGHBORHOOD_CK_H