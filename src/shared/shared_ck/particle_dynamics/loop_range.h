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
 * @file 	loop_range.h
 * @brief 	Here, we define loop ranges for parallel computing.
 * @author	Xiangyu Hu
 */

#ifndef LOOP_RANGE_H
#define LOOP_RANGE_H

#include "base_body.h"
#include "base_body_part.h"
#include "base_particles.hpp"
#include "reduce_functors.h"

namespace SPH
{
template <typename...>
class LoopRangeCK;

template <class ExecutionPolicy>
class LoopRangeCK<ExecutionPolicy, SPHBody>
{
  public:
    LoopRangeCK(SPHBody &sph_body)
        : loop_bound_(sph_body.getBaseParticles().svTotalRealParticles()->DelegatedData(ExecutionPolicy{})) {};
    LoopRangeCK(SingularVariable<UnsignedInt> *sv_total_particles)
        : loop_bound_(sv_total_particles->DelegatedData(ExecutionPolicy{})) {};
    template <class UnaryFunc>
    void computeUnit(const UnaryFunc &f, UnsignedInt i) const { f(i); };
    template <class ReturnType, class BinaryFunc, class UnaryFunc>
    ReturnType computeUnit(ReturnType temp, const BinaryFunc &bf, const UnaryFunc &uf, UnsignedInt i) const { return uf(i); };
    UnsignedInt LoopBound() const { return *loop_bound_; };

  protected:
    UnsignedInt *loop_bound_;
};

template <class ExecutionPolicy>
class LoopRangeCK<ExecutionPolicy, BodyPartByParticle>
{
  public:
    LoopRangeCK(BodyPartByParticle &body_part)
        : particle_list_(body_part.dvParticleList()->DelegatedData(ExecutionPolicy{})),
          loop_bound_(body_part.svRangeSize()->DelegatedData(ExecutionPolicy{})) {};
    template <class UnaryFunc>
    void computeUnit(const UnaryFunc &f, UnsignedInt i) const { f(particle_list_[i]); };
    template <class ReturnType, class BinaryFunc, class UnaryFunc>
    ReturnType computeUnit(ReturnType temp, const BinaryFunc &bf, const UnaryFunc &uf, UnsignedInt i) const
    {
        return uf(particle_list_[i]);
    };
    UnsignedInt LoopBound() const { return *loop_bound_; };

  protected:
    UnsignedInt *particle_list_;
    UnsignedInt *loop_bound_;
};

template <class ExecutionPolicy>
class LoopRangeCK<ExecutionPolicy, BodyPartByCell>
{
  public:
    LoopRangeCK(BodyPartByCell &body_part)
        : cell_list_(body_part.dvCellList()->DelegatedData(ExecutionPolicy{})),
          loop_bound_(body_part.svRangeSize()->DelegatedData(ExecutionPolicy{})),
          particle_index_(body_part.dvParticleIndex()->DelegatedData(ExecutionPolicy{})),
          cell_offset_(body_part.dvCellOffset()->DelegatedData(ExecutionPolicy{})) {};
    template <class UnaryFunc>
    void computeUnit(const UnaryFunc &uf, UnsignedInt i) const
    {
        UnsignedInt cell_index = cell_list_[i];
        for (size_t k = cell_offset_[cell_index]; k != cell_offset_[cell_index + 1]; ++k)
        {
            uf(particle_index_[k]);
        }
    };

    template <class ReturnType, class BinaryFunc, class UnaryFunc>
    ReturnType computeUnit(ReturnType temp, const BinaryFunc &bf, const UnaryFunc &uf, UnsignedInt i) const
    {
        UnsignedInt cell_index = cell_list_[i];
        for (size_t k = cell_offset_[cell_index]; k != cell_offset_[cell_index + 1]; ++k)
        {
            temp = bf(temp, uf(particle_index_[k]));
        }
        return temp;
    };
    UnsignedInt LoopBound() const { return *loop_bound_; };

  protected:
    UnsignedInt *cell_list_;
    UnsignedInt *loop_bound_;
    UnsignedInt *particle_index_;
    UnsignedInt *cell_offset_;
};
} // namespace SPH
#endif // LOOP_RANGE_H
