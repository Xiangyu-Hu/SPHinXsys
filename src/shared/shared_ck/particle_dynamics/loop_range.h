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
 * @file 	loop_range.h
 * @brief 	Here, we define loop ranges for parallel computing.
 * @author	Xiangyu Hu
 */

#ifndef LOOP_RANGE_H
#define LOOP_RANGE_H

#include "base_body.h"
#include "base_body_part.h"
#include "base_particles.hpp"

namespace SPH
{
template <typename...>
class LoopRangeCK;

template <class ExecutionPolicy>
class LoopRangeCK<ExecutionPolicy, SPHBody>
{
  public:
    LoopRangeCK(SPHBody &sph_body)
        : loop_bound_(sph_body.getBaseParticles().svTotalRealParticles()->DelegatedData(ExecutionPolicy{})){};
    virtual ~LoopRangeCK(){};
    UnsignedInt LoopBound() { return *loop_bound_; };

  protected:
    UnsignedInt *loop_bound_;
};

template <class ExecutionPolicy>
class LoopRangeCK<ExecutionPolicy, BodyPartByParticle>
{
  public:
    LoopRangeCK(BodyPartByParticle &body_part_by_particle)
        : particle_indexes_(body_part_by_particle.dvParticleIndexes()->DelegatedDataField(ExecutionPolicy{})),
          loop_bound_(body_part_by_particle.svRangeSize()->DelegatedData(ExecutionPolicy{})){};

    UnsignedInt *ParticleIndexes() { return particle_indexes_; };
    UnsignedInt LoopBound() { return *loop_bound_; };

  protected:
    UnsignedInt *particle_indexes_;
    UnsignedInt *loop_bound_;
};
} // namespace SPH
#endif // LOOP_RANGE_H