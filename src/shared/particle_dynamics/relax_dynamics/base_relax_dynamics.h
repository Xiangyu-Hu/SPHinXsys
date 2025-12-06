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
 * @file 	base_relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_RELAX_DYNAMICS_H
#define BASE_RELAX_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
class LevelSetCorrection
{
};

namespace relax_dynamics
{
/**
 * @class RandomizeParticlePosition
 * @brief Randomize the initial particle position
 */
class RandomizeParticlePosition
    : public LocalDynamics
{
  protected:
    Vecd *pos_;
    Real randomize_scale_;

  public:
    explicit RandomizeParticlePosition(SPHBody &sph_body);
    virtual ~RandomizeParticlePosition() {};

    void update(size_t index_i, Real dt = 0.0);
};

} // namespace relax_dynamics
} // namespace SPH
#endif // BASE_RELAX_DYNAMICS_H
