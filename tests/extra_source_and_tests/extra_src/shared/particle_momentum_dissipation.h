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
 * @file 	particle_dynamics_dissipation.h
 * @brief 	Here are the classes for damping the magnitude of
 * 			any variables.
 * 			Note that, currently, these classes works only in single resolution.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_MOMENTUM_DISSIPATION_H
#define PARTICLE_MOMENTUM_DISSIPATION_H

#include "all_particle_dynamics.h"
#include "porous_media_dynamics.h"

namespace SPH
{
namespace multi_species_continuum
{
/**
 * @class PorousMediaDampingPairwiseInner
 * @brief A quantity damping by a pairwise splitting scheme
 * this method modifies the quantity directly
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */

template <typename VariableType>
class PorousMediaDampingPairwiseInner : public LocalDynamics, public DataDelegateInner
{
  public:
    PorousMediaDampingPairwiseInner(BaseInnerRelation &inner_relation, const std::string &variable_name, Real eta);
    virtual ~PorousMediaDampingPairwiseInner() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_, *mass_;

    VariableType *variable_;
    Real eta_; /**< damping coefficient */
};
} // namespace multi_species_continuum
} // namespace SPH
#endif // PARTICLE_MOMENTUM_DISSIPATION_H