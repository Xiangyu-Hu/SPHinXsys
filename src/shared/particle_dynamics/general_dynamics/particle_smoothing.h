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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file particle_smoothing.h
 * @brief Methods on smoothing discrete values with its neighbors
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_SMOOTHING_H
#define PARTICLE_SMOOTHING_H

#include "base_general_dynamics.h"

namespace SPH
{
/**
 * @class ParticleSmoothing
 * @brief computing smoothed variable field by averaging with neighbors
 */
template <typename VariableType>
class ParticleSmoothing : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit ParticleSmoothing(BaseInnerRelation &inner_relation, const std::string &variable_name);
    virtual ~ParticleSmoothing(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    const Real W0_;
    VariableType *smoothed_, *temp_;
};

/**
 * @class ParticleSnapshotAverage
 * @brief TBD.
 */
template <typename VariableType>
class ParticleSnapshotAverage : public LocalDynamics
{
  public:
    explicit ParticleSnapshotAverage(SPHBody &sph_body, const std::string &variable_name);
    virtual ~ParticleSnapshotAverage(){};
    virtual void setupDynamics(Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    size_t number_of_snapshot_ = 0;
    VariableType *target_variable_, *averaged_variable_;
};
} // namespace SPH
#endif // PARTICLE_SMOOTHING_H
