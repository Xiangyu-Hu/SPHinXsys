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
 * @file solid_to_shell_coupling.h
 * @brief Constraints for solid to shell coupling.
 * @author	Xiangyu Hu, Weiyi Kong (Virtonomy GmbH)
 */

#ifndef SOLID_TO_SHELL_CONSTRAINT_H
#define SOLID_TO_SHELL_CONSTRAINT_H

#include "force_prior.h"
#include "general_constraint.h"

namespace SPH
{
namespace solid_dynamics
{
/**@class ComputeWeight
 * @brief Compute the total weight for each solid particle from the coupled shell particles.
 * This weight is used for interpolation of velocity and force
 */
template <class DynamicsIdentifier>
class TotalWeightComputation : public BaseLocalDynamics<DynamicsIdentifier>,
                               public DataDelegateContact
{
  private:
    Real *total_weight_; // \sum_i W_ij * Vol_j for j of all contact bodies
    StdVec<Real *> contact_Vol_;

  public:
    TotalWeightComputation(DynamicsIdentifier &identifier, BaseContactRelation &contact_relation);
    void update(size_t index_i, Real dt = 0.0);
};

/**@class ConsistentMapping
 * @brief Constrain the variable by the interpolation from coupled particle velocities with consistent mapping.
 * Usually used for velocity.
 * When there is more than one contact body, interpolate from all contact bodies.
 */
template <class DynamicsIdentifier, typename DataType>
class ConsistentMapping : public BaseLocalDynamics<DynamicsIdentifier>,
                          public DataDelegateContact
{
  private:
    DataType *interpolated_quantities_;
    Real *total_weight_;
    StdVec<Real *> contact_Vol_;
    StdVec<DataType *> contact_data_;

  public:
    ConsistentMapping(DynamicsIdentifier &identifier,
                      BaseContactRelation &contact_relation,
                      const std::string &variable_name,          // variable name of the source body, e.g., "Velocity"
                      const std::string &contact_variable_name); // variable name of the contact body, e.g., "Velocity"
    void update(size_t index_i, Real dt = 0.0);
};

/**@class ConservativeMapping
 * @brief Constrain the variable by the interpolation from coupled particle with conservative mapping.
 * Usually used for force. Note that for prior force, the force needs to be registered as a new prior force.
 * When there is more than one contact body, interpolate force from each contact body and sum up.
 */
template <class DynamicsIdentifier, typename DataType>
class ConservativeMapping : public BaseLocalDynamics<DynamicsIdentifier>,
                            public DataDelegateContact
{
  private:
    DataType *interpolated_quantities_;
    Real *Vol_;
    StdVec<Real *> contact_total_weight_;
    StdVec<DataType *> contact_data_;

  public:
    ConservativeMapping(DynamicsIdentifier &identifier,
                        BaseContactRelation &contact_relation,
                        const std::string &variable_name,          // variable name of the source body, e.g., "CouplingForce"
                        const std::string &contact_variable_name); // variable name of the contact body, e.g., "Force"
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace solid_dynamics
} // namespace SPH

#endif // SOLID_TO_SHELL_CONSTRAINT_H