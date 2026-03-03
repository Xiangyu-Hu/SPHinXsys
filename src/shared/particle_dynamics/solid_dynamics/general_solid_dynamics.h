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
 * @file 	general_solid_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_SOLID_DYNAMICS_H
#define GENERAL_SOLID_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_general_dynamics.h"
#include "base_kernel.h"
#include "elastic_solid.h"

namespace SPH
{
namespace solid_dynamics
{
/**
 * @class DistributingPointForces
 * @brief Distribute a series of point forces to its contact shell bodies.
 */
class DistributingPointForces : public LocalDynamics
{
  protected:
    std::vector<Vecd> point_forces_, reference_positions_, time_dependent_point_forces_;
    Real time_to_full_external_force_;
    Real particle_spacing_ref_, h_spacing_ratio_;
    Vecd *pos_, *force_prior_;
    Real *thickness_;
    std::vector<Real *> weight_;
    std::vector<Real> sum_of_weight_;
    Real *physical_time_;

    void getWeight();

  public:
    DistributingPointForces(SPHBody &sph_body, std::vector<Vecd> point_forces,
                            std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                            Real particle_spacing_ref, Real h_spacing_ratio = 1.6);
    virtual ~DistributingPointForces() {};

    virtual void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace solid_dynamics
} // namespace SPH
#endif // GENERAL_SOLID_DYNAMICS_H