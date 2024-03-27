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
 * @file 	general_life_time_dynamics.h
 * @brief 	Classes on life time related events.
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "base_general_dynamics.h"

namespace SPH
{
/**
 * @class BaseLifeTimeDynamics
 * @brief Base class for particle life time events.
 */
class BaseLifeTimeDynamics : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    explicit BaseLifeTimeDynamics(SPHBody &sph_body);
    virtual ~BaseLifeTimeDynamics(){};

  protected:
    ParticleSplitAndMerge &particle_split_merge_;
    Real inv_rho0_;
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Real> &mass_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class BaseSplitDynamics
 * @brief Base class for particle split.
 */
template <typename SplitParameters>
class BaseSplitDynamics : public BaseLifeTimeDynamics
{
  public:
    BaseSplitDynamics(SPHBody &sph_body, size_t body_buffer_width)
        : BaseLifeTimeDynamics(sph_body)
    {
        if (particles_->total_real_particles_ == particles_->real_particles_bound_)
        {
            std::cout << "EmitterInflowBoundaryCondition constructor: \n"
                      << "No buffer particles have been generated! "
                      << "Please check the particle generator."
                      << "\n";
            exit(1);
        }
    };
    virtual ~BaseSplitDynamics(){};

  protected:
    virtual bool checkSplit(size_t index_i) = 0;
    virtual SplitParameters execFirstSplit(size_t index_i) = 0;
    virtual void execOtherSplit(size_t index_i, const SplitParameters &split_parameters) = 0;
};

/**
 * @class RefinementInPrescribedRegion
 * @brief particle split in prescribed region.
 */
class RefinementInPrescribedRegion : public BaseSplitDynamics<Vecd>
{
  public:
    RefinementInPrescribedRegion(SPHBody &sph_body, size_t body_buffer_width, Shape &refinement_region);
    virtual ~RefinementInPrescribedRegion(){};
    virtual void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);

  protected:
    std::mutex mutex_split_; /**< mutex exclusion for memory conflict */
    BoundingBox refinement_region_bounds_;
    std::random_device random_device_;
    std::mt19937 random_seed_;
    std::normal_distribution<Real> normal_distribution_;

    virtual bool checkSplit(size_t index_i) override;
    virtual Vecd execFirstSplit(size_t index_i) override;
    virtual void execOtherSplit(size_t index_i, const Vecd &split_shift) override;
    virtual bool checkLocation(const BoundingBox &refinement_region_bounds, Vecd position, Real volume);
};
} // namespace SPH
