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
 * @file surface_particles.h
 * @brief This is the derived class of base particles.
 * @author Chi Zhang, Dong Wu and Xiangyu Hu
 */

#ifndef SURFACE_PARTICLES_H
#define SURFACE_PARTICLES_H

#include "base_particles.hpp"

namespace SPH
{
/**
 * @class SurfaceParticles
 * @brief A group of particles with shell particle data.
 */
class SurfaceParticles : public BaseParticles
{
  public:
    SurfaceParticles(SPHBody &sph_body, BaseMaterial *base_material);
    virtual ~SurfaceParticles() {};

    Vecd *n_;                      /**< normal direction */
    Real *thickness_;              /**< shell thickness */
    Matd *transformation_matrix0_; /**< initial transformation matrix from global to local coordinates */

    void registerSurfaceProperties(StdVec<Vecd> &n, StdVec<Real> &thickness);
    void registerSurfacePropertiesFromReload();
    virtual Real ParticleVolume(size_t index_i) override { return Vol_[index_i] * thickness_[index_i]; }
    virtual void registerTransformationMatrix();
    virtual void initializeBasicParticleVariables() override;
};

} // namespace SPH
#endif // SURFACE_PARTICLES_H
