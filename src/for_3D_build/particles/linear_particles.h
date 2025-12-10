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
 * @file linear_particles.h
 * @brief This is the derived class of shell particles.
 * @author Xipeng Lyu and Xiangyu Hu
 */

#ifndef LINEAR_PARTICLES_H
#define LINEAR_PARTICLES_H

#include "surface_particles.h"

namespace SPH
{

class LinearParticles : public SurfaceParticles
{
  public:
    LinearParticles(SPHBody &sph_body, BaseMaterial *base_material);
    virtual ~LinearParticles() {};

    Vecd *b_n_; /**< binormal direction */
    Real *width_;

    /** get particle volume. */
    virtual Real ParticleVolume(size_t index_i) override { return Vol_[index_i] * thickness_[index_i] * width_[index_i]; }
    void registerLineProperties(StdVec<Vecd> &b_n, StdVec<Real> &width);
    void registerLinePropertiesFromReload();
    virtual void registerTransformationMatrix() override;
};

} // namespace SPH
#endif // SOLID_PARTICLES_H
