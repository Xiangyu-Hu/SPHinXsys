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
 * @file 	line_particle_generator.h
 * @brief generate particle for linear structures.
 * @author	Xipeng Lyu and Xiangyu Hu
 */

#ifndef LINE_PARTICLE_GENERATOR_H
#define LINE_PARTICLE_GENERATOR_H

#include "base_particle_generator.h"

namespace SPH
{

template <> // Generate a smooth line-shape using particles
class ParticleGenerator<LinearParticles> : public ParticleGenerator<BaseParticles>
{
    StdVec<Vecd> line_normal_;    /**< line normal */
    StdVec<Real> line_thickness_; /**< line thickness */
    StdVec<Vecd> line_binormal_;  /**< line binormal */
    StdVec<Real> line_width_;     /**< line width */

  public:
    explicit ParticleGenerator(SPHBody &sph_body, LinearParticles &linear_particles);
    virtual ~ParticleGenerator() {};

  protected:
    LinearParticles &linear_particles_;
    virtual void addLineProperties(
        const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width);
    virtual void initializeParticleVariables() override;
    virtual void initializeParticleVariablesFromReload() override;
};
} // namespace SPH
#endif // LINE_PARTICLE_GENERATOR_H
