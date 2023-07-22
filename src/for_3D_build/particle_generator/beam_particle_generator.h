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
 * @file 	beam_particle_generator.h
 * @brief generate particle for linear structures.
 * @author	Xipeng Lyu
 */

#ifndef BEAM_PARTICLE_GENERATOR_H
#define BEAM_PARTICLE_GENERATOR_H

#include "base_data_package.h"
#include "base_particle_generator.h"
#include "large_data_containers.h"
#include "sph_data_containers.h"

namespace SPH
{

class SPHBody;
class BaseParticles;
class IOEnvironment;

/**
 * @class LineParticleGenerator
 * @brief Generate volumetric particles by initialize extra Line variables.
 */
class LineParticleGenerator : public ParticleGenerator
{
  public:
    explicit LineParticleGenerator(SPHBody &sph_body);
    virtual ~LineParticleGenerator(){};

  protected:
    StdLargeVec<Vecd> &n_;         /**< line normal */
    StdLargeVec<Real> &thickness_; /**< line thickness */
    StdLargeVec<Vecd> &b_n_;       /**< line binormal */
    StdLargeVec<Real> &width_;     /**< line width */

    /** Initialize line particle. */
    virtual void initializeLineProperties(const Vecd &line_normal, const Vecd &line_binormal, Real thickness, Real width);
};

} // namespace SPH
#endif // BEAM_PARTICLE_GENERATOR_H
