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
 * @file 	base_particle_generator.h
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLE_GENERATOR_H
#define BASE_PARTICLE_GENERATOR_H

#include "base_data_package.h"
#include "large_data_containers.h"
#include "sph_data_containers.h"

namespace SPH
{

class SPHBody;
class BaseParticles;
class IOEnvironment;

/**
 * @class BaseParticleGenerator
 * @brief Abstract base particle generator.
 */
class BaseParticleGenerator
{
  public:
    explicit BaseParticleGenerator(SPHBody &sph_body);
    virtual ~BaseParticleGenerator(){};
    /** Initialize geometric parameters. */
    virtual void initializeGeometricVariables() = 0;
    virtual void generateParticlesWithBasicVariables();

  protected:
    BaseParticles &base_particles_;
    BaseMaterial &base_material_;
    StdLargeVec<Vecd> &pos_;           /**< current position */
    StdLargeVec<size_t> &unsorted_id_; /**< original particle ids */
    /** Initialize particle position. */
    virtual void initializePosition(const Vecd &position);
};

/**
 * @class ParticleGenerator
 * @brief Generate volumetric particles by initialize position and volume.
 */
class ParticleGenerator : public BaseParticleGenerator
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body);
    virtual ~ParticleGenerator(){};

  protected:
    StdLargeVec<Real> &Vol_; /**< particle volume */
    /** Initialize particle position and measured volume. */
    virtual void initializePositionAndVolumetricMeasure(const Vecd &position, Real volumetric_measure);
};

/**
 * @class SurfaceParticleGenerator
 * @brief Generate volumetric particles by initialize extra surface variables.
 */
class SurfaceParticleGenerator : public ParticleGenerator
{
  public:
    explicit SurfaceParticleGenerator(SPHBody &sph_body);
    virtual ~SurfaceParticleGenerator(){};

  protected:
    StdLargeVec<Vecd> &n_;         /**< surface normal */
    StdLargeVec<Real> &thickness_; /**< surface thickness */
    /** Initialize surface particle. */
    virtual void initializeSurfaceProperties(const Vecd &surface_normal, Real thickness);
};


/**
 * @class ObserverParticleGenerator
 * @brief Generate particle directly from position-and-volume data.
 * @details The values of positions will be given in the derived class.
 */
class ObserverParticleGenerator : public ParticleGenerator
{
  public:
    explicit ObserverParticleGenerator(SPHBody &sph_body)
        : ParticleGenerator(sph_body){};
    ObserverParticleGenerator(SPHBody &sph_body, const StdVec<Vecd> &positions)
        : ParticleGenerator(sph_body), positions_(positions){};
    virtual ~ObserverParticleGenerator(){};
    /** Initialize geometrical variable for observe particles. */
    virtual void initializeGeometricVariables() override;

  protected:
    StdVec<Vecd> positions_;
};

/**
 * @class ParticleGeneratorReload
 * @brief Generate particle by reloading particle position and volume.
 */
class ParticleGeneratorReload : public ParticleGenerator
{
    std::string file_path_;

  public:
    ParticleGeneratorReload(SPHBody &sph_body, IOEnvironment &io_environment, const std::string &reload_body_name);
    virtual ~ParticleGeneratorReload(){};
    /** Initialize geometrical variable for reload particles. */
    virtual void initializeGeometricVariables() override;
    virtual void generateParticlesWithBasicVariables() override;
};
} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_H
