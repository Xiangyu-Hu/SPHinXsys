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
 * @file base_particle_generator.h
 * @brief This is the base class of particle generator, which generates particles
 * with given positions and volumes. The direct generator simply generate
 * particle with given position and volume.
 * The particle generators naming are using template partial specialization,
 * with three keywords. The first indicates volume (default), surface and line particles,
 * the second generating methods, and the third the control parameters,
 * such as adaptive for adaptive resolution.
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
//---------------------------------------------------------------------------
// Geometric types of particles. The default type is volume metric particles.
//---------------------------------------------------------------------------
class Surface;
class ThickSurface; // Surface thickness equal or larger than the particle spacing
class Observer;
class Reload;

template <typename... Parameters>
class ParticleGenerator;

template <typename... Parameters>
class GeneratingMethod;

template <> // default volume metric particle generator
class ParticleGenerator<Base>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body);
    virtual ~ParticleGenerator(){};
    virtual void prepareGeometricData() = 0;
    void generateParticlesWithGeometricVariables();

  protected:
    BaseParticles &base_particles_;
    Real particle_spacing_ref_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<size_t> &unsorted_id_;
    virtual void prepareParticlePosition(const Vecd &position);
    virtual void preparePositionAndVolumetricMeasure(const Vecd &position, Real volumetric_measure);
};

template <> // generate surface particles
class ParticleGenerator<Surface> : public ParticleGenerator<Base>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body);
    virtual ~ParticleGenerator(){};

  protected:
    StdLargeVec<Vecd> &n_;         /**< surface normal */
    StdLargeVec<Real> &thickness_; /**< surface thickness */
    virtual void prepareSurfaceProperties(const Vecd &surface_normal, Real thickness);
};

template <> // generate observer particles
class ParticleGenerator<Observer> : public ParticleGenerator<Base>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body)
        : ParticleGenerator<Base>(sph_body){};
    ParticleGenerator(SPHBody &sph_body, const StdVec<Vecd> &positions)
        : ParticleGenerator<Base>(sph_body), positions_(positions){};
    virtual ~ParticleGenerator(){};
    virtual void prepareGeometricData() override;

  protected:
    StdVec<Vecd> positions_;
};

template <> // generate particles by reloading dynamically relaxed particles
class ParticleGenerator<Reload> : public ParticleGenerator<Base>
{
    BaseMaterial &base_material_;
    std::string file_path_;

  public:
    ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name);
    virtual ~ParticleGenerator(){};
    virtual void prepareGeometricData() override;
};

} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_H
