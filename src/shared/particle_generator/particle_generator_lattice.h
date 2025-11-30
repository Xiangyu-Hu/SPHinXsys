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
 * @file particle_generator_lattice.h
 * @brief The lattice generator generates particles
 * at lattice position by check whether the position is contained by a SPH body.
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_GENERATOR_LATTICE_H
#define PARTICLE_GENERATOR_LATTICE_H

#include "base_particle_generator.h"
#include "particle_generator_reserve.h"
namespace SPH
{

class Shape;
class AdaptiveByShape;
class SurfaceParticles;

template <> // Base class for generating particles from lattice positions
class GeneratingMethod<Lattice>
{
  public:
    explicit GeneratingMethod(SPHBody &sph_body);
    virtual ~GeneratingMethod() {};

  protected:
    Real lattice_spacing_;      /**< Initial particle spacing. */
    BoundingBoxd domain_bounds_; /**< Domain bounds. */
    Shape &initial_shape_;      /**< Geometry shape for body. */
};

template <>
class ParticleGenerator<BaseParticles, Lattice>
    : public ParticleGenerator<BaseParticles>, public GeneratingMethod<Lattice>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles);
    virtual ~ParticleGenerator() {};
    virtual void prepareGeometricData() override;
};

template <> // For generating particles with adaptive resolution from lattice positions
class ParticleGenerator<BaseParticles, Lattice, AdaptiveByShape> : public ParticleGenerator<BaseParticles, Lattice>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, Shape &target_shape);
    explicit ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles);
    virtual ~ParticleGenerator() {};

  protected:
    Shape &target_shape_;
    AdaptiveByShape *particle_adaptation_;
    virtual void addPositionAndVolumetricMeasure(const Vecd &position, Real volume) override;
};

template <> // For generating surface particles from lattice positions using reduced order approach
class ParticleGenerator<SurfaceParticles, Lattice>
    : public ParticleGenerator<SurfaceParticles>, public GeneratingMethod<Lattice>
{
  public:
    ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles, Real thickness);
    virtual ~ParticleGenerator() {};
    virtual void prepareGeometricData() override;

  protected:
    Real total_volume_;                  /**< Total volume of body calculated from level set. */
    Real thickness_;                     /**< Global average thickness. */
    Real particle_spacing_;              /**< Particle spacing. */
    Real avg_particle_volume_;           /**< Average particle volume. */
    size_t all_cells_;                   /**< Number of cells enclosed by the volume. */
    size_t planned_number_of_particles_; /**< Number of particles in planned manner. */
};
} // namespace SPH
#endif // PARTICLE_GENERATOR_LATTICE_H
