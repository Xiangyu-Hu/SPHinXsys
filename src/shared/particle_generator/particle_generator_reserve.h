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
 * @file particle_generator_reserve.h
 * @brief TBD
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_GENERATOR_RESERVE_H
#define PARTICLE_GENERATOR_RESERVE_H

#include "base_particle_generator.h"

namespace SPH
{
class BufferReservation;

template <class T, class = void>
struct has_ghost_particles : std::false_type
{
};

template <class T>
struct has_ghost_particles<T, std::void_t<decltype(&T::reserveGhostParticles)>> : std::true_type
{
};

template <typename... OtherParameters>
class ParticleGenerator<BufferReservation, OtherParameters...>
    : public ParticleGenerator<OtherParameters...>
{
  public:
    template <typename... Args>
    ParticleGenerator(SPHBody &sph_body, Real buffer_size_factor, Args &&...args)
        : ParticleGenerator<OtherParameters...>(sph_body, std::forward<Args>(args)...)
    {
        static_assert(!has_ghost_particles<ParticleGenerator<OtherParameters...>>::value,
                      "ParticleGenerator: GhostReservation is not allowed ahead of BufferReservation.");
        buffer_size_ = std::ceil(Real(this->base_particles_.total_real_particles_) * buffer_size_factor);
    };
    virtual ~ParticleGenerator(){};

    virtual void generateParticlesWithBasicVariables() override
    {
        ParticleGenerator<OtherParameters...>::generateParticlesWithBasicVariables();
        reserveBufferParticles();
    };

  protected:
    void reserveBufferParticles()
    {
        Real this->base_particles_.addBufferParticles(buffer_size);
    };

  private:
    Real buffer_size_;
};

template <typename BoundaryType, typename... OtherParameters>
class ParticleGenerator<Ghost<BoundaryType>, OtherParameters...>
    : public ParticleGenerator<OtherParameters...>
{
  public:
    template <typename... Args>
    ParticleGenerator(SPHBody &sph_body, Ghost<BoundaryType> &ghost_boundary, Args &&...args)
        : ParticleGenerator<OtherParameters...>(sph_body, std::forward<Args>(args)...){};
    virtual ~ParticleGenerator(){};

    virtual void generateParticlesWithBasicVariables() override
    {
        ParticleGenerator<OtherParameters...>::generateParticlesWithBasicVariables();
        reserveGhostParticles();
    };

  protected:
    Ghost<BoundaryType> &ghost_boundary_;

    void reserveGhostParticles() :
    {
        ghost_boundary_.reserveGhostParticle(this->base_particles_);
    };
};
} // namespace SPH
#endif // PARTICLE_GENERATOR_RESERVE_H
