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
 * @brief Buffer particles are saved behind real particles.
 * They may be switched from real particles or switch to real particles.
 * Ghost particles whose states are updated according to
 * boundary condition if their indices are included in the neighbor particle list.
 * The ghost particles are saved behind the buffer particles in the form of one or more ghost bounds.
 * @details Note that, due to present arrangement of memories,
 * the generation of buffer particles must be done before the generation of ghost particles.
 * Therefore, type traits are used to check the order of the particle generation.
 * @author Xiangyu Hu
 */

#ifndef PARTICLE_GENERATOR_RESERVE_H
#define PARTICLE_GENERATOR_RESERVE_H

#include "base_particle_generator.h"
#include "particle_reserve.h"
namespace SPH
{
template <class T, class = void>
struct has_ghost_particles : std::false_type
{
};

template <class T>
struct has_ghost_particles<T, std::void_t<decltype(&T::reserveGhostParticles)>> : std::true_type
{
};

template <class ParticlesType, class ReserveSizeEstimator, typename... OtherParameters>
class ParticleGenerator<ParticlesType, RealParticleReserve<ReserveSizeEstimator>, OtherParameters...>
    : public ParticleGenerator<ParticlesType, OtherParameters...>
{
  public:
    template <typename... Args>
    ParticleGenerator(SPHBody &sph_body, ParticlesType &particles,
                      RealParticleReserve<ReserveSizeEstimator> &particle_reserve, Args &&...args)
        : ParticleGenerator<ParticlesType, OtherParameters...>(sph_body, particles, std::forward<Args>(args)...),
          particle_reserve_(particle_reserve)
    {
        static_assert(!has_ghost_particles<ParticleGenerator<ParticlesType, OtherParameters...>>::value,
                      "ParticleGenerator: GhostReservation is not allowed ahead of BufferReservation.");
    };
    virtual ~ParticleGenerator(){};

    virtual void setAllParticleBounds() override
    {
        ParticleGenerator<ParticlesType, OtherParameters...>::setAllParticleBounds();
        particle_reserve_.reserveRealParticles(this->base_particles_, this->particle_spacing_ref_);
    };

  private:
    RealParticleReserve<ReserveSizeEstimator> &particle_reserve_;
};

template <class ParticlesType, typename GhostParameter, typename... OtherParameters>
class ParticleGenerator<ParticlesType, Ghost<GhostParameter>, OtherParameters...>
    : public ParticleGenerator<ParticlesType, OtherParameters...>
{
  public:
    template <typename... Args>
    ParticleGenerator(SPHBody &sph_body, ParticlesType &particles,
                      Ghost<GhostParameter> &ghost_boundary, Args &&...args)
        : ParticleGenerator<ParticlesType, OtherParameters...>(sph_body, particles, std::forward<Args>(args)...),
          ghost_boundary_(ghost_boundary){};
    virtual ~ParticleGenerator(){};

    virtual void setAllParticleBounds() override
    {
        ParticleGenerator<ParticlesType, OtherParameters...>::setAllParticleBounds();
        ghost_boundary_.reserveGhostParticles(this->base_particles_, this->particle_spacing_ref_);
    };

  protected:
    Ghost<GhostParameter> &ghost_boundary_;
};
} // namespace SPH
#endif // PARTICLE_GENERATOR_RESERVE_H
