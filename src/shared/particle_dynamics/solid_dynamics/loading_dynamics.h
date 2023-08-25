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
 * @file 	loading_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef LOADING_DYNAMICS_H
#define LOADING_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "general_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
namespace solid_dynamics
{
//----------------------------------------------------------------------
//		for general solid dynamics
//----------------------------------------------------------------------
typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
typedef DataDelegateInner<SolidParticles> SolidDataInner;

/**@class ImposeExternalForce
 * @brief impose external force on a solid body part
 * by add extra acceleration
 */
class ImposeExternalForce : public LocalDynamics, public SolidDataSimple
{
  public:
    ImposeExternalForce(SPHBody &sph_body);
    virtual ~ImposeExternalForce(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos0_, &vel_;
    /**
     * @brief acceleration will be specified by the application
     */
    virtual Vecd getAcceleration(Vecd &pos) = 0;
};

/**
 * @class SpringDamperConstraintParticleWise
 * @brief Exerts spring force and damping force in the form of acceleration to each particle.
 * The spring force is calculated based on the difference from the particle's initial position.
 * The damping force is calculated based on the particle's current velocity.
 * Only for 3D applications
 */
class SpringDamperConstraintParticleWise : public LocalDynamics, public SolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_, &pos0_, &vel_, &acc_prior_;
    Vecd stiffness_;
    Vecd damping_coeff_; // damping component parallel to the spring force component

    virtual Vecd getSpringForce(size_t index_i, Vecd &disp);
    virtual Vecd getDampingForce(size_t index_i);

  public:
    SpringDamperConstraintParticleWise(SPHBody &sph_body, Vecd stiffness, Real damping_ratio = 0.05);

    void update(size_t index_i, Real dt = 0.0);
};
/**
 * @class SpringNormalOnSurfaceParticles
 * @brief Exerts spring force force on the surface in normal direction in the form of acceleration to each particle.
 * The input stiffness should be defined in Pa/m. The stiffness is scaled by the surface area of the particle to get N/m
 * The force is applied to all the surface particles that can be seen (outer_surface = false)
 * or cannot be seen (outer_surface = true) from the source point.
 * Can be used for outer or inner surface of a shell structure ofr example.
 * The spring force is calculated based on the difference from the particle's initial position.
 * Only for 3D applications
 * Only for uniform surface particle size.
 */
class SpringNormalOnSurfaceParticles : public LocalDynamics, public SolidDataSimple
{
  public:
    SpringNormalOnSurfaceParticles(SPHBody &sph_body, bool outer_surface,
                                   Vecd source_point, Real stiffness, Real damping_ratio = 0.05);

    StdLargeVec<bool> &GetApplySpringForceToParticle() { return apply_spring_force_to_particle_; }
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_, &pos0_, &n_, &n0_, &vel_, &acc_prior_;
    StdLargeVec<Real> &mass_;
    Real stiffness_;
    Real damping_coeff_; // damping component parallel to the spring force component
    StdLargeVec<bool> apply_spring_force_to_particle_;

    virtual Vecd getSpringForce(size_t index_i, Vecd disp);
    virtual Vecd getDampingForce(size_t index_i);
};
/**
 * @class SpringOnSurfaceParticles
 * @brief Exerts spring force force on the surface in the form of acceleration to each particle.
 * The input stiffness should be defined in Pa/m. The stiffness is scaled by the surface area of the particle to get N/m
 * The force is applied to all the surface particles.
 * The spring force is calculated based on the difference from the particle's initial position.
 * Only for 3D applications
 * BodyPartByParticle define the ody part that the spring is applied to.
 * Only for uniform surface particle size.
 */
class SpringOnSurfaceParticles : public LocalDynamics, public SolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_, &pos0_, &vel_, &acc_prior_;
    StdLargeVec<Real> &mass_;
    Real stiffness_;
    Real damping_coeff_; // damping component parallel to the spring force component
    StdLargeVec<bool> apply_spring_force_to_particle_;

  public:
    SpringOnSurfaceParticles(SPHBody &sph_body, Real stiffness, Real damping_ratio = 0.05);

    const StdLargeVec<bool> &GetApplySpringForceToParticle() { return apply_spring_force_to_particle_; }
    void update(size_t index_i, Real dt = 0.0);
};
/**
 * @class AccelerationForBodyPartInBoundingBox
 * @brief Adds acceleration to the part of the body that's inside a bounding box
 */
class AccelerationForBodyPartInBoundingBox : public LocalDynamics, public SolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_, &acc_prior_;
    BoundingBox bounding_box_;
    Vecd acceleration_;

  public:
    AccelerationForBodyPartInBoundingBox(SPHBody &sph_body, BoundingBox &bounding_box, Vecd acceleration);
    virtual ~AccelerationForBodyPartInBoundingBox(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ForceInBodyRegion
 * @brief ForceInBodyRegion, distributes the force vector as acceleration among the particles in a given body part
 */
class ForceInBodyRegion : public BaseLocalDynamics<BodyPartByParticle>, public SolidDataSimple
{
  public:
    ForceInBodyRegion(BodyPartByParticle &body_part, Vecd force, Real end_time);
    virtual ~ForceInBodyRegion(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos0_, &acc_prior_;
    Vecd acceleration_;
    Real end_time_;
};

/**
 * @class SurfacePressureFromSource
 * @brief SurfacePressureFromSource, applies pressure on the surface particles coming from a source point
 */
class SurfacePressureFromSource : public BaseLocalDynamics<BodyPartByParticle>, public SolidDataSimple
{
  public:
    SurfacePressureFromSource(BodyPartByParticle &body_part,
                              Vecd source_point, StdVec<std::array<Real, 2>> pressure_over_time);

    StdLargeVec<bool> &GetApplyPressureToParticle() { return apply_pressure_to_particle_; }
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos0_, &n_, &acc_prior_;
    StdLargeVec<Real> &mass_;
    StdVec<std::array<Real, 2>> pressure_over_time_;
    StdLargeVec<bool> apply_pressure_to_particle_;
    Real getPressure();
};
} // namespace solid_dynamics
} // namespace SPH
#endif // LOADING_DYNAMICS_H