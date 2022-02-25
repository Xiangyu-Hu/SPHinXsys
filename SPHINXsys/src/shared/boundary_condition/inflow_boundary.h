/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file 	inflow_boundary.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the inlet region.
 * @details 	These classes allow the boundary conditon to be constructed in any direction.
 * @author	Huiqiang Yue
 */

#ifndef INFLOW_BOUNDARY_H
#define INFLOW_BOUNDARY_H

#include "boundary_face.h"

namespace SPH
{
    /**
     * @brief A body part with the collection of particles containing a boundary face.
     */
    class BodyRegionByParticleWithFace : public BodyRegionWithFace
    {
    public:
        IndexVector body_part_particles_; /**< Collection particle in this body part. */
        BodyRegionByParticleWithFace(RealBody& real_body, SegmentFace& segment_face, Real scale = 1.0);

    protected:
        void tagParticles();

    protected:
        BaseParticles* base_particles_;
    };

    class PartDynamicsByParticleWithFace : public ParticleDynamics<void>
    {
    public:
        PartDynamicsByParticleWithFace(RealBody& real_body, BodyRegionByParticleWithFace& body_part)
            : ParticleDynamics<void>(real_body),
            body_part_particles_(body_part.body_part_particles_) {};
        virtual ~PartDynamicsByParticleWithFace() {}

        virtual void exec(Real dt = 0.0) override;
        virtual void parallel_exec(Real dt = 0.0) override;

    protected:
        IndexVector& body_part_particles_;
        ParticleFunctor particle_functor_;
    };

    class PartSimpleDynamicsByParticleWithFace : public PartDynamicsByParticleWithFace
    {
    public:
        PartSimpleDynamicsByParticleWithFace(RealBody& real_body, BodyRegionByParticleWithFace& body_part);
        virtual ~PartSimpleDynamicsByParticleWithFace() {}

    protected:
        virtual void Update(size_t index_i, Real dt = 0.0) = 0;
    };

    /**
     * @brief Injecting particle to computational domain at inlet.
     */
    class InflowInjectingWithFace : public PartSimpleDynamicsByParticleWithFace,
        public DataDelegateSimple<FluidBody, FluidParticles, Fluid>
    {
    public:
        explicit InflowInjectingWithFace(FluidBody& fluid_body, BodyRegionByParticleWithFace& body_part, size_t body_buffer_width);
        virtual ~InflowInjectingWithFace() {}
        /** This class is only implemented in sequential due to memory conflicts. */
        virtual void parallel_exec(Real dt = 0.0) override { exec(); };

    protected:
        BodyRegionByParticleWithFace& body_part_;
        StdLargeVec<Vecd>& pos_n_;
        StdLargeVec<Real>& rho_n_, & p_;
        Real periodic_translation_;
        virtual void checking_bound_(size_t unsorted_index_i, Real dt = 0.0);
        virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override
        {
            checking_bound_(unsorted_index_i, dt);
        }
    };

    /**
     * @brief Inflow boundary condition which impose target velocity profile.
     *  The body part region is no more limited to having parallel coordinate axis surfaces.
     */
    class InflowConditionWithFace : public PartSimpleDynamicsByCellsWithFace,
        public DataDelegateSimple<FluidBody, FluidParticles, Fluid>
    {
    public:
        explicit InflowConditionWithFace(FluidBody& fluid_body, BodyRegionByCellsWithFace& body_part);
        virtual ~InflowConditionWithFace() {}

    protected:
        BodyRegionByCellsWithFace& body_part_;
        StdLargeVec<Vecd>& vel_n_, & pos_n_;
    };

    class VelocityInflowConditionWithFace : public InflowConditionWithFace
    {
    public:
        explicit VelocityInflowConditionWithFace(FluidBody& fluid_body, BodyRegionByCellsWithFace& body_part);
        virtual ~VelocityInflowConditionWithFace() {}

    protected:
        /** inflow velocity profile to be defined in applications */
        virtual Vecd getTargetVelocity(Vecd& position, Vecd& velocity);
        virtual Vecd defineVelocityProfile(Vecd& position, Vecd& velocity) = 0;
        virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override;
    };

} // namespace SPH


#endif // !INFLOW_BOUNDARY_H


