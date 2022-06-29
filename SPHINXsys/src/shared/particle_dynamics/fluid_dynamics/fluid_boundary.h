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
 * @file 	fluid_boundary.h
 * @brief 	Here, we define the boundary condition classes for fluid dynamics.
 * @details The boundary conditions very often based on different types of buffers.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_BOUNDARY_H
#define FLUID_BOUNDARY_H

#include "fluid_dynamics_inner.h"

namespace SPH
{
    namespace fluid_dynamics
    {
        /**
         * @class FlowRelaxationBuffer
         * @brief Flow buffer in which the particles relaxes to a given target velocity profile.
         * This technique will be used for applying several boundary conditions,
         * such as freestream, inflow, damping boundary conditions.
         */
        class FlowRelaxationBuffer : public PartDynamicsByCell, public FluidDataSimple
        {
        public:
            FlowRelaxationBuffer(FluidBody &fluid_body, BodyPartByCell &body_part);
            virtual ~FlowRelaxationBuffer(){};

        protected:
            StdLargeVec<Vecd> &pos_n_, &vel_n_;
            /** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
            Real relaxation_rate_;

            /** inflow profile to be defined in applications,
             * argument parameters and return value are in frame (local) coordinate */
            virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) = 0;
            virtual void Update(size_t index_i, Real dt = 0.0) override;
        };

        /**
         * @class InflowBoundaryCondition
         * @brief Inflow boundary condition which imposes directly to a given velocity profile.
         */
        class InflowBoundaryCondition : public FlowRelaxationBuffer
        {
        public:
            InflowBoundaryCondition(FluidBody &fluid_body, BodyAlignedBoxByCell &aligned_box_part);
            virtual ~InflowBoundaryCondition(){};

        protected:
            Transformd &transform_;
            Vecd halfsize_;

            virtual void Update(size_t index_i, Real dt = 0.0) override;
        };

        /**
         * @class DampingBoundaryCondition
         * @brief damping boundary condition which relaxes
         * the particles to zero velocity profile.
         * TODO: one can using aligned box shape and generalize the damping factor along 
         * one axis direction.
         */
        class DampingBoundaryCondition : public PartDynamicsByCell, public FluidDataSimple
        {
        public:
            DampingBoundaryCondition(FluidBody &fluid_body, BodyRegionByCell &body_part);
            virtual ~DampingBoundaryCondition(){};

        protected:
            StdLargeVec<Vecd> &pos_n_, &vel_n_;
            /** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
            Real strength_;
            BoundingBox damping_zone_bounds_;
            virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
        };

        /**
         * @class EmitterInflowCondition
         * @brief Inflow boundary condition imposed on an emitter, in which pressure and density profile are imposed too.
         * The body part region is required to have parallel lower- and upper-bound surfaces.
         */
        class EmitterInflowCondition : public PartSimpleDynamicsByParticle, public FluidDataSimple
        {
        public:
            explicit EmitterInflowCondition(FluidBody &fluid_body, BodyAlignedBoxByParticle &aligned_box_part);
            virtual ~EmitterInflowCondition(){};

        protected:
            StdLargeVec<Vecd> &pos_n_, &vel_n_;
            StdLargeVec<Real> &rho_n_, &p_;
            /** inflow pressure condition */
            Real inflow_pressure_;
            Real rho0_;
            AlignedBoxShape &aligned_box_;
            Transformd &updated_transform_, old_transform_;

            /** no transform by default */
            virtual void updateTransform(){};
            virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) = 0;
            virtual void setupDynamics(Real dt = 0.0) override { updateTransform(); };
            virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override;
        };

        /**
         * @class EmitterInflowInjecting
         * @brief Inject particles into the computational domain.
         */
        class EmitterInflowInjecting : public PartSimpleDynamicsByParticle, public FluidDataSimple
        {
        public:
            explicit EmitterInflowInjecting(FluidBody &fluid_body, BodyAlignedBoxByParticle &aligned_box_part,
                                            size_t body_buffer_width, int axis_direction, bool positive);
            virtual ~EmitterInflowInjecting(){};

            /** This class is only implemented in sequential due to memory conflicts. */
            virtual void parallel_exec(Real dt = 0.0) override { exec(); };
            AlignedBoxShape &getBodyPartByParticle(){};

        protected:
            StdLargeVec<Vecd> &pos_n_;
            StdLargeVec<Real> &rho_n_, &p_;
            const int axis_; /**< the axis direction for bounding*/
            size_t body_buffer_width_;
            AlignedBoxShape &aligned_box_;

            virtual void checkLowerBound(size_t unsorted_index_i, Real dt = 0.0);
            virtual void checkUpperBound(size_t unsorted_index_i, Real dt = 0.0);
            ParticleFunctor checking_bound_;

            virtual void Update(size_t unsorted_index_i, Real dt = 0.0) override
            {
                checking_bound_(unsorted_index_i, dt);
            };
        };

        /**
         * @class StaticConfinementDensity
         * @brief static confinement condition for density summation
         */
        class StaticConfinementDensity : public PartDynamicsByCell, public FluidDataSimple
        {
        public:
            StaticConfinementDensity(FluidBody &fluid_body, NearShapeSurface &near_surface);
            virtual ~StaticConfinementDensity(){};

        protected:
            Real rho0_, inv_sigma0_;
            StdLargeVec<Real> &mass_, &rho_sum_;
            StdLargeVec<Vecd> &pos_n_;
            LevelSetShape *level_set_shape_;

            virtual void Update(size_t index_i, Real dt = 0.0) override;
        };

        /**
         * @class StaticConfinementPressureRelaxation
         * @brief static confinement condition for pressure relaxation
         */
        class StaticConfinementPressureRelaxation : public PartDynamicsByCell, public FluidDataSimple
        {
        public:
            StaticConfinementPressureRelaxation(FluidBody &fluid_body, NearShapeSurface &near_surface);
            virtual ~StaticConfinementPressureRelaxation(){};

        protected:
            StdLargeVec<Real> &rho_n_, &p_;
            StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_;
            LevelSetShape *level_set_shape_;
            AcousticRiemannSolver riemann_solver_;

            virtual void Update(size_t index_i, Real dt = 0.0) override;
        };

        /**
         * @class StaticConfinementDensityRelaxation
         * @brief static confinement condition for density relaxation
         */
        class StaticConfinementDensityRelaxation : public PartDynamicsByCell, public FluidDataSimple
        {
        public:
            StaticConfinementDensityRelaxation(FluidBody &fluid_body, NearShapeSurface &near_surface);
            virtual ~StaticConfinementDensityRelaxation(){};

        protected:
            StdLargeVec<Real> &rho_n_, &p_, &drho_dt_;
            StdLargeVec<Vecd> &pos_n_, &vel_n_;
            LevelSetShape *level_set_shape_;
            AcousticRiemannSolver riemann_solver_;

            virtual void Update(size_t index_i, Real dt = 0.0) override;
        };

        /**
         * @class StaticConfinement
         * @brief Static confined boundary condition for complex structures.
         */
        class StaticConfinement
        {
        public:
            StaticConfinementDensity density_summation_;
            StaticConfinementPressureRelaxation pressure_relaxation_;
            StaticConfinementDensityRelaxation density_relaxation_;

            StaticConfinement(FluidBody &fluid_body, NearShapeSurface &near_surface);
            virtual ~StaticConfinement(){};
        };

    }
}
#endif // FLUID_BOUNDARY_H