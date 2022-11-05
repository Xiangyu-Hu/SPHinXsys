/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	fluid_boundary.h
 * @brief 	Here, we define the boundary condition classes for fluid dynamics.
 * @details The boundary conditions very often based on different types of buffers.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_BOUNDARY_H
#define FLUID_BOUNDARY_H

#include "fluid_dynamics_inner.h"

#include <mutex>
namespace SPH
{
    namespace fluid_dynamics
    {
        /**
         * @class BaseFlowBoundaryCondition
         * @brief Base class for all boundary conditions.
         */
        class BaseFlowBoundaryCondition : public LocalDynamics, public FluidDataSimple
        {
        public:
            BaseFlowBoundaryCondition(BodyPartByCell &body_part)
                : LocalDynamics(body_part.getSPHBody()), FluidDataSimple(sph_body_),
                  rho_(particles_->rho_), p_(particles_->p_),
                  pos_(particles_->pos_), vel_(particles_->vel_){};
            virtual ~BaseFlowBoundaryCondition(){};

        protected:
            StdLargeVec<Real> &rho_, &p_;
            StdLargeVec<Vecd> &pos_, &vel_;
        };

        /**
         * @class FlowVelocityBuffer
         * @brief Flow buffer in which the particle velocity relaxes to a given target profile.
         * This technique will be used for applying several boundary conditions,
         * such as freestream, inflow, damping boundary conditions.
         */
        template <typename TargetVelocity>
        class FlowVelocityBuffer : public BaseFlowBoundaryCondition
        {
        public:
            FlowVelocityBuffer(BodyPartByCell &body_part, Real relaxation_rate = 0.3)
                : BaseFlowBoundaryCondition(body_part),
                  relaxation_rate_(relaxation_rate), target_velocity(){};
            virtual ~FlowVelocityBuffer(){};

            void update(size_t index_i, Real dt = 0.0)
            {
                vel_[index_i] += relaxation_rate_ * (target_velocity(pos_[index_i], vel_[index_i]) - vel_[index_i]);
            };

        protected:
            /** default value is 0.3 suggests reaching target profile in several time steps */
            Real relaxation_rate_;
            TargetVelocity target_velocity;
        };

        /**
         * @class InflowVelocityCondition
         * @brief Inflow boundary condition which imposes directly to a given velocity profile.
         */
        template <typename TargetVelocity>
        class InflowVelocityCondition : public BaseFlowBoundaryCondition
        {
        public:
            InflowVelocityCondition(BodyAlignedBoxByCell &aligned_box_part)
                : BaseFlowBoundaryCondition(aligned_box_part),
                  transform_(aligned_box_part.aligned_box_.getTransform()),
                  halfsize_(aligned_box_part.aligned_box_.HalfSize()),
                  target_velocity(){};
            virtual ~InflowVelocityCondition(){};

            void update(size_t index_i, Real dt = 0.0)
            {
                Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
                Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
                Vecd prescribed_velocity =
                    transform_.xformFrameVecToBase(target_velocity(frame_position, frame_velocity));
                vel_[index_i] = prescribed_velocity;
            };

        protected:
            Transformd &transform_;
            Vecd halfsize_;
            TargetVelocity target_velocity;
        };

        /**
         * @class FarFieldVelocityCorrection
         * @brief this function is applied to freestream flows
         * @brief modify the velocity of free surface particles with far-field velocity
         */
        template <typename TargetVelocity>
        class FarFieldVelocityCorrection : public LocalDynamics, public FluidDataSimple
        {
        protected:
            Real u_ref_, t_ref_, rho_ref_;
            StdLargeVec<Real> &rho_sum;
            StdLargeVec<Vecd> &pos_, &vel_;
            StdLargeVec<int> &surface_indicator_;
            TargetVelocity target_velocity;

        public:
            explicit FarFieldVelocityCorrection(SPHBody &sph_body)
                : LocalDynamics(sph_body), FluidDataSimple(sph_body),
                  u_ref_(1.0), t_ref_(2.0), rho_ref_(particles_->fluid_.ReferenceDensity()),
                  rho_sum(particles_->rho_sum_), pos_(particles_->pos_), vel_(particles_->vel_),
                  surface_indicator_(*particles_->getVariableByName<int>("SurfaceIndicator")),
                  target_velocity(){};
            virtual ~FarFieldVelocityCorrection(){};

            void update(size_t index_i, Real dt = 0.0)
            {
                if (surface_indicator_[index_i] == 1)
                {
                    Vecd u_freestream = target_velocity(pos_[index_i], vel_[index_i]);
                    vel_[index_i] = u_freestream + SMIN(rho_sum[index_i], rho_ref_) * (vel_[index_i] - u_freestream) / rho_ref_;
                }
            };
        };

        /**
         * @class DampingBoundaryCondition
         * @brief damping boundary condition which relaxes
         * the particles to zero velocity profile.
         * TODO: one can using aligned box shape and generalize the damping factor along
         * one axis direction.
         */
        class DampingBoundaryCondition : public BaseFlowBoundaryCondition
        {
        public:
            DampingBoundaryCondition(BodyRegionByCell &body_part);
            virtual ~DampingBoundaryCondition(){};
            void update(size_t index_particle_i, Real dt = 0.0);

        protected:
            /** default value is 0.1 suggests reaching  target inflow velocity in about 10 time steps */
            Real strength_;
            BoundingBox damping_zone_bounds_;
        };

        /**
         * @class EmitterInflowCondition
         * @brief Inflow boundary condition imposed on an emitter, in which pressure and density profile are imposed too.
         * The body part region is required to have parallel lower- and upper-bound surfaces.
         */
        class EmitterInflowCondition : public LocalDynamics, public FluidDataSimple
        {
        public:
            explicit EmitterInflowCondition(BodyAlignedBoxByParticle &aligned_box_part);
            virtual ~EmitterInflowCondition(){};

            virtual void setupDynamics(Real dt = 0.0) override { updateTransform(); };
            void update(size_t unsorted_index_i, Real dt = 0.0);

        protected:
            Fluid &fluid_;
            StdLargeVec<Vecd> &pos_, &vel_, &acc_;
            StdLargeVec<Real> &rho_, &p_, &drho_dt_;
            /** inflow pressure condition */
            Real inflow_pressure_;
            Real rho0_;
            AlignedBoxShape &aligned_box_;
            Transformd &updated_transform_, old_transform_;

            /** no transform by default */
            virtual void updateTransform(){};
            virtual Vecd getTargetVelocity(Vecd &position, Vecd &velocity) = 0;
        };

        /**
         * @class EmitterInflowInjection
         * @brief Inject particles into the computational domain.
         * Note that the axis is at the local coordinate and upper bound direction is
         * the local positive direction.
         */
        class EmitterInflowInjection : public LocalDynamics, public FluidDataSimple
        {
        public:
            EmitterInflowInjection(BodyAlignedBoxByParticle &aligned_box_part,
                                   size_t body_buffer_width, int axis);
            virtual ~EmitterInflowInjection(){};

            void update(size_t unsorted_index_i, Real dt = 0.0);

        protected:
            std::mutex mutex_switch_to_real_; /**< mutex exclusion for memory conflict */
            Fluid &fluid_;
            StdLargeVec<Vecd> &pos_;
            StdLargeVec<Real> &rho_, &p_;
            const int axis_; /**< the axis direction for bounding*/
            AlignedBoxShape &aligned_box_;
        };

        /**
         * @class DisposerOutflowDeletion
         * @brief Delete particles who ruing out the computational domain.
         */
        class DisposerOutflowDeletion : public LocalDynamics, public FluidDataSimple
        {
        public:
            DisposerOutflowDeletion(BodyAlignedBoxByCell &aligned_box_part, int axis);
            virtual ~DisposerOutflowDeletion(){};

            void update(size_t index_i, Real dt = 0.0);

        protected:
            std::mutex mutex_switch_to_buffer_; /**< mutex exclusion for memory conflict */
            StdLargeVec<Vecd> &pos_;
            const int axis_; /**< the axis direction for bounding*/
            AlignedBoxShape &aligned_box_;
        };

        /**
         * @class StaticConfinementDensity
         * @brief static confinement condition for density summation
         */
        class StaticConfinementDensity : public LocalDynamics, public FluidDataSimple
        {
        public:
            StaticConfinementDensity(NearShapeSurface &near_surface);
            virtual ~StaticConfinementDensity(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real rho0_, inv_sigma0_;
            StdLargeVec<Real> &mass_, &rho_sum_;
            StdLargeVec<Vecd> &pos_;
            LevelSetShape *level_set_shape_;
        };

        /**
         * @class StaticConfinementIntegration1stHalf
         * @brief static confinement condition for pressure relaxation
         */
        class StaticConfinementIntegration1stHalf : public LocalDynamics, public FluidDataSimple
        {
        public:
            StaticConfinementIntegration1stHalf(NearShapeSurface &near_surface);
            virtual ~StaticConfinementIntegration1stHalf(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Fluid &fluid_;
            StdLargeVec<Real> &rho_, &p_;
            StdLargeVec<Vecd> &pos_, &vel_, &acc_;
            LevelSetShape *level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };

        /**
         * @class StaticConfinementIntegration2ndHalf
         * @brief static confinement condition for density relaxation
         */
        class StaticConfinementIntegration2ndHalf : public LocalDynamics, public FluidDataSimple
        {
        public:
            StaticConfinementIntegration2ndHalf(NearShapeSurface &near_surface);
            virtual ~StaticConfinementIntegration2ndHalf(){};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Fluid &fluid_;
            StdLargeVec<Real> &rho_, &p_, &drho_dt_;
            StdLargeVec<Vecd> &pos_, &vel_;
            LevelSetShape *level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };

        /**
         * @class StaticConfinement
         * @brief Static confined boundary condition for complex structures.
         */
        class StaticConfinement
        {
        public:
            SimpleDynamics<StaticConfinementDensity, NearShapeSurface> density_summation_;
            SimpleDynamics<StaticConfinementIntegration1stHalf, NearShapeSurface> pressure_relaxation_;
            SimpleDynamics<StaticConfinementIntegration2ndHalf, NearShapeSurface> density_relaxation_;

            StaticConfinement(NearShapeSurface &near_surface);
            virtual ~StaticConfinement(){};
        };

    }
}
#endif // FLUID_BOUNDARY_H