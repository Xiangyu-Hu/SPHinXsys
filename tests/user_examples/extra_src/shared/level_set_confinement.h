/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	static_confinement.h
 * @brief 	Here, we define the static confinement boundary condition classes for fluid dynamics.
 * @details     This boundary condition is based on Level-set filed.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef LEVEL_SET_COFINEMENT_H
#define LEVEL_SET_COFINEMENT_H


#include "shape_confinement.h"
#include "body_part_by_cell_tracing.h"
#include "base_fluid_dynamics.h"
#include "general_constraint.h"
#include "riemann_solver.h"
#include <mutex>

namespace SPH
{
    namespace fluid_dynamics
    {

        /**
         * @class StaticConfinementTransportVelocity
         * @brief static confinement condition for transport velocity
         */
        class StaticConfinementTransportVelocity : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementTransportVelocity(NearShapeSurface& near_surface, Real coefficient = 0.2);
            virtual ~StaticConfinementTransportVelocity() {};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<int>& surface_indicator_;
            const Real coefficient_;
            Real smoothing_length_sqr_;
            LevelSetShape* level_set_shape_;
            StdLargeVec<Vecd> &transport_acc_;
        };

        /**
         * @class StaticConfinementViscousAcceleration
         * @brief static confinement condition for viscous acceleration
         */
        class StaticConfinementViscousAcceleration : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementViscousAcceleration(NearShapeSurface& near_surface);
            virtual ~StaticConfinementViscousAcceleration() {};
            void update(size_t index_i, Real dt = 0.0);
            //StdLargeVec<Vecd> &getForceFromFluid() { return force_from_fluid_; };
        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<Real>& rho_, &mass_;
            StdLargeVec<Vecd>& vel_, &force_prior_;
            Real mu_;
            LevelSetShape* level_set_shape_;
            StdLargeVec<Vecd> force_from_fluid_;
            StdLargeVec<Real> kernel_gradient_rij_;

 
        };

        /**
         * @class StaticConfinementViscousAcceleration
         * @brief static confinement condition for viscous acceleration
         */
        class StaticConfinementViscousAccelerationForce : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementViscousAccelerationForce(NearShapeSurface& near_surface)
                : BaseLocalDynamics<BodyPartByCell>(near_surface), FluidDataSimple(sph_body_),
			pos_(particles_->pos_), mass_(particles_->mass_), force_prior_(particles_->force_prior_), rho_(particles_->rho_),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), vel_(particles_->vel_),
			level_set_shape_(&near_surface.level_set_shape_)
		{}
            virtual ~StaticConfinementViscousAccelerationForce() {};
            void interaction(size_t index_i, Real dt = 0.0)
            {
                Vecd force = Vecd::Zero();
			    Vecd vel_derivative = Vecd::Zero();
			    Vecd vel_level_set_cell_j = Vecd::Zero();
			    Real rho_i = rho_[index_i];
			    /*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			    Real phi_r_ij = abs(level_set_shape_->findSignedDistance(pos_[index_i]));
			    vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			    Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
			    force = 2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
			    force_from_fluid_[index_i] += force;
            }
            StdLargeVec<Vecd> &getForceFromFluid() { return force_prior_; };
        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<Real>& rho_, &mass_;
            StdLargeVec<Vecd>& vel_, &force_prior_;
            Real mu_;
            LevelSetShape* level_set_shape_;
            StdLargeVec<Vecd> force_from_fluid_;
 
        };


        class BaseForceFromFluidStaticConfinement : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public: 
            explicit BaseForceFromFluidStaticConfinement(NearShapeSurface& near_surface);
            virtual ~BaseForceFromFluidStaticConfinement() {};
            StdLargeVec<Vecd> &getForceFromFluid() { return force_from_fluid_; };

        protected:
            StdLargeVec<Vecd> force_from_fluid_;
            LevelSetShape* level_set_shape_;
            //StdLargeVec<Real> &Vol_;
        };

        class ViscousForceFromFluidStaticConfinement : public BaseForceFromFluidStaticConfinement
        {
        public:
            explicit ViscousForceFromFluidStaticConfinement(NearShapeSurface& near_surface);
            virtual ~ViscousForceFromFluidStaticConfinement() {};
            inline void interaction(size_t index_i, Real dt = 0.0)
            {
			    Vecd acceleration = Vecd::Zero();
			    Vecd vel_derivative = Vecd::Zero();
			    Vecd vel_level_set_cell_j = Vecd::Zero();
			    Real rho_i = rho_[index_i];
			    /*Here we give the Level-set boundary velocity as zero, but later we need a vector to set the velocity of each level-set cell*/
			    vel_derivative = 2.0 * (vel_[index_i] - vel_level_set_cell_j);
			    Real kernel_gradient_divide_Rij = level_set_shape_->computeKernelGradientDivideRijIntegral(pos_[index_i]);
			    force_from_fluid_[index_i]= -2.0 * mu_ * mass_[index_i] * kernel_gradient_divide_Rij * vel_derivative /rho_i;
                //force_from_fluid_[index_i]= -2.0 * mu_ * kernel_gradient_divide_Rij * vel_derivative ;
            }
        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<Real>& rho_, &mass_;
            StdLargeVec<Vecd>& vel_;
            Real mu_;
        };
        /**
        * @class StaticConfinementIntegration1stHalf
        * @brief static confinement condition for pressure relaxation
        */
        class StaticConfinementExtendIntegration1stHalf : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementExtendIntegration1stHalf(NearShapeSurface& near_surface, Real penalty_strength = 2.0);
            virtual ~StaticConfinementExtendIntegration1stHalf() {};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real penalty_strength_;
            Fluid& fluid_;
            StdLargeVec<Real>&rho_, &p_, &mass_;
            StdLargeVec<Vecd>&pos_, &vel_, &force_;
            LevelSetShape* level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };

        /**
        * @class StaticConfinementIntegration1stHalf
        * @brief static confinement condition for pressure relaxation
        */
        class StaticConfinementIntegration1stHalfPenaltyVelocity : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementIntegration1stHalfPenaltyVelocity(NearShapeSurface& near_surface, Real  sound_speed, Real penalty_strength = 2.0);
            virtual ~StaticConfinementIntegration1stHalfPenaltyVelocity() {};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real penalty_strength_, c_0_;
            Fluid& fluid_;
            StdLargeVec<Real>& rho_, & p_, &mass_;
            StdLargeVec<Vecd>& pos_, & vel_, &force_;
            LevelSetShape* level_set_shape_;
            AcousticRiemannSolver riemann_solver_;
        };

         /**
        * @class StaticConfinementFreeSurfaceIndication
        * @brief static confinement condition for free surface particle indicate
        */
        class StaticConfinementFreeSurfaceIndication : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementFreeSurfaceIndication(NearShapeSurface& near_surface);
            virtual ~StaticConfinementFreeSurfaceIndication() {};
            void interaction(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd>& pos_;
            StdLargeVec<Real>& pos_div_;
            StdLargeVec<int> &surface_indicator_;
            LevelSetShape* level_set_shape_;
        };

        /**
        * @class StaticConfinementIntegration1stHalf
        * @brief static confinement condition for pressure relaxation
        */
        class StaticConfinementBounding : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
        {
        public:
            StaticConfinementBounding(NearShapeSurface& near_surface);
            virtual ~StaticConfinementBounding() {};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd>& pos_;
            LevelSetShape* level_set_shape_;
            Real constrained_distance_;
        };

        /**
        * @class StaticConfinement
        * @brief Static confined boundary condition for complex structures with bounding.
        */
         class StaticConfinementWithBounding
         {
         public:

             SimpleDynamics<StaticConfinementDensity> density_summation_;
             SimpleDynamics<StaticConfinementIntegration1stHalf> pressure_relaxation_;
             SimpleDynamics<StaticConfinementIntegration2ndHalf> density_relaxation_;
             SimpleDynamics<StaticConfinementBounding> surface_bounding_;


             StaticConfinementWithBounding(NearShapeSurface& near_surface);
             virtual ~StaticConfinementWithBounding() {};
         };


        /**
        * @class StaticConfinement
        * @brief Static confined boundary condition for complex structures with penalty force for light phase.
        */
        class StaticConfinementWithPenalty
        {
        public:

            SimpleDynamics<StaticConfinementDensity> density_summation_;
            SimpleDynamics<StaticConfinementIntegration1stHalf> pressure_relaxation_;
            SimpleDynamics<StaticConfinementIntegration2ndHalf> density_relaxation_;
            SimpleDynamics<StaticConfinementTransportVelocity> transport_velocity_;
            SimpleDynamics<StaticConfinementViscousAcceleration, SequencedPolicy> viscous_acceleration_;
            SimpleDynamics<StaticConfinementExtendIntegration1stHalf> extend_intergration_1st_half_;
            SimpleDynamics<StaticConfinementIntegration1stHalfPenaltyVelocity> extend_intergration_1st_half_Velocity;
            SimpleDynamics<StaticConfinementBounding> surface_bounding_;

            StaticConfinementWithPenalty(NearShapeSurface& near_surface, Real sound_speed, Real penalty_strength);
            virtual ~StaticConfinementWithPenalty() {};
        };

        class StaticConfinementGeneral
        {
        public:
            SimpleDynamics<StaticConfinementDensity, SequencedPolicy> density_summation_;
            SimpleDynamics<StaticConfinementIntegration1stHalf> pressure_relaxation_;
            SimpleDynamics<StaticConfinementIntegration2ndHalf> density_relaxation_;
            SimpleDynamics<StaticConfinementTransportVelocity, SequencedPolicy> transport_velocity_;
            SimpleDynamics<StaticConfinementViscousAcceleration, SequencedPolicy> viscous_acceleration_;
            InteractionDynamics<StaticConfinementFreeSurfaceIndication> free_surface_indication_;
            SimpleDynamics<StaticConfinementBounding> surface_bounding_;

            StaticConfinementGeneral(NearShapeSurface &near_surface);
            virtual ~StaticConfinementGeneral(){};
        };

    }
}
#endif  LEVEL_SET_COFINEMENT_H
