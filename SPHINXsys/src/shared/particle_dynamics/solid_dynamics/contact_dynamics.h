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
 * @file 	contact_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid contact dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef CONTACT_DYNAMICS_H
#define CONTACT_DYNAMICS_H

#include "general_solid_dynamics.h"

namespace SPH
{
	class SPHBody;
	class Kernel;
	namespace solid_dynamics
	{
		typedef DataDelegateContact<SolidParticles, SolidParticles> ContactDynamicsData;
		typedef DataDelegateContact<SolidParticles, SolidParticles> ContactWithWallData;

		/**
		 * @class SelfContactDensitySummation
		 * @brief Computing the summation density due to solid self-contact model.
		 */
		class SelfContactDensitySummation : public LocalDynamics, public SolidDataInner
		{
		public:
			explicit SelfContactDensitySummation(SelfSurfaceContactRelation &self_contact_relation);
			virtual ~SelfContactDensitySummation(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real>* self_contact_density_;
			StdLargeVec<Real> &mass_;
			Real offset_W_ij_;
		};

		/**
		 * @class ContactDensitySummation
		 * @brief Computing the summation density due to solid-solid contact model.
		 */
		class ContactDensitySummation : public LocalDynamics, public ContactDynamicsData
		{
		public:
			explicit ContactDensitySummation(SurfaceContactRelation &solid_body_contact_relation);
			virtual ~ContactDensitySummation(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real>* contact_density_;
			StdLargeVec<Real> &mass_;
			StdVec<StdLargeVec<Real> *> contact_mass_;
			StdVec<Real> offset_W_ij_;
		};

		/**
		 * @class ShellContactDensity
		 * @brief Computing the contact density due to shell contact using a
		 * 		 surface integral being solved by Gauss-Legendre quadrature integration.
		 */
		class ShellContactDensity : public LocalDynamics, public ContactDynamicsData
		{
		public:
			explicit ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation);
			virtual ~ShellContactDensity(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Solid &solid_;
			Kernel *kernel_;

			StdLargeVec<Vecd> &pos_;
			Real particle_spacing_;
			StdVec<Real> calibration_factor_;
			StdVec<Real> contact_h_ratio_;
			StdVec<Real> offset_W_ij_;
			StdVec<Real> contact_particle_spacing_;
			StdLargeVec<Real>* contact_density_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;
			StdVec<StdLargeVec<Vecd> *> contact_n_;
			StdVec<StdLargeVec<Vecd> *> contact_pos_;

			/** Abscissas and weights for Gauss-Legendre quadrature integration with n=3 nodes */
			const StdVec<Real> three_gaussian_points_ = { -0.7745966692414834, 0.0, 0.7745966692414834 };
			const StdVec<Real> three_gaussian_weights_ = { 0.5555555555555556, 0.8888888888888889, 0.5555555555555556 };

		};

		/**
		 * @class SelfContactForce
		 * @brief Computing the self-contact force.
		 */
		class SelfContactForce : public LocalDynamics, public SolidDataInner
		{
		public:
			explicit SelfContactForce(SelfSurfaceContactRelation &self_contact_relation);
			virtual ~SelfContactForce(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Solid &solid_;
			StdLargeVec<Real> &mass_, &self_contact_density_, &Vol_;
			StdLargeVec<Vecd> &acc_prior_, &vel_;
			Real contact_impedance_;
		};

		/**
		 * @class ContactForce
		 * @brief Computing the contact force.
		 */
		class ContactForce : public LocalDynamics, public ContactDynamicsData
		{
		public:
			explicit ContactForce(SurfaceContactRelation &solid_body_contact_relation);
			virtual ~ContactForce(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Solid &solid_;
			StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
			StdLargeVec<Vecd> &acc_prior_;
			StdVec<Solid *> contact_solids_;
			StdVec<StdLargeVec<Real> *> contact_contact_density_;
		};

		/**
		 * @class ContactForceFromWall
		 * @brief Computing the contact force from a rigid wall.
		 *  Note that the body surface of the wall should be
		 *  updated before computing the contact force.
		 */
		class ContactForceFromWall : public LocalDynamics, public ContactWithWallData
		{
		public:
			explicit ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation);
			virtual ~ContactForceFromWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Solid &solid_;
			StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
			StdLargeVec<Vecd> &acc_prior_;
		};

		/**
		 * @class ContactForceToWall
		 * @brief Computing contact force acting on a rigid wall.
		 */
		class ContactForceToWall : public LocalDynamics, public ContactDynamicsData
		{
		public:
			explicit ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation);
			virtual ~ContactForceToWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &acc_prior_;
			StdVec<Solid *> contact_solids_;
			StdVec<StdLargeVec<Real> *> contact_contact_density_;
		};

		/**
		 * @class PairwiseFrictionFromWall
		 * @brief Damping to wall by which the wall velocity is not updated
		 * and the mass of wall particle is not considered.
		 * Note that, currently, this class works only when the contact 
		 * bodies have the same resolution.
		 */
		class PairwiseFrictionFromWall : public LocalDynamics, public ContactWithWallData
		{
		public:
			PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta);
			virtual ~PairwiseFrictionFromWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real eta_; /**< friction coefficient */
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_n_, wall_n_;
		};

		/**
		 * @class DynamicContactForceWithWall
		 * @brief Computing the contact force with a rigid wall.
		 *  Note that the body surface of the wall should be
		 *  updated before computing the contact force.
		 */
		class DynamicContactForceWithWall : public LocalDynamics, public ContactDynamicsData
		{
		public:
			explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0);
			virtual ~DynamicContactForceWithWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Solid &solid_;
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_, &acc_prior_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_, contact_n_;
			Real penalty_strength_;
			Real impedance_, reference_pressure_;
		};
	}
}
#endif // CONTACT_DYNAMICS_H
