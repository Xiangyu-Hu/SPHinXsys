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
 * @file 	contact_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid contact dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
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
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, SolidBody, SolidParticles, Solid> ContactDynamicsData;
		typedef DataDelegateContact<SolidBody, SolidParticles, ElasticSolid, SolidBody, SolidParticles, Solid> ContactWithWallData;

		/**
		 * @class SelfContactDensitySummation
		 * @brief Computing the summation density due to solid self-contact model.
		 */
		class SelfContactDensitySummation : public PartInteractionDynamicsByParticle, public SolidDataInner
		{
		public:
			explicit SelfContactDensitySummation(SolidBodyRelationSelfContact &self_contact_relation);
			virtual ~SelfContactDensitySummation(){};

		protected:
			StdLargeVec<Real> self_contact_density_;
			StdLargeVec<Real> &mass_;
			Real offset_W_ij_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ContactDensitySummation
		 * @brief Computing the summation density due to solid-solid contact model.
		 */
		class ContactDensitySummation : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactDensitySummation(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactDensitySummation(){};

		protected:
			StdLargeVec<Real> contact_density_;
			StdLargeVec<Real> &mass_;
			StdVec<StdLargeVec<Real> *> contact_mass_;
			StdVec<Real> offset_W_ij_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ShellContactDensity
		 * @brief Computing the contact density due to shell contact using a
		 * 		 surface integral being solved by Gauss-Legendre quadrature integration.
		 */
		class ShellContactDensity : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ShellContactDensity(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ShellContactDensity(){};

		protected:
			StdLargeVec<Real> contact_density_;
			StdVec<StdLargeVec<Vecd> *> contact_pos_;
			StdLargeVec<Vecd> &pos_n_;

			Kernel *kernel_;
			Real spacing_ref_, boundary_factor_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;

			/** Abscissas and weights for Gauss-Legendre quadrature integration with n=3 nodes */
			Real x_0 = 0.774596669241483377035853079956;
			Real x_1 = 0.000000000000000000000000000000;
			Real x_2 = -x_0;

			Real w_0 = 0.555555555555555555555555555556;
			Real w_1 = 0.888888888888888888888888888889;
			Real w_2 = w_0;
		};

		/**
		 * @class SelfContactForce
		 * @brief Computing the self-contact force.
		 */
		class SelfContactForce : public PartInteractionDynamicsByParticle, public SolidDataInner
		{
		public:
			explicit SelfContactForce(SolidBodyRelationSelfContact &self_contact_relation);
			virtual ~SelfContactForce(){};

		protected:
			StdLargeVec<Real> &mass_, &self_contact_density_, &Vol_;
			StdLargeVec<Vecd> &dvel_dt_prior_, &vel_n_;
			Real contact_impedance_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ContactForce
		 * @brief Computing the contact force.
		 */
		class ContactForce : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactForce(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactForce(){};

		protected:
			StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
			StdLargeVec<Vecd> &dvel_dt_prior_;
			StdVec<StdLargeVec<Real> *> contact_contact_density_, contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ContactForceFromWall
		 * @brief Computing the contact force from a rigid wall.
		 *  Note that the body surface of the wall should be
		 *  updated before computing the contact force.
		 */
		class ContactForceFromWall : public PartInteractionDynamicsByParticle, public ContactWithWallData
		{
		public:
			explicit ContactForceFromWall(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactForceFromWall(){};

		protected:
			StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
			StdLargeVec<Vecd> &dvel_dt_prior_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ContactForceToWall
		 * @brief Computing contact force acting on a rigid wall.
		 */
		class ContactForceToWall : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactForceToWall(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactForceToWall(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &dvel_dt_prior_;
			StdVec<StdLargeVec<Real> *> contact_contact_density_, contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class PairwiseFrictionFromWall
		 * @brief Damping to wall by which the wall velocity is not updated
		 * and the mass of wall particle is not considered.
		 */
		class PairwiseFrictionFromWall : public InteractionDynamicsSplitting, public ContactWithWallData
		{
		public:
			PairwiseFrictionFromWall(BaseBodyRelationContact &contact_relation, Real eta);
			virtual ~PairwiseFrictionFromWall(){};

		protected:
			Real eta_; /**< friction coefficient */
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_n_;
			StdVec<StdLargeVec<Real> *> wall_Vol_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_n_, wall_n_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class DynamicContactForceWithWall
		 * @brief Computing the contact force with a rigid wall.
		 *  Note that the body surface of the wall should be
		 *  updated before computing the contact force.
		 */
		class DynamicContactForceWithWall : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit DynamicContactForceWithWall(SolidBodyRelationContact &solid_body_contact_relation, Real penalty_strength = 1.0);
			virtual ~DynamicContactForceWithWall(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_, contact_n_;
			Real penalty_strength_;
			Real impedance_, reference_pressure_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif // CONTACT_DYNAMICS_H
