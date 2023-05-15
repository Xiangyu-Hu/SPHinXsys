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
 * @file 	eulerian_weakly_compressible_fluid_dynamics_complex.h
 * @brief 	Here, we define the algorithm classes for complex weakly compressible fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Zhentong Wang, Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"
#include "eulerian_weakly_compressible_fluid_dynamics_inner.hpp"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		typedef DataDelegateContact<WeaklyCompressibleFluidParticles, SolidParticles, DataDelegateEmptyBase>
			WCFluidWallData;
		typedef DataDelegateContact<WeaklyCompressibleFluidParticles, BaseParticles, DataDelegateEmptyBase>
			WCFluidContactData;
		typedef DataDelegateContact<WeaklyCompressibleFluidParticles, SolidParticles> WCFSIContactData;
		/**
		 * @class InteractionWithWall
		 * @brief Abstract base class for general relaxation algorithms with wall
		 */
		template <class BaseIntegrationType>
		class InteractionWithWall : public BaseIntegrationType, public WCFluidWallData
		{
		public:
			template <class BaseBodyRelationType>
			InteractionWithWall(BaseBodyRelationType &base_body_relation,
							   BaseContactRelation &wall_contact_relation);
			virtual ~InteractionWithWall(){};

		protected:
			StdVec<Real> wall_inv_rho0_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;
		};

		/**
		 * @class ViscousWithWall
		 * @brief  template class viscous acceleration with wall boundary
		 */
		template <class BaseViscousAccelerationType>
		class ViscousWithWall : public InteractionWithWall<BaseViscousAccelerationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			ViscousWithWall(BaseBodyRelationType &base_body_relation,
							BaseContactRelation &wall_contact_relation);
			virtual ~ViscousWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/** template interface class for different pressure relaxation with wall schemes */
		template <class BaseViscousAccelerationType>
		class BaseViscousAccelerationWithWall : public BaseViscousAccelerationType
		{
		public:
			explicit BaseViscousAccelerationWithWall(ComplexRelation &fluid_wall_relation);
			BaseViscousAccelerationWithWall(BaseInnerRelation &fluid_inner_relation,
											BaseContactRelation &wall_contact_relation);
			BaseViscousAccelerationWithWall(ComplexRelation &fluid_complex_relation,
											BaseContactRelation &wall_contact_relation);
		};
		using ViscousAccelerationWithWall = BaseViscousAccelerationWithWall<ViscousWithWall<ViscousAccelerationInner>>;

		/**
		 * @class Integration1stHalf
		 * @brief  template class pressure relaxation scheme with wall boundary
		 */
		template <class BaseIntegration1stHalfType>
		class BaseIntegration1stHalfWithWall : public InteractionWithWall<BaseIntegration1stHalfType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			BaseIntegration1stHalfWithWall(BaseBodyRelationType &base_body_relation,
										   BaseContactRelation &wall_contact_relation);
			explicit BaseIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
				: BaseIntegration1stHalfWithWall(fluid_wall_relation.getInnerRelation(),
												 fluid_wall_relation.getContactRelation()){};
			virtual ~BaseIntegration1stHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		using Integration1stHalfHLLCRiemannWithLimiterWithWall = BaseIntegration1stHalfWithWall<Integration1stHalfHLLCWithLimiterRiemann>;

		/**
		 * @class Integration2ndHalf
		 * @brief template density relaxation scheme without using different Riemann solvers.
		 * The difference from the free surface version is that no Riemann problem is applied
		 */
		template <class BaseIntegration2ndHalfType>
		class BaseIntegration2ndHalfWithWall : public InteractionWithWall<BaseIntegration2ndHalfType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			BaseIntegration2ndHalfWithWall(BaseBodyRelationType &base_body_relation,
							   BaseContactRelation &wall_contact_relation);
			explicit BaseIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
				: BaseIntegration2ndHalfWithWall(fluid_wall_relation.getInnerRelation(),
												 fluid_wall_relation.getContactRelation()){};
			virtual ~BaseIntegration2ndHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};
		using Integration2ndHalfHLLCRiemannWithLimiterWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalfHLLCWithLimiterRiemann>;
	}
}
