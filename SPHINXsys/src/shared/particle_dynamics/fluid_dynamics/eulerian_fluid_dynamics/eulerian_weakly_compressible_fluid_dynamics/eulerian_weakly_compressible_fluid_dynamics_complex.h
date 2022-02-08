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
* @file eulerian_weakly_compressible_fluid_dynamics_complex.h
* @brief Here, we define the algorithm classes for complex weakly compressible fluid dynamics, 
* which is involving with either solid walls (with suffix WithWall) 
* or/and other bodies treated as wall for the fluid (with suffix Complex).   
* @author	Zhentong Wang,Chi Zhang and Xiangyu Hu
*/

#pragma once

#include "eulerian_weakly_compressible_fluid_dynamics_inner.h"
namespace SPH
{
	namespace eulerian_weakly_compressible_fluid_dynamics
	{
		typedef DataDelegateContact<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid,
			SolidBody, SolidParticles, Solid, DataDelegateEmptyBase>
			WCFluidWallData;
		typedef DataDelegateContact<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid,
			SPHBody, BaseParticles, BaseMaterial, DataDelegateEmptyBase> 
			WCFluidContactData;
		typedef DataDelegateContact<EulerianFluidBody, WeaklyCompressibleFluidParticles, Fluid,
			SolidBody, SolidParticles, Solid> 
			WCFSIContactData;
		/**
		* @class RelaxationWithWall
		* @brief Abstract base class for general relaxation algorithms with wall
		*/
		template <class BaseRelaxationType>
		class RelaxationWithWall : public BaseRelaxationType, public WCFluidWallData
		{
		public:
			template <class BaseBodyRelationType>
			RelaxationWithWall(BaseBodyRelationType &base_body_relation,
				BaseBodyRelationContact &wall_contact_relation);
			virtual ~RelaxationWithWall() {};

		protected:
			StdVec<Real> wall_inv_rho0_;
			StdVec<StdLargeVec<Real> *> wall_mass_, wall_Vol_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_dvel_dt_ave_, wall_n_;
		};

		/**
		 * @class ViscousWithWall
		 * @brief  template class viscous acceleration with wall boundary
		 */
		template <class BaseViscousAccelerationType>
		class ViscousWithWall : public RelaxationWithWall<BaseViscousAccelerationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			ViscousWithWall(BaseBodyRelationType &base_body_relation,
				BaseBodyRelationContact &wall_contact_relation);
			virtual ~ViscousWithWall() {};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/** template interface class for different pressure relaxation with wall schemes */
		template <class BaseViscousAccelerationType>
		class BaseViscousAccelerationWithWall : public BaseViscousAccelerationType
		{
		public:
			explicit BaseViscousAccelerationWithWall(ComplexBodyRelation &fluid_wall_relation);
			BaseViscousAccelerationWithWall(BaseBodyRelationInner &fluid_inner_relation,
				BaseBodyRelationContact &wall_contact_relation);
			BaseViscousAccelerationWithWall(ComplexBodyRelation &fluid_complex_relation,
				BaseBodyRelationContact &wall_contact_relation);
		};
		using ViscousAccelerationWithWall = BaseViscousAccelerationWithWall<ViscousWithWall<ViscousAccelerationInner>>;

		/**
		 * @class PressureRelaxation
		 * @brief  template class pressure relaxation scheme with wall boundary
		 */
		template <class BasePressureRelaxationType>
		class PressureRelaxation : public RelaxationWithWall<BasePressureRelaxationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			PressureRelaxation(BaseBodyRelationType &base_body_relation,
				BaseBodyRelationContact &wall_contact_relation);
			virtual ~PressureRelaxation() {};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/** template interface class for different pressure relaxation with wall schemes */
		template <class BasePressureRelaxationType>
		class BasePressureRelaxationWithWall : public BasePressureRelaxationType
		{
		public:
			explicit BasePressureRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation);
			BasePressureRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation,
				BaseBodyRelationContact &wall_contact_relation);
			BasePressureRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation,
				BaseBodyRelationContact &wall_contact_relation);
		};
		using PressureRelaxationWithWall = BasePressureRelaxationWithWall<PressureRelaxation<PressureRelaxationInner>>;
		using PressureRelaxationAcousticRiemannWithWall
			= BasePressureRelaxationWithWall<PressureRelaxation<PressureRelaxationAcousticRiemannInner>>;
		using PressureRelaxationHLLCRiemannWithWall
			= BasePressureRelaxationWithWall<PressureRelaxation<PressureRelaxationHLLCRiemannInner>>;
		using PressureRelaxationHLLCRiemannWithLimiterWithWall
			= BasePressureRelaxationWithWall<PressureRelaxation<PressureRelaxationHLLCWithLimiterRiemannInner>>;

		/**
		* @class DensityRelaxation
		* @brief template density relaxation scheme without using different Riemann solvers.
		* The difference from the free surface version is that no Riemann problem is applied
		*/
		template <class BaseDensityAndEnergyRelaxationType>
		class DensityAndEnergyRelaxation : public RelaxationWithWall<BaseDensityAndEnergyRelaxationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			DensityAndEnergyRelaxation(BaseBodyRelationType &base_body_relation,
				BaseBodyRelationContact &wall_contact_relation);
			virtual ~DensityAndEnergyRelaxation() {};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/** template interface class for different density relaxation schemes */
		template <class BaseDensityAndEnergyRelaxationType>
		class BaseDensityAndEnergyRelaxationWithWall : public DensityAndEnergyRelaxation<BaseDensityAndEnergyRelaxationType>
		{
		public:
			explicit BaseDensityAndEnergyRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation);
			BaseDensityAndEnergyRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation,
				BaseBodyRelationContact &wall_contact_relation);
			BaseDensityAndEnergyRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation,
				BaseBodyRelationContact &wall_contact_relation);
		};
		using DensityAndEnergyRelaxationAcousticRiemannWithWall = BaseDensityAndEnergyRelaxationWithWall<DensityAndEnergyRelaxationAcousticRiemannInner>;
		using DensityAndEnergyRelaxationHLLCRiemannWithWall = BaseDensityAndEnergyRelaxationWithWall<DensityAndEnergyRelaxationHLLCRiemannInner>;
		using DensityAndEnergyRelaxationHLLCRiemannWithLimiterWithWall = BaseDensityAndEnergyRelaxationWithWall<DensityAndEnergyRelaxationHLLCWithLimiterRiemannInner>;

		/**
		* @class FreeSurfaceIndicationComplex
		* @brief indicate the particles near the free fluid surface.
		*/
		class FreeSurfaceIndicationComplex : public FreeSurfaceIndicationInner, public WCFluidContactData
		{
		public:
			FreeSurfaceIndicationComplex(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation, Real thereshold = 0.75);
			explicit FreeSurfaceIndicationComplex(ComplexBodyRelation &complex_relation, Real thereshold = 0.75);
			virtual ~FreeSurfaceIndicationComplex() {};

		protected:
			StdVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
