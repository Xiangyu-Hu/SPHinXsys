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
* @file fluid_dynamics_complex.h
* @brief Here, we define the algorithm classes for complex fluid dynamics, 
* which is involving with either solid walls (with suffix WithWall) 
* or/and other bodies treated as wall for the fluid (with suffix Complex).   
* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef FLUID_DYNAMICS_COMPLEX_H
#define FLUID_DYNAMICS_COMPLEX_H

#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"

namespace SPH
{
	namespace fluid_dynamics
	{
		typedef DataDelegateContact<FluidBody, FluidParticles, Fluid,
									SolidBody, SolidParticles, Solid, DataDelegateEmptyBase>
			FluidWallData;
		typedef DataDelegateContact<FluidBody, FluidParticles, Fluid,
									SPHBody, BaseParticles, BaseMaterial, DataDelegateEmptyBase>
			FluidContactData;
		typedef DataDelegateContact<FluidBody, FluidParticles, Fluid,
									SolidBody, SolidParticles, Solid>
			FSIContactData;
		/**
		* @class RelaxationWithWall
		* @brief Abstract base class for general relaxation algorithms with wall
		*/
		template <class BaseRelaxationType>
		class RelaxationWithWall : public BaseRelaxationType, public FluidWallData
		{
		public:
			template <class BaseBodyRelationType>
			RelaxationWithWall(BaseBodyRelationType &base_body_relation,
							   BaseBodyRelationContact &wall_contact_relation);
			virtual ~RelaxationWithWall(){};

		protected:
			StdVec<Real> wall_inv_rho0_;
			StdVec<StdLargeVec<Real> *> wall_mass_, wall_Vol_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_dvel_dt_ave_, wall_n_;
		};

		/**
		* @class FreeSurfaceIndicationComplex
		* @brief indicate the particles near the free fluid surface.
		*/
		class FreeSurfaceIndicationComplex : public FreeSurfaceIndicationInner, public FluidContactData
		{
		public:
			FreeSurfaceIndicationComplex(BaseBodyRelationInner &inner_relation,
										 BaseBodyRelationContact &contact_relation, Real thereshold = 0.75);
			explicit FreeSurfaceIndicationComplex(ComplexBodyRelation &complex_relation, Real thereshold = 0.75);
			virtual ~FreeSurfaceIndicationComplex(){};

		protected:
			StdVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using SpatialTemporalFreeSurfaceIdentificationComplex =
			SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIndicationComplex>;

		/**
		* @class DensitySummation
		* @brief computing density by summation considering  contribution from contact bodies
		*/
		template <class DensitySummationInnerType>
		class DensitySummation : public ParticleDynamicsComplex<DensitySummationInnerType, FluidContactData>
		{
		public:
			DensitySummation(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation);
			explicit DensitySummation(ComplexBodyRelation &complex_relation);
			DensitySummation(ComplexBodyRelation &complex_relation, BaseBodyRelationContact &extra_contact_relation);
			virtual ~DensitySummation(){};

		protected:
			StdVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;

			virtual void prepareContactData() override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		/** the case without free surface */
		using DensitySummationComplex = DensitySummation<DensitySummationInner>;
		/** the case with free surface */
		using DensitySummationFreeSurfaceComplex = DensitySummation<DensitySummationFreeSurfaceInner>;
		/** the case with free stream */
		using DensitySummationFreeStreamComplex = DensitySummation<DensitySummationFreeStreamInner>;

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
			virtual ~ViscousWithWall(){};

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
		 * @class TransportVelocityCorrectionComplex
		 * @brief  transport velocity correction consdiering  the contribution from contact bodies
		 */
		class TransportVelocityCorrectionComplex
			: public ParticleDynamicsComplex<TransportVelocityCorrectionInner, FluidContactData>
		{
		public:
			TransportVelocityCorrectionComplex(BaseBodyRelationInner &inner_relation,
											   BaseBodyRelationContact &contact_relation);

			explicit TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation);

			TransportVelocityCorrectionComplex(ComplexBodyRelation &complex_relation,
											   BaseBodyRelationContact &extra_contact_relation);
			virtual ~TransportVelocityCorrectionComplex(){};

		protected:
			StdVec<StdLargeVec<Real> *> contact_Vol_;

			virtual void prepareContactData() override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

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
			virtual ~PressureRelaxation(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		/**
		 * @class ExtendPressureRelaxation
		 * @brief template class for pressure relaxation scheme with wall boundary
		 * and considering non-conservative acceleration term and wall penalty to prevent 
		 * particle penetration.
		 */
		template <class BasePressureRelaxationType>
		class ExtendPressureRelaxation : public PressureRelaxation<BasePressureRelaxationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			ExtendPressureRelaxation(BaseBodyRelationType &base_body_relation,
									 BaseBodyRelationContact &wall_contact_relation, Real penalty_strength = 1.0);

			virtual ~ExtendPressureRelaxation(){};

		protected:
			Real penalty_strength_;
			StdLargeVec<Vecd> non_cnsrv_dvel_dt_;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
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
		using PressureRelaxationRiemannWithWall = BasePressureRelaxationWithWall<PressureRelaxation<PressureRelaxationRiemannInner>>;

		/** template interface class for the extended pressure relaxation with wall schemes */
		template <class BasePressureRelaxationType>
		class ExtendPressureRelaxationWithWall : public BasePressureRelaxationType
		{
		public:
			explicit ExtendPressureRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation, Real penalty_strength = 1.0);
			ExtendPressureRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation,
											 BaseBodyRelationContact &wall_contact_relation, Real penalty_strength = 1.0);
			ExtendPressureRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation,
											 BaseBodyRelationContact &wall_contact_relation, Real penalty_strength = 1.0);
		};
		using ExtendPressureRelaxationRiemannWithWall = ExtendPressureRelaxationWithWall<ExtendPressureRelaxation<PressureRelaxationRiemannInner>>;

		/**
		* @class DensityRelaxation
		* @brief template density relaxation scheme without using different Riemann solvers.
		* The difference from the free surface version is that no Riemann problem is applied
		*/
		template <class BaseDensityRelaxationType>
		class DensityRelaxation : public RelaxationWithWall<BaseDensityRelaxationType>
		{
		public:
			// template for different combination of constructing body relations
			template <class BaseBodyRelationType>
			DensityRelaxation(BaseBodyRelationType &base_body_relation,
							  BaseBodyRelationContact &wall_contact_relation);
			virtual ~DensityRelaxation(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/** template interface class for different density relaxation schemes */
		template <class BaseDensityRelaxationType>
		class BaseDensityRelaxationWithWall : public DensityRelaxation<BaseDensityRelaxationType>
		{
		public:
			explicit BaseDensityRelaxationWithWall(ComplexBodyRelation &fluid_wall_relation);
			BaseDensityRelaxationWithWall(BaseBodyRelationInner &fluid_inner_relation,
										  BaseBodyRelationContact &wall_contact_relation);
			BaseDensityRelaxationWithWall(ComplexBodyRelation &fluid_complex_relation,
										  BaseBodyRelationContact &wall_contact_relation);
		};
		using DensityRelaxationWithWall = BaseDensityRelaxationWithWall<DensityRelaxationInner>;
		using DensityRelaxationRiemannWithWall = BaseDensityRelaxationWithWall<DensityRelaxationRiemannInner>;

		/**
		* @class PressureRelaxationWithWallOldroyd_B
		* @brief  first half of the pressure relaxation scheme using Riemann solver.
		*/
		class PressureRelaxationWithWallOldroyd_B : public PressureRelaxation<PressureRelaxationInnerOldroyd_B>
		{
		public:
			explicit PressureRelaxationWithWallOldroyd_B(ComplexBodyRelation &fluid_wall_relation)
				: PressureRelaxation<PressureRelaxationInnerOldroyd_B>(fluid_wall_relation.inner_relation_,
																			  fluid_wall_relation.contact_relation_){};

			virtual ~PressureRelaxationWithWallOldroyd_B(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DensityRelaxationWithWallOldroyd_B
		* @brief  second half of the pressure relaxation scheme using Riemann solver.
		*/
		class DensityRelaxationWithWallOldroyd_B : public DensityRelaxation<DensityRelaxationInnerOldroyd_B>
		{
		public:
			explicit DensityRelaxationWithWallOldroyd_B(ComplexBodyRelation &fluid_wall_relation)
				: DensityRelaxation<DensityRelaxationInnerOldroyd_B>(fluid_wall_relation.inner_relation_,
																			fluid_wall_relation.contact_relation_){};

			virtual ~DensityRelaxationWithWallOldroyd_B(){};

		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ColorFunctionGradientComplex
		* @brief indicate the particles near the free fluid surface.
		*/
		class ColorFunctionGradientComplex : public ColorFunctionGradientInner, public FluidContactData
		{
		public:
			ColorFunctionGradientComplex(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation);
			ColorFunctionGradientComplex(ComplexBodyRelation &complex_relation);
			virtual ~ColorFunctionGradientComplex(){};

		protected:
			StdVec<StdLargeVec<Real> *> contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class 	SurfaceNormWithWall
		 * @brief  Modify surface norm when contact with wall
		 */
		class SurfaceNormWithWall : public InteractionDynamics, public FSIContactData
		{
		public:
			SurfaceNormWithWall(BaseBodyRelationContact &contact_relation, Real contact_angle);
			virtual ~SurfaceNormWithWall(){};

		protected:
			Real contact_angle_;
			Real smoothing_length_;
			Real particle_spacing_;
			StdLargeVec<int> &surface_indicator_;
			StdLargeVec<Vecd> &surface_norm_;
			StdLargeVec<Real> &pos_div_;
			StdVec<StdLargeVec<Vecd> *> wall_n_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //FLUID_DYNAMICS_COMPLEX_H