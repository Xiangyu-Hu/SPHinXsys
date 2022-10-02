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
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		typedef DataDelegateContact<FluidParticles, SolidParticles, DataDelegateEmptyBase>
			FluidWallData;
		typedef DataDelegateContact<FluidParticles, BaseParticles, DataDelegateEmptyBase>
			FluidContactData;
		typedef DataDelegateContact<FluidParticles, SolidParticles> FSIContactData;
		/**
		 * @class RelaxationWithWall
		 * @brief Base class adding interaction with wall to general relaxation process
		 */
		template <class BaseRelaxationType>
		class RelaxationWithWall : public BaseRelaxationType, public FluidWallData
		{
		public:
			template <class BaseBodyRelationType, typename... Args>
			RelaxationWithWall(BaseContactRelation &wall_contact_relation,
							   BaseBodyRelationType &base_body_relation, Args &&...args);
			template <typename... Args>
			RelaxationWithWall(ComplexRelation &fluid_wall_relation, Args &&...args)
				: RelaxationWithWall(fluid_wall_relation.contact_relation_,
									 fluid_wall_relation.inner_relation_, std::forward<Args>(args)...) {}
			virtual ~RelaxationWithWall(){};

		protected:
			StdVec<Real> wall_inv_rho0_;
			StdVec<StdLargeVec<Real> *> wall_mass_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;
		};
		/**
		 * @class DensitySummation
		 * @brief computing density by summation considering contribution from contact bodies
		 */
		template <class DensitySummationInnerType>
		class DensitySummation : public BaseInteractionComplex<DensitySummationInnerType, FluidContactData>
		{
		public:
			template <typename... Args>
			DensitySummation(Args &&...args);
			virtual ~DensitySummation(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;
		};
		/** the instance without considering free surface */
		using DensitySummationComplex = DensitySummation<DensitySummationInner>;

		/**
		 * @class BaseViscousAccelerationWithWall
		 * @brief template class viscous acceleration together with wall boundary condition
		 */
		template <class ViscousAccelerationInnerType>
		class BaseViscousAccelerationWithWall : public RelaxationWithWall<ViscousAccelerationInnerType>
		{
		public:
			template <typename... Args>
			BaseViscousAccelerationWithWall(Args &&...args)
				: RelaxationWithWall<ViscousAccelerationInnerType>(std::forward<Args>(args)...){};
			virtual ~BaseViscousAccelerationWithWall(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		using ViscousAccelerationWithWall = BaseViscousAccelerationWithWall<ViscousAccelerationInner>;

		/**
		 * @class TransportVelocityCorrectionComplex
		 * @brief  transport velocity correction considering the contribution from contact bodies
		 */
		class TransportVelocityCorrectionComplex
			: public BaseInteractionComplex<TransportVelocityCorrectionInner, FluidContactData>
		{
		public:
			template <typename... Args>
			TransportVelocityCorrectionComplex(Args &&...args)
				: BaseInteractionComplex<TransportVelocityCorrectionInner, FluidContactData>(
					  std::forward<Args>(args)...){};
			virtual ~TransportVelocityCorrectionComplex(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class BasePressureRelaxationWithWall
		 * @brief  template class pressure relaxation scheme together with wall boundary
		 */
		template <class BasePressureRelaxationType>
		class BasePressureRelaxationWithWall : public RelaxationWithWall<BasePressureRelaxationType>
		{
		public:
			template <typename... Args>
			BasePressureRelaxationWithWall(Args &&...args)
				: RelaxationWithWall<BasePressureRelaxationType>(std::forward<Args>(args)...){};
			virtual ~BasePressureRelaxationWithWall(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using PressureRelaxationWithWall = BasePressureRelaxationWithWall<PressureRelaxationInner>;
		using PressureRelaxationRiemannWithWall = BasePressureRelaxationWithWall<PressureRelaxationRiemannInner>;

		/**
		 * @class BaseExtendPressureRelaxationWithWall
		 * @brief template class for pressure relaxation scheme with wall boundary
		 * and considering non-conservative acceleration term and wall penalty to prevent
		 * particle penetration.
		 */
		template <class BasePressureRelaxationType>
		class BaseExtendPressureRelaxationWithWall : public BasePressureRelaxationWithWall<BasePressureRelaxationType>
		{
		public:
			template <class BaseBodyRelationType, typename... Args>
			BaseExtendPressureRelaxationWithWall(BaseContactRelation &wall_contact_relation,
												 BaseBodyRelationType &base_body_relation,
												 Args &&...args, Real penalty_strength = 1.0)
				: BasePressureRelaxationWithWall<BasePressureRelaxationType>(
					  wall_contact_relation, base_body_relation, std::forward<Args>(args)...),
				  penalty_strength_(penalty_strength)
			{
				this->particles_->registerVariable(non_cnsrv_acc_, "NonConservativeAcceleration");
			};
			template <typename... Args>
			BaseExtendPressureRelaxationWithWall(ComplexRelation &fluid_wall_relation,
												 Args &&...args, Real penalty_strength = 1.0)
				: BaseExtendPressureRelaxationWithWall(fluid_wall_relation.contact_relation_,
													   fluid_wall_relation.inner_relation_,
													   std::forward<Args>(args)..., penalty_strength){};
			virtual ~BaseExtendPressureRelaxationWithWall(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real penalty_strength_;
			StdLargeVec<Vecd> non_cnsrv_acc_;

			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using ExtendPressureRelaxationRiemannWithWall = BaseExtendPressureRelaxationWithWall<PressureRelaxationRiemannInner>;

		/**
		 * @class BaseDensityRelaxationWithWall
		 * @brief template density relaxation scheme without using different Riemann solvers.
		 * The difference from the free surface version is that no Riemann problem is applied
		 */
		template <class BaseDensityRelaxationType>
		class BaseDensityRelaxationWithWall : public RelaxationWithWall<BaseDensityRelaxationType>
		{
		public:
			template <typename... Args>
			BaseDensityRelaxationWithWall(Args &&...args)
				: RelaxationWithWall<BaseDensityRelaxationType>(std::forward<Args>(args)...){};
			virtual ~BaseDensityRelaxationWithWall(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		using DensityRelaxationWithWall = BaseDensityRelaxationWithWall<DensityRelaxationInner>;
		using DensityRelaxationRiemannWithWall = BaseDensityRelaxationWithWall<DensityRelaxationRiemannInner>;

		/**
		 * @class PressureRelaxationWithWallOldroyd_B
		 * @brief  first half of the pressure relaxation scheme using Riemann solver.
		 */
		class PressureRelaxationWithWallOldroyd_B : public BasePressureRelaxationWithWall<PressureRelaxationInnerOldroyd_B>
		{
		public:
			explicit PressureRelaxationWithWallOldroyd_B(ComplexRelation &fluid_wall_relation)
				: BasePressureRelaxationWithWall<PressureRelaxationInnerOldroyd_B>(fluid_wall_relation){};

			virtual ~PressureRelaxationWithWallOldroyd_B(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class DensityRelaxationWithWallOldroyd_B
		 * @brief  second half of the pressure relaxation scheme using Riemann solver.
		 */
		class DensityRelaxationWithWallOldroyd_B : public BaseDensityRelaxationWithWall<DensityRelaxationInnerOldroyd_B>
		{
		public:
			explicit DensityRelaxationWithWallOldroyd_B(ComplexRelation &fluid_wall_relation)
				: BaseDensityRelaxationWithWall<DensityRelaxationInnerOldroyd_B>(fluid_wall_relation){};

			virtual ~DensityRelaxationWithWallOldroyd_B(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif // FLUID_DYNAMICS_COMPLEX_H