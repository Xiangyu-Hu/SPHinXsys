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
 * @file 	fluid_dynamics_complex.h
 * @brief 	Here, we define the algorithm classes for complex fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
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
		 * @class InteractionWithWall
		 * @brief Base class adding interaction with wall to general relaxation process
		 */
		template <class BaseIntegrationType>
		class InteractionWithWall : public BaseIntegrationType, public FluidWallData
		{
		public:
			template <class BaseBodyRelationType, typename... Args>
			InteractionWithWall(BaseContactRelation &wall_contact_relation,
								BaseBodyRelationType &base_body_relation, Args &&...args);
			template <typename... Args>
			InteractionWithWall(ComplexRelation &fluid_wall_relation, Args &&...args)
				: InteractionWithWall(fluid_wall_relation.getContactRelation(),
									  fluid_wall_relation.getInnerRelation(), std::forward<Args>(args)...) {}
			virtual ~InteractionWithWall(){};

		protected:
			StdVec<Real> wall_inv_rho0_;
			StdVec<StdLargeVec<Real> *> wall_mass_;
			StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;
		};

        class BaseDensitySummationComplexKernel {
        public:
            BaseDensitySummationComplexKernel(Real *contactInvRho0, Real **contactMass,
                                              NeighborhoodDevice **contactConfiguration,
                                              size_t contactConfigurationSize) :
                                              contact_inv_rho0_(contactInvRho0),
                                              contact_mass_(contactMass),
                                              contact_configuration_(contactConfiguration),
                                              contact_configuration_size_(contactConfigurationSize) {}

            template<class ContactMassFunc, class ContactConfigFunc>
            static Real ContactSummation(size_t index_i, std::size_t contact_configuration_size,
                                         const Real* contact_inv_rho0, ContactMassFunc&& contactMassFunc,
                                         ContactConfigFunc&& contactConfigFunc)
            {
                Real sigma(0.0);
                for (size_t k = 0; k < contact_configuration_size; ++k)
                {
                    Real* contact_mass_k = contactMassFunc(k);
                    Real contact_inv_rho0_k = contact_inv_rho0[k];
                    auto& contact_neighborhood = contactConfigFunc(k, index_i);
                    for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                        sigma += contact_neighborhood.W_ij_[n] * contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
                }
                return sigma;
            }

            Real ContactSummation(size_t index_i)
            {
                return BaseDensitySummationComplexKernel::ContactSummation(index_i, contact_configuration_size_,
                           contact_inv_rho0_, [&](auto k){ return contact_mass_[k]; },
                           [&](auto k, auto index_i) -> NeighborhoodDevice& { return contact_configuration_[k][index_i]; });
            };

        private:
            Real *contact_inv_rho0_, **contact_mass_;
            NeighborhoodDevice** contact_configuration_;
            std::size_t contact_configuration_size_;
        };

		/**
		 * @class DensitySummation
		 * @brief computing density by summation considering contribution from contact bodies
		 */
		template <class DensitySummationInnerType>
		class BaseDensitySummationComplex
			: public BaseInteractionComplex<DensitySummationInnerType, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit BaseDensitySummationComplex(Args &&...args);
			virtual ~BaseDensitySummationComplex(){};

		protected:
			StdSharedVec<Real> contact_inv_rho0_;
			StdVec<StdLargeVec<Real> *> contact_mass_;
            StdSharedVec<Real*> contact_mass_device_;

            ExecutionProxy<BaseDensitySummationComplex, BaseDensitySummationComplexKernel> device_proxy;

			Real ContactSummation(size_t index_i);

            auto& getDeviceProxy() {
                return device_proxy;
            }
		};

        class DensitySummationComplexKernel : public BaseDensitySummationComplexKernel,
                                              public DensitySummationInnerKernel {
        public:
            DensitySummationComplexKernel(const BaseDensitySummationComplexKernel& complexKernel,
                                          const DensitySummationInnerKernel& innerKernel)
            : BaseDensitySummationComplexKernel(complexKernel), DensitySummationInnerKernel(innerKernel) {}

            template<class InteractionFunc, class ContactSummationFunc>
            static void interaction(size_t index_i, Real dt, Real* rho_sum, Real rho0, Real inv_sigma0, const Real* mass,
                                    InteractionFunc&& interaction, ContactSummationFunc&& ContactSummation)
            {
                interaction(index_i, dt);
                Real sigma = ContactSummation(index_i);
                rho_sum[index_i] += sigma * rho0 * rho0 * inv_sigma0 / mass[index_i];
            }

            void interaction(size_t index_i, Real dt)
            {
                interaction(index_i, dt, rho_sum_, rho0_, inv_sigma0_, mass_,
                            [&](auto idx, auto delta) { DensitySummationInnerKernel::interaction(idx, delta); },
                            [&](auto idx) { return BaseDensitySummationComplexKernel::ContactSummation(idx); });
            }

        };

		/**
		 * @class DensitySummationComplex
		 * @brief computing density by summation considering contribution from contact bodies
		 */
		class DensitySummationComplex
			: public BaseDensitySummationComplex<DensitySummationInner>
		{
		public:
			template <typename... Args>
			explicit DensitySummationComplex(Args &&...args)
				: BaseDensitySummationComplex<DensitySummationInner>(std::forward<Args>(args)...),
                  device_proxy(this, *BaseDensitySummationComplex<DensitySummationInner>::getDeviceProxy().getKernel(),
                               *DensitySummationInner::getDeviceProxy().getKernel()){};
			virtual ~DensitySummationComplex(){};

			inline void interaction(size_t index_i, Real dt = 0.0);

            auto& getDeviceProxy() {
                return device_proxy;
            }

        private:
            ExecutionProxy<DensitySummationComplex, DensitySummationComplexKernel> device_proxy;
		};

		/**
		 * @class DensitySummationComplexAdaptive
		 * @brief computing density by summation considering  contribution from contact bodies
		 */
		class DensitySummationComplexAdaptive
			: public BaseDensitySummationComplex<DensitySummationInnerAdaptive>
		{
		public:
			template <typename... Args>
			explicit DensitySummationComplexAdaptive(Args &&...args)
			: BaseDensitySummationComplex<DensitySummationInnerAdaptive>(std::forward<Args>(args)...){};
			virtual ~DensitySummationComplexAdaptive(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class ViscousWithWall
		 * @brief  template class viscous acceleration with wall boundary
		 */
		template <class ViscousAccelerationInnerType>
		class BaseViscousAccelerationWithWall : public InteractionWithWall<ViscousAccelerationInnerType>
		{
		public:
			template <typename... Args>
			BaseViscousAccelerationWithWall(Args &&...args)
				: InteractionWithWall<ViscousAccelerationInnerType>(std::forward<Args>(args)...){};
			virtual ~BaseViscousAccelerationWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
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
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class TransportVelocityCorrectionComplexAdaptive
		 * @brief  transport velocity correction considering the contribution from contact bodies
		 */
		class TransportVelocityCorrectionComplexAdaptive
			: public BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactData>
		{
		public:
			template <typename... Args>
			TransportVelocityCorrectionComplexAdaptive(Args &&...args)
				: BaseInteractionComplex<TransportVelocityCorrectionInnerAdaptive, FluidContactData>(
					  std::forward<Args>(args)...){};
			virtual ~TransportVelocityCorrectionComplexAdaptive(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class BaseIntegration1stHalfWithWall
		 * @brief  template class pressure relaxation scheme together with wall boundary
		 */
		template <class BaseIntegration1stHalfType>
		class BaseIntegration1stHalfWithWall : public InteractionWithWall<BaseIntegration1stHalfType>
		{
		public:
			template <typename... Args>
			BaseIntegration1stHalfWithWall(Args &&...args)
				: InteractionWithWall<BaseIntegration1stHalfType>(std::forward<Args>(args)...){};
			virtual ~BaseIntegration1stHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);

		protected:
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using Integration1stHalfWithWall = BaseIntegration1stHalfWithWall<Integration1stHalf>;
		using Integration1stHalfRiemannWithWall = BaseIntegration1stHalfWithWall<Integration1stHalfRiemann>;

		/**
		 * @class BaseExtendIntegration1stHalfWithWall
		 * @brief template class for pressure relaxation scheme with wall boundary
		 * and considering non-conservative acceleration term and wall penalty to prevent
		 * particle penetration.
		 */
		template <class BaseIntegration1stHalfType>
		class BaseExtendIntegration1stHalfWithWall : public BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>
		{
		public:
			template <class BaseBodyRelationType, typename... Args>
			BaseExtendIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation,
												 BaseBodyRelationType &base_body_relation,
												 Args &&...args, Real penalty_strength = 1.0)
				: BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>(
					  wall_contact_relation, base_body_relation, std::forward<Args>(args)...),
				  penalty_strength_(penalty_strength)
			{
				this->particles_->registerVariable(non_cnsrv_acc_, "NonConservativeAcceleration");
			};
			template <typename... Args>
			BaseExtendIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation,
												 Args &&...args, Real penalty_strength = 1.0)
				: BaseExtendIntegration1stHalfWithWall(fluid_wall_relation.getContactRelation(),
													   fluid_wall_relation.getInnerRelation(),
													   std::forward<Args>(args)..., penalty_strength){};
			virtual ~BaseExtendIntegration1stHalfWithWall(){};
			void initialization(size_t index_i, Real dt = 0.0);
			
			inline void interaction(size_t index_i, Real dt = 0.0);

		protected:
			Real penalty_strength_;
			StdLargeVec<Vecd> non_cnsrv_acc_;

			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using ExtendIntegration1stHalfRiemannWithWall = BaseExtendIntegration1stHalfWithWall<Integration1stHalfRiemann>;

		/**
		 * @class BaseIntegration2ndHalfWithWall
		 * @brief template density relaxation scheme without using different Riemann solvers.
		 * The difference from the free surface version is that no Riemann problem is applied
		 */
		template <class BaseIntegration2ndHalfType>
		class BaseIntegration2ndHalfWithWall : public InteractionWithWall<BaseIntegration2ndHalfType>
		{
		public:
			template <typename... Args>
			BaseIntegration2ndHalfWithWall(Args &&...args)
				: InteractionWithWall<BaseIntegration2ndHalfType>(std::forward<Args>(args)...){};
			virtual ~BaseIntegration2ndHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		using Integration2ndHalfWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalf>;
		using Integration2ndHalfRiemannWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalfRiemann>;

		/**
		 * @class Oldroyd_BIntegration1stHalfWithWall
		 * @brief  first half of the pressure relaxation scheme using Riemann solver.
		 */
		class Oldroyd_BIntegration1stHalfWithWall : public BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>
		{
		public:
			explicit Oldroyd_BIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
				: BaseIntegration1stHalfWithWall<Oldroyd_BIntegration1stHalf>(fluid_wall_relation){};

			virtual ~Oldroyd_BIntegration1stHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class Oldroyd_BIntegration2ndHalfWithWall
		 * @brief  second half of the pressure relaxation scheme using Riemann solver.
		 */
		class Oldroyd_BIntegration2ndHalfWithWall : public BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>
		{
		public:
			explicit Oldroyd_BIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
				: BaseIntegration2ndHalfWithWall<Oldroyd_BIntegration2ndHalf>(fluid_wall_relation){};

			virtual ~Oldroyd_BIntegration2ndHalfWithWall(){};
			
			inline void interaction(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif // FLUID_DYNAMICS_COMPLEX_H