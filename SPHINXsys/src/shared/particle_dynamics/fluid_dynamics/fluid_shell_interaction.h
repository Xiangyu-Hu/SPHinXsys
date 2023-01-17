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
 *  HU1527/12-1 and HU1527/12-4												*
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
* @file     fluid_shell_interaction.h
* @brief    Here, we define the algorithm classes for the interaction between fluid and
*			thin structure, plate and shell.
* @author	Chi ZHang and Xiangyu Hu
*/


#ifndef FLUID_SHELL_INTERACTION_H
#define FLUID_SHELL_INTERACTION_H

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	namespace fluid_dynamics
	{		
		typedef DataDelegateContact<FluidParticles, ShellParticles> FluidShellContactData;
		typedef DataDelegateContact<FluidParticles, ShellParticles, DataDelegateEmptyBase> FluidShellData;
		
		/**
		* @class InteractionWithShell
		* @brief Abstract base class for general fluid-shell interaction model. 
		* @note  Here, 
		*/
		template<class BaseIntegrationType>
		class InteractionWithShell : public BaseIntegrationType, public FluidShellData
		{
		public:
			template <class BaseBodyRelationType, typename... Args>
			InteractionWithShell(BaseContactRelation &contact_relation, BaseBodyRelationType &base_body_relation, Args &&...args);

			template <typename... Args>
			InteractionWithShell(ComplexRelation &fluid_shell_relation, Args &&...args)
				: InteractionWithShell(fluid_shell_relation.getContactRelation(), fluid_shell_relation.getInnerRelation(), std::forward<Args>(args)...) 
			{}

			virtual ~InteractionWithShell(){};

		protected:
			Real spacing_ref_;
			StdVec<StdLargeVec<Real> *> shell_thickness_;
			StdVec<StdLargeVec<Vecd> *> shell_n_, shell_vel_ave_, shell_acc_ave_;
		};

		/**
		 * @class BaseFluidShellIntegration1stHalf
		 * @brief  template class for fluid-shell pressure relaxation scheme
		 */
		template <class Integration1stHalfType>
		class BaseFluidShellIntegration1stHalf : public InteractionWithShell<Integration1stHalfType>
		{
		public:
			template <typename... Args>
			explicit BaseFluidShellIntegration1stHalf(Args &&...args) 
				: InteractionWithShell<Integration1stHalfType>(std::forward<Args>(args)...){};
			virtual ~BaseFluidShellIntegration1stHalf(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using FluidShellIntegration1stHalf = BaseFluidShellIntegration1stHalf<Integration1stHalf>;
		using FluidShellIntegration1stHalfRiemann = BaseFluidShellIntegration1stHalf<Integration1stHalfRiemann>;

		using FluidShellandWallIntegration1stHalf = BaseIntegration1stHalfWithWall<FluidShellIntegration1stHalf>;
		using FluidShellandWallIntegration1stHalfRiemann = BaseIntegration1stHalfWithWall<FluidShellIntegration1stHalfRiemann>;
		using ExtendFluidShellandWallIntegration1stHalfRiemann = BaseExtendIntegration1stHalfWithWall<FluidShellIntegration1stHalfRiemann>;

		/**
		 * @class BaseFluidShellIntegration2ndHalf
		 * @brief  template class pressure relaxation scheme with wall boundary
		 */
		template <class Integration2ndHalfType>
		class BaseFluidShellIntegration2ndHalf : public InteractionWithShell<Integration2ndHalfType>
		{
		public:
			template <typename... Args>
			explicit BaseFluidShellIntegration2ndHalf(Args &&...args) 
				: InteractionWithShell<Integration2ndHalfType>(std::forward<Args>(args)...)
			{};
			virtual ~BaseFluidShellIntegration2ndHalf(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};
		using FluidShellIntegration2ndHalf = BaseFluidShellIntegration2ndHalf<Integration2ndHalf>;
		using FluidShellIntegration2ndHalfRiemann = BaseFluidShellIntegration2ndHalf<Integration2ndHalfRiemann>;

		using FluidShellandWallIntegration2ndHalf = BaseIntegration2ndHalfWithWall<FluidShellIntegration2ndHalf>;
		using FluidShellandWallIntegration2ndHalfRiemann = BaseIntegration2ndHalfWithWall<FluidShellIntegration2ndHalfRiemann>;

		/**
		 * @class ViscousWithWall
		 * @brief  template class viscous acceleration with wall boundary
		 */
		template <class ViscousAccelerationInnerType>
		class BaseViscousAccelerationWithShell : public InteractionWithShell<ViscousAccelerationInnerType>
		{
		public:
			template <typename... Args>
			BaseViscousAccelerationWithShell(Args &&...args)
				: InteractionWithShell<ViscousAccelerationInnerType>(std::forward<Args>(args)...)
			{};
			virtual ~BaseViscousAccelerationWithShell(){};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		using ViscousAccelerationWithShell = BaseViscousAccelerationWithShell<ViscousAccelerationInner>;
		using ViscousAccelerationWithShellandWall = BaseViscousAccelerationWithWall<ViscousAccelerationWithShell>;

    }
}
#endif //FLUID_SHELL_INTERACTION_H