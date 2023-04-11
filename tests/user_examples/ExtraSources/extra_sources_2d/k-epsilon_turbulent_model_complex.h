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
 * @file 	k-epsilon_turbulent_model.h
 * @brief 	Here,.
 * @details     T.
 * @author Xiangyu Hu
 */

#ifndef K_EPSILON_TURBULENT_MODEL_COMPLEX_H
#define K_EPSILON_TURBULENT_MODEL_COMPLEX_H

#include "k-epsilon_turbulent_model_inner.h"
#include "k-epsilon_turbulent_model_inner.hpp"

namespace SPH
{
    namespace fluid_dynamics
    {
		/**
		 * @class K_TurtbulentModelRelaxationWithWall
		 * @brief .
		 * The
		 */
		class K_TurtbulentModelComplex 
			: public BaseInteractionComplex<K_TurtbulentModelInner, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit K_TurtbulentModelComplex(Args &&...args)
				: BaseInteractionComplex<K_TurtbulentModelInner, FluidContactData>
				(std::forward<Args>(args)...) 
			{
				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					contact_vel_ave_.push_back(&(this->contact_particles_[k]->vel_));
				}
			};
			virtual ~K_TurtbulentModelComplex() {};
			inline void interaction(size_t index_i, Real dt = 0.0);
		protected:
			StdVec<StdLargeVec<Vecd>*> contact_vel_ave_;
		};

		/**
		 * @class E_TurtbulentModelRelaxationWithWall
		 * @brief .
		 */
		class E_TurtbulentModelComplex
			: public BaseInteractionComplex<E_TurtbulentModelInner, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit E_TurtbulentModelComplex(Args &&...args)
				: BaseInteractionComplex<E_TurtbulentModelInner, FluidContactData>(
					std::forward<Args>(args)...) {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class TKEnergyAccComplex
		 * @brief .
		 */
		class TKEnergyAccComplex
			: public BaseInteractionComplex<TKEnergyAccInner, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit TKEnergyAccComplex(Args &&...args)
				: BaseInteractionComplex<TKEnergyAccInner, FluidContactData>(
					std::forward<Args>(args)...) {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class ViscousWithWall
		 * @brief  template class viscous acceleration with wall boundary
		 */
		template <class TurbuViscousAccInnerType>
		class BaseTurbuViscousAccWithWall : public InteractionWithWall<TurbuViscousAccInnerType>
		{
		public:
			template <typename... Args>
			BaseTurbuViscousAccWithWall(Args &&...args)
				: InteractionWithWall<TurbuViscousAccInnerType>(std::forward<Args>(args)...) {};
			virtual ~BaseTurbuViscousAccWithWall() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		using TurbulentViscousAccelerationWithWall = BaseTurbuViscousAccWithWall<TurbuViscousAccInner>;


    }
}
#endif // K_EPSILON_TURBULENT_MODEL_COMPLEX_H