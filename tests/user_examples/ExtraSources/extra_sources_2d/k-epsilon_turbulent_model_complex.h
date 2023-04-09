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
		 * @class Base_K_TurtbulentModelComplex
		 * @brief Base_K_TurtbulentModelComplex
		 */
		template <class K_TurtbulentModelInnerType>
		class Base_K_TurtbulentModelComplex
			: public BaseInteractionComplex<K_TurtbulentModelInnerType, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit Base_K_TurtbulentModelComplex(Args &&...args);
			virtual ~Base_K_TurtbulentModelComplex() {};

		protected:
			StdVec<StdLargeVec<Vecd>*> contact_vel_ave_;

		};
		/**
		 * @class K_TurtbulentModelRelaxationWithWall
		 * @brief .
		 * The
		 */
		class K_TurtbulentModelWithWall 
			: public Base_K_TurtbulentModelComplex<K_TurtbulentModelInner>
		{
		public:
			template <typename... Args>
			explicit K_TurtbulentModelWithWall(Args &&...args)
				: Base_K_TurtbulentModelComplex<K_TurtbulentModelInner>(std::forward<Args>(args)...) {};
			virtual ~K_TurtbulentModelWithWall() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class Base_E_TurtbulentModelComplex
		 * @brief Base_E_TurtbulentModelComplex
		 */
		template <class E_TurtbulentModelInnerType>
		class Base_E_TurtbulentModelComplex
			: public BaseInteractionComplex<E_TurtbulentModelInnerType, FluidContactData>
		{
		public:
			template <typename... Args>
			explicit Base_E_TurtbulentModelComplex(Args &&...args);
			virtual ~Base_E_TurtbulentModelComplex() {};

		};
		/**
		 * @class E_TurtbulentModelRelaxationWithWall
		 * @brief .
		 */
		class E_TurtbulentModelWithWall
			: public Base_E_TurtbulentModelComplex<E_TurtbulentModelInner>
		{
		public:
			template <typename... Args>
			explicit E_TurtbulentModelWithWall(Args &&...args)
				: Base_E_TurtbulentModelComplex<E_TurtbulentModelInner>(std::forward<Args>(args)...) {};
			virtual ~E_TurtbulentModelWithWall() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};


    }
}
#endif // K_EPSILON_TURBULENT_MODEL_COMPLEX_H