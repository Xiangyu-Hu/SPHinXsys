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
 * @brief 	Here, we define the static confinement boundary condition classes for fluid dynamics.
 * @details     This boundary condition is based on Level-set filed.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef K_EPSILON_TURBULENT_MODEL_H
#define K_EPSILON_TURBULENT_MODEL_H

#include "fluid_dynamics_inner.h"
#include "fluid_boundary.h"
#include "fluid_dynamics_complex.h"
#include <mutex>

namespace SPH
{
    namespace fluid_dynamics
    {
		/**
		* @class BaseTurbulentClosureCoefficient
		* @brief  BaseTurbulentClosureCoefficient
		*/
		class BaseTurbulentClosureCoefficient
		{
		public:
			explicit BaseTurbulentClosureCoefficient();
			virtual ~BaseTurbulentClosureCoefficient() {};

		protected:
			Real Karman;
			Real turbu_const_E;
			Real C_mu;
			Real TurbulentIntensity;
			Real TurbulentLength;

			//K equation
			Real sigma_k;

			//closure coefficients for epsilon model
			Real C_l, C_2;
			Real sigma_E;


		};

		/**
		* @class BaseTurtbulentData
		* @brief  BaseTurtbulentData
		*/
		class BaseTurtbulentData : public FluidDataInner, public BaseTurbulentClosureCoefficient
		{
		public:
			explicit BaseTurtbulentData(BaseInnerRelation &inner_relation);
			virtual ~BaseTurtbulentData() {};

		protected:

			Real particle_spacing_min_;
			Real mu_;
			StdLargeVec<Real>& Vol_, & rho_, & p_;
			StdLargeVec<Vecd>& vel_, & acc_prior_;
			int dimension_;
		};

		/**
		* @class BaseTurtbulentModel
		* @brief  BaseTurtbulentModel
		*/
		class BaseTurtbulentModel : public BaseTurtbulentData, public LocalDynamics
		{
		public:
			explicit BaseTurtbulentModel(BaseInnerRelation &inner_relation);
			virtual ~BaseTurtbulentModel() {};

		protected:
			StdLargeVec<Real> turbu_mu_;
			StdLargeVec<Real> turbu_k_;
			StdLargeVec<Real> turbu_epsilon_;
			Real ProductionTermBase;
			Real smoothing_length_;
			StdLargeVec<Matd>& grad_calculated_;
			StdLargeVec<Matd> Rij;
			StdLargeVec<Matd> grad_vel_ij, transpose_grad_vel_ij;
		};



		/**
		* @class K_TutbulentModel
		* @brief  K_TutbulentModel
		*/
		class K_TurtbulentModelInner : public BaseTurtbulentModel
		{
		public:
			explicit K_TurtbulentModelInner(BaseInnerRelation &inner_relation);
			virtual ~K_TurtbulentModelInner() {};

		protected:
			StdLargeVec<Real> dk_dt_;
			StdLargeVec<Real> production_k_;
			StdLargeVec<Real> lap_k_, lap_k_term_;
			StdLargeVec<Real> temp_dW_, temp_rij_;

			//only for convinience for output U
			StdLargeVec<Real> vel_x_n_;

			inline void interaction(size_t index_i, Real dt = 0.0);

			void update(size_t index_i, Real dt = 0.0);
		};


		/**
		 * @class K_TurtbulentModelRelaxationWithWall
		 * @brief .
		 * The
		 */
		template <class BaseIntegration2ndHalfType>
		class Base_K_TurtbulentModelWithWall : public InteractionWithWall<BaseIntegration2ndHalfType>
		{
		public:
			template <typename... Args>
			Base_K_TurtbulentModelWithWall(Args &&...args)
				: InteractionWithWall<BaseIntegration2ndHalfType>(std::forward<Args>(args)...) {};
			virtual ~Base_K_TurtbulentModelWithWall() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		};

        using K_TurtbulentModelRelaxationWithWall = Base_K_TurtbulentModelWithWall<K_TurtbulentModelInner>;


    }
}
#endif // FLUID_BOUNDARY_STATIC_COFINEMENT_H