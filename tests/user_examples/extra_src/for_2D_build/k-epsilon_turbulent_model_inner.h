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
 * @brief 	
 * @details     
 * @author Xiangyu Hu
 */

#ifndef K_EPSILON_TURBULENT_MODEL_INNER_H
#define K_EPSILON_TURBULENT_MODEL_INNER_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
    namespace fluid_dynamics
    {
		/**
		* @class BaseTurbuClosureCoeffInner
		* @brief  Some turbulent empirical parameters
		*/
		class BaseTurbuClosureCoeffInner
		{
		public:
			explicit BaseTurbuClosureCoeffInner();
			virtual ~BaseTurbuClosureCoeffInner() {};

		protected:
			Real Karman;
			Real turbu_const_E;
			Real C_mu;
			Real TurbulentIntensity;
			//** Closure coefficients for K *
			Real sigma_k;

			//** Closure coefficients for Epsilon *
			Real C_l, C_2;
			Real sigma_E;
		};

		/**
		 * @class BaseTurtbulentModelInner
		 * @brief BaseTurtbulentModelInner
		 */
		class BaseTurtbulentModelInner : public LocalDynamics, public FluidDataInner, public BaseTurbuClosureCoeffInner
		{
		public:
			explicit BaseTurtbulentModelInner(BaseInnerRelation& inner_relation);
			virtual ~BaseTurtbulentModelInner() {};

		protected:
			StdLargeVec<Real> turbu_mu_;
			StdLargeVec<Real> turbu_k_;
			StdLargeVec<Real> turbu_epsilon_;
			Real smoothing_length_;
			Real particle_spacing_min_;
			Real mu_;
			StdLargeVec<Real> & rho_;
			StdLargeVec<Vecd>& vel_;
			int dimension_;
		};

		/**
		 * @class K_TurtbulentModelInner
		 * @brief  K_TurtbulentModelInner
		 */
		class K_TurtbulentModelInner : public BaseTurtbulentModelInner
		{
		public:
			explicit K_TurtbulentModelInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values);
			virtual ~K_TurtbulentModelInner() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real> dk_dt_;
			StdLargeVec<Matd> velocity_gradient_;
			//StdLargeVec<Matd> B_;
			StdLargeVec<Real> k_production_;

			StdLargeVec<int> is_near_wall_P1_; //** This is used to specially treat near wall region  *
			Real turbu_k_initial_, turbu_ep_initial_, turbu_mu_initial_;

			//** for test */
			StdLargeVec<Real>  k_diffusion_, vel_x_;
			//StdLargeVec<Matd> velocity_gradient_wall;

		};

		/**
		 * @class E_TurtbulentModelInner
		 * @brief  E_TurtbulentModelInner
		 */
		class E_TurtbulentModelInner : public BaseTurtbulentModelInner
		{
		public:
			explicit E_TurtbulentModelInner(BaseInnerRelation& inner_relation);
			virtual ~E_TurtbulentModelInner() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real> dE_dt_, ep_production, ep_dissipation_, ep_diffusion_;
			StdLargeVec<Real>& turbu_mu_;
			StdLargeVec<Real>& turbu_k_;
			StdLargeVec<Real>& turbu_epsilon_;
			StdLargeVec<Real> &k_production_;
		};

		/**
		 * @class TurbuViscousAccInner
		 * @brief  the turbulent viscosity force induced acceleration
		 */
		class TurbulentEddyViscosity : public LocalDynamics, 
			public FluidDataSimple, public BaseTurbuClosureCoeffInner
		{
		public:
			explicit TurbulentEddyViscosity(SPHBody& sph_body);
			virtual ~TurbulentEddyViscosity() {};

			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real>& rho_;
			StdLargeVec<Real>& turbu_mu_;
			StdLargeVec<Real>& turbu_k_;
			StdLargeVec<Real>& turbu_epsilon_;
			StdLargeVec<Real>& wall_Y_plus_, & wall_Y_star_;
			Real mu_;
		};

		/**
		 * @class TurbulentAdvectionTimeStepSize
		 * @brief Computing the turbulent advection time step size
		 */
		class TurbulentAdvectionTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>,
			public FluidDataSimple
		{
		public:
			explicit TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL = 0.25);
			virtual ~TurbulentAdvectionTimeStepSize() {};
			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;
		protected:
			StdLargeVec<Vecd>& vel_;
			Real smoothing_length_min_;
			Real advectionCFL_;
			StdLargeVec<Real>& turbu_mu_;
			Fluid& fluid_;
		};

		/**
		 * @class   InflowTurbulentCondition
		 * @brief   Inflow boundary condition which imposes directly to a given velocity profile.
		 *          TargetVelocity gives the velocity profile along the inflow direction,
		 *          i.e. x direction in local frame.
		 */
		class InflowTurbulentCondition :public BaseFlowBoundaryCondition, public BaseTurbuClosureCoeffInner
		{
		public:
			explicit InflowTurbulentCondition(BodyPartByCell& body_part,
				Real CharacteristicLength, Real relaxation_rate=1.0);
			virtual ~InflowTurbulentCondition() {};
			void update(size_t index_i, Real dt = 0.0);
		protected:
			Real relaxation_rate_;
			StdLargeVec<Real>& turbu_k_;
			StdLargeVec<Real>& turbu_epsilon_;
			Real TurbulentLength_;
			Real CharacteristicLength_;

			virtual Real getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k);
			virtual Real getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E);
		};
    }
}
#endif // K_EPSILON_TURBULENT_MODEL_INNER_H