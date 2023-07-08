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

#ifndef K_EPSILON_TURBULENT_MODEL_INNER_H
#define K_EPSILON_TURBULENT_MODEL_INNER_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
    namespace fluid_dynamics
    {
		/**
		* @class BaseTurbuClosureCoeff
		* @brief  BaseTurbuClosureCoeff
		*/
		class BaseTurbuClosureCoeff
		{
		public:
			explicit BaseTurbuClosureCoeff();
			virtual ~BaseTurbuClosureCoeff() {};

		protected:
			Real Karman;
			Real turbu_const_E;
			Real C_mu;
			Real TurbulentIntensity;
			//K equation
			Real sigma_k;

			//closure coefficients for epsilon model
			Real C_l, C_2;
			Real sigma_E;

		};

		/**
		 * @class BaseTurtbulentModelInner
		 * @brief BaseTurtbulentModelInner
		 */
		class BaseTurtbulentModelInner : public LocalDynamics, public FluidDataInner, public BaseTurbuClosureCoeff
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
			explicit K_TurtbulentModelInner(BaseInnerRelation& inner_relation);
			virtual ~K_TurtbulentModelInner() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real> dk_dt_;
			StdLargeVec<Matd>  velocity_gradient;

			StdLargeVec<Real> k_production_;

			//** for test */
			StdLargeVec<Real> lap_k_, lap_k_term_, vel_x_;
			StdLargeVec<Matd>  velocity_gradient_inner, velocity_gradient_wall;
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
			StdLargeVec<Real> dE_dt_;
			StdLargeVec<Real>& turbu_mu_;
			StdLargeVec<Real>& turbu_k_;
			StdLargeVec<Real>& turbu_epsilon_;
			StdLargeVec<Real> &k_production_;
		};
		/**
		 * @class E_TurtbulentModelInner
		 * @brief  E_TurtbulentModelInner
		 */
		class TKEnergyAccInner : public BaseTurtbulentModelInner
		{
		public:
			explicit TKEnergyAccInner(BaseInnerRelation& inner_relation);
			virtual ~TKEnergyAccInner() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real>& turbu_k_;
			StdLargeVec<Vecd>& acc_prior_;

			StdLargeVec<Vecd>& pos_;
			StdLargeVec<int>& surface_indicator_;
			StdLargeVec<Vecd> tke_acc_inner_, tke_acc_wall_;
			StdLargeVec<Vecd> test_k_grad_rslt_;
		};

		/**
		 * @class TurbuViscousAccInner
		 * @brief  the turbulent viscosity force induced acceleration
		 */
		class TurbuViscousAccInner : public BaseViscousAccelerationInner
		{
		public:
			explicit TurbuViscousAccInner(BaseInnerRelation& inner_relation) ;
			virtual ~TurbuViscousAccInner() {};

			inline void interaction(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Real>& turbu_mu_;
			StdLargeVec<Real>& wall_Y_plus_;
			StdLargeVec<Vecd>& velo_friction_;
			StdLargeVec<Vecd> visc_acc_inner_, visc_acc_wall_;
			StdLargeVec<Real>& distance_to_wall_;
		};

		/**
		 * @class TurbuViscousAccInner
		 * @brief  the turbulent viscosity force induced acceleration
		 */
		class TurbulentEddyViscosity : public LocalDynamics, 
			public FluidDataSimple, public BaseTurbuClosureCoeff
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
		class InflowTurbulentCondition :public BaseFlowBoundaryCondition, public BaseTurbuClosureCoeff
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