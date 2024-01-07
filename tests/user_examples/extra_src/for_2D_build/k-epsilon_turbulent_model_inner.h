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

		/** Note this is a temporary treatment *
		* @class BaseGetTimeAverageData
		* @brief  BaseGetTimeAverageData
		*/
		class BaseGetTimeAverageData : public BaseTurtbulentModelInner
		{
		public:
			explicit BaseGetTimeAverageData(BaseInnerRelation& inner_relation, int num_observer_points);
			virtual ~BaseGetTimeAverageData() {};

			//void update(size_t index_i, Real dt = 0.0);
			void output_time_history_data(Real cutoff_time);
			void get_time_average_data(Real cutoff_time);
		protected:
			PltEngine plt_engine_;

			StdLargeVec<Vecd>& pos_;
			StdLargeVec<Real>& turbu_mu_, & turbu_k_, & turbu_epsilon_;
			//std::vector<std::vector<Real>>  data_sto_;
			StdLargeVec<std::vector<Real>> data_sto_, data_loaded_;
			StdLargeVec<Real>  data_time_aver_sto_;
			//ConcurrentVec<ConcurrentVec<Real>> data_sto_;
			StdLargeVec<int> num_in_cell_;
			int num_cell, num_data;
			StdLargeVec<std::string> file_name_;
			std::string file_path_output_, file_path_input_;
		};

		/** Note this is a temporary treatment *
		* @class GetTimeAverageCrossSectionData
		* @brief  GetTimeAverageCrossSectionData
		*/
		class GetTimeAverageCrossSectionData : public BaseGetTimeAverageData
		{
		public:
			explicit GetTimeAverageCrossSectionData(BaseInnerRelation& inner_relation,int num_observer_points, const StdVec<Real>& bound_x, Real offset_dist_y = 0.0);
			virtual ~GetTimeAverageCrossSectionData() {};

			void update(size_t index_i, Real dt = 0.0);
		protected:
			Real x_min_, x_max_;
			Real offset_dist_y_;
			StdVec<Real> monitor_cellcenter_y;
		};
		/** Note this is a temporary treatment *
		* @class GetTimeAverageCenterLineData
		* @brief  GetTimeAverageCenterLineData
		*/
		class GetTimeAverageCenterLineData : public BaseGetTimeAverageData
		{
		public:
			explicit GetTimeAverageCenterLineData(BaseInnerRelation& inner_relation, int num_observer_points,Real observe_x_ratio, 
				const StdVec<Real>& bound_y, const StdVec<Real>& bound_x_f, const StdVec<Real>& bound_x_b);
			virtual ~GetTimeAverageCenterLineData() {};

			void update(size_t index_i, Real dt = 0.0);
			void output_monitor_x_coordinate();
		protected:
			Real observe_x_ratio_, observe_x_spacing_;
			StdVec<Real> bound_x_f_, bound_x_b_, bound_y_;
		};

		/**
		 * @class GetVelocityGradientInner
		 * @brief  GetVelocityGradientInner
		 */
		//class GetVelocityGradientInner : public LocalDynamics, public FluidDataInner
		//{
		//public:
		//	explicit GetVelocityGradientInner(BaseInnerRelation& inner_relation);
		//	virtual ~GetVelocityGradientInner() {};

		//	inline void interaction(size_t index_i, Real dt = 0.0);
		//protected:
		//	StdLargeVec<Vecd>& vel_, & pos_;
		//	StdLargeVec<Matd>& velocity_gradient_;
		//	StdLargeVec<int>& is_near_wall_P1_; //** This is used to specially treat near wall region  *

		//	//**For test*
		//	StdLargeVec<Matd> velocity_gradient_wall;
		//};

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
		 * @class E_TurtbulentModelInner
		 * @brief  E_TurtbulentModelInner
		 */
		//class TKEnergyAccInner : public BaseTurtbulentModelInner
		//{
		//public:
		//	explicit TKEnergyAccInner(BaseInnerRelation& inner_relation);
		//	virtual ~TKEnergyAccInner() {};

		//	inline void interaction(size_t index_i, Real dt = 0.0);
		//protected:
		//	StdLargeVec<Real>& turbu_k_;
		//	StdLargeVec<Vecd>& acc_prior_;
		//	StdLargeVec<Vecd>& pos_;
		//	StdLargeVec<int>& indicator_;
		//	StdLargeVec<Vecd> tke_acc_inner_, tke_acc_wall_;
		//	StdLargeVec<Vecd> test_k_grad_rslt_;
		//};

		/**
		 * @class TurbuViscousAccInner
		 * @brief  the turbulent viscosity force induced acceleration
		 */
		//class TurbuViscousAccInner : public BaseViscousAccelerationInner, public BaseTurbuClosureCoeffInner
		//{
		//public:
		//	explicit TurbuViscousAccInner(BaseInnerRelation& inner_relation) ;
		//	virtual ~TurbuViscousAccInner() {};

		//	inline void interaction(size_t index_i, Real dt = 0.0);
		//protected:
		//	StdLargeVec<Real>& turbu_mu_;
		//	StdLargeVec<Real>& wall_Y_plus_;
		//	StdLargeVec<Vecd>& velo_friction_;
		//	StdLargeVec<Vecd> visc_acc_inner_, visc_acc_wall_;
		//	StdLargeVec<Real>& y_p_;
		//};

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