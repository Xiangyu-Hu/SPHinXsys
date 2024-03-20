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

#ifndef K_EPSILON_TURBULENT_MODEL_H
#define K_EPSILON_TURBULENT_MODEL_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
	class BaseTurbuClosureCoeff
	{
	public:
		explicit BaseTurbuClosureCoeff();
		virtual ~BaseTurbuClosureCoeff() {};

	protected:
		Real Karman_;
		Real turbu_const_E_;
		Real C_mu_, C_mu_25_, C_mu_75_;
		Real turbulent_intensity_;

		//** Closure coefficients for K *
		Real sigma_k_;

		//** Closure coefficients for Epsilon *
		Real C_l_, C_2_;
		Real sigma_E_;
		Real mixing_length_for_epsilon_inlet_;

		//** Start time for laminar law *
		Real start_time_laminar_;
		Real y_star_threshold_;
	};
//=================================================================================================//
	class WallFunction : public BaseTurbuClosureCoeff
	{
	public:
		explicit WallFunction() {};
		virtual ~WallFunction() {};

		Real get_dimensionless_velocity(Real y_star);
		Real get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity);

		Real log_law_wall_functon(Real y_star);
		Real laminar_law_wall_functon(Real y_star);
		Real log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law);
		Real laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity);
	};
//=================================================================================================//
	template <typename... InteractionTypes>
	class GetVelocityGradient;

	template <class DataDelegationType>
	class GetVelocityGradient<Base, DataDelegationType>
		: public LocalDynamics, public DataDelegationType
	{
	public:
		template <class BaseRelationType>
		explicit GetVelocityGradient(BaseRelationType& base_relation);
		virtual ~GetVelocityGradient() {};
	protected:
		StdLargeVec<Vecd>& vel_, & pos_;
		StdLargeVec<int>& is_near_wall_P1_; //** This is used to specially treat near wall region  *

		StdLargeVec<Matd> velocity_gradient_;
		//**For test*
		StdLargeVec<Matd> velocity_gradient_wall;
	};
	//** Inner part *
	template <>
	class GetVelocityGradient<Inner<>> : public GetVelocityGradient<Base, FluidDataInner>
	{
	public:
		explicit GetVelocityGradient(BaseInnerRelation& inner_relation);
			//: GetVelocityGradient<Base, FluidDataInner>(inner_relation) {};
		virtual ~GetVelocityGradient() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Matd>& velocity_gradient_;
	};
	using GetVelocityGradientInner = GetVelocityGradient<Inner<>>;
	
	//** Wall part *
	template <>
	class GetVelocityGradient<Contact<>> : public GetVelocityGradient<Base, FSIContactData>
	{
	public:
		explicit GetVelocityGradient(BaseContactRelation& contact_relation);
			//: GetVelocityGradient<Base, FluidContactData>(contact_relation) {};
		virtual ~GetVelocityGradient() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Matd>& velocity_gradient_;
		StdVec<StdLargeVec<Vecd>*> wall_vel_ave_;
	};
	
	//** Interface part *
	template <class InnerInteractionType, class... ContactInteractionTypes>
	using BaseGetVelocityGradientComplex = ComplexInteraction<GetVelocityGradient<InnerInteractionType, ContactInteractionTypes...>>;

	using GetVelocityGradientComplex = BaseGetVelocityGradientComplex<Inner<>, Contact<>>;
//=================================================================================================//
	template <typename... T>
	class BaseTurtbulentModel;
	
	template <class DataDelegationType>
	class BaseTurtbulentModel<Base, DataDelegationType>
		: public LocalDynamics, public DataDelegationType, public BaseTurbuClosureCoeff
	{
	public:
		template <class BaseRelationType>
		explicit BaseTurtbulentModel(BaseRelationType& base_relation);
		virtual ~BaseTurtbulentModel() {};

	protected:
		
		StdLargeVec<Real> turbu_k_, turbu_k_prior_;
		StdLargeVec<Real> turbu_epsilon_, turbu_epsilon_prior_;
		StdLargeVec<Real> turbu_mu_;
		
		Real mu_, smoothing_length_, particle_spacing_min_;
		StdLargeVec<Real>& rho_;
		StdLargeVec<Vecd>& vel_;
		int dimension_;
	};
//=================================================================================================//
	/**
	 * @class K_TurtbulentModelInner
	 * @brief  K_TurtbulentModelInner
	 */
	class K_TurtbulentModelInner : public BaseTurtbulentModel<Base, FluidDataInner>
	{
	public:
		explicit K_TurtbulentModelInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values);
		virtual ~K_TurtbulentModelInner() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
		void update_prior_turbulent_value();
	protected:
		StdLargeVec<Real> dk_dt_;
		StdLargeVec<Real> k_production_, k_production_prior_;

		StdLargeVec<int>& is_near_wall_P1_; //** This is used to specially treat near wall region  *
		StdLargeVec<Matd>& velocity_gradient_;
		StdLargeVec<Real>& turbu_k_;
		StdLargeVec<Real>& turbu_k_prior_;
		StdLargeVec<Real>& turbu_epsilon_;
		StdLargeVec<Real>& turbu_epsilon_prior_;
		StdLargeVec<Real>& turbu_mu_;

		//** for test */
		StdLargeVec<Real>  k_diffusion_, vel_x_;
	};




//=================================================================================================//
	/**
	 * @class E_TurtbulentModelInner
	 * @brief  E_TurtbulentModelInner
	 */
	class E_TurtbulentModelInner : public BaseTurtbulentModel<Base, FluidDataInner>
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
		StdLargeVec<Real>& turbu_k_prior_;
		StdLargeVec<Real>& turbu_epsilon_;
		StdLargeVec<Real>& turbu_epsilon_prior_;
		StdLargeVec<Real>& k_production_;
		StdLargeVec<Real>& k_production_prior_;
		StdLargeVec<int>& is_near_wall_P1_;
	};
//=================================================================================================//
	template <typename... InteractionTypes>
	class TKEnergyForce;

	template <class DataDelegationType>
	class TKEnergyForce<Base, DataDelegationType>
		: public BaseTurtbulentModel<Base, DataDelegationType>
	{
	public:
		template <class BaseRelationType>
		explicit TKEnergyForce(BaseRelationType& base_relation);
		virtual ~TKEnergyForce() {};
	protected:
		StdLargeVec<Real>& turbu_k_;
		StdLargeVec<Vecd>& force_;
		StdLargeVec<Vecd>& pos_;
		StdLargeVec<Real>& mass_;
		StdLargeVec<int>& indicator_;
		StdLargeVec<Vecd> tke_acc_inner_, tke_acc_wall_;
		StdLargeVec<Vecd> test_k_grad_rslt_;
	};
	//** Inner part *
	template <>
	class TKEnergyForce<Inner<>> : public TKEnergyForce<Base, FluidDataInner>
	{
	public:
		explicit TKEnergyForce(BaseInnerRelation& inner_relation);
		virtual ~TKEnergyForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd>  &test_k_grad_rslt_;
	};
	//** Wall part *
	template <>
	class TKEnergyForce<Contact<>> : public TKEnergyForce<Base, FluidContactData>
	{
	public:
		explicit TKEnergyForce(BaseContactRelation& contact_relation);
			//: TKEnergyForce<Base, FluidContactData>(contact_relation) {};
		virtual ~TKEnergyForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd>& test_k_grad_rslt_;
	};

	//** Interface part *
	template <class InnerInteractionType, class... ContactInteractionTypes>
	using BaseTKEnergyForceComplex = ComplexInteraction<TKEnergyForce<InnerInteractionType, ContactInteractionTypes...>>;

	using TKEnergyForceComplex = BaseTKEnergyForceComplex<Inner<>, Contact<>>;
//=================================================================================================//
	template <typename... InteractionTypes>
	class TurbuViscousForce;

	template <class DataDelegationType>
	class TurbuViscousForce<DataDelegationType>: public ViscousForce<DataDelegationType>
	{
	public:
		template <class BaseRelationType>
		explicit TurbuViscousForce(BaseRelationType& base_relation);
		virtual ~TurbuViscousForce() {};
		
	protected:
		StdLargeVec<Real>& turbu_k_;
		StdLargeVec<Real>& turbu_mu_;
		StdLargeVec<Real>& wall_Y_plus_;
		StdLargeVec<Real>& wall_Y_star_;
		StdLargeVec<Vecd>& velo_friction_;
		StdLargeVec<Real>& y_p_;
		Real molecular_viscosity_;
		StdLargeVec<int>& is_near_wall_P2_;
		//** For test *
		//StdLargeVec<Matd> visc_direction_matrix_;
		StdLargeVec<Vecd> visc_acc_inner_, visc_acc_wall_;

	};

	//** Inner part *
	template <>
	class TurbuViscousForce<Inner<>> : public TurbuViscousForce<FluidDataInner>, public ForcePrior
	{
	public:
		explicit TurbuViscousForce(BaseInnerRelation& inner_relation);
		virtual ~TurbuViscousForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	};

	//** Wall part *
	using BaseTurbuViscousForceWithWall = InteractionWithWall<TurbuViscousForce>;
	template <>
	class TurbuViscousForce<Contact<Wall>> : public BaseTurbuViscousForceWithWall, public WallFunction
	{
	public:
		explicit TurbuViscousForce(BaseContactRelation& wall_contact_relation);
		virtual ~TurbuViscousForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		Real wall_particle_spacing_;
	};

	//** Interface part *
	using TurbulentViscousForceWithWall = ComplexInteraction<TurbuViscousForce<Inner<>, Contact<Wall>>>;
//=================================================================================================//
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
//=================================================================================================//
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
//=================================================================================================//
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
			Real CharacteristicLength, Real relaxation_rate = 1.0);
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
//=================================================================================================//
	class JudgeIsNearWall : public LocalDynamics, public FSIContactData,
		public BaseTurbuClosureCoeff
	{
	public:
		JudgeIsNearWall(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation, NearShapeSurface& near_surface);
		virtual ~JudgeIsNearWall() {};
		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Real> distance_to_dummy_interface_levelset_;
		StdLargeVec<Real> distance_to_dummy_interface_;
		StdLargeVec<Real> distance_to_dummy_interface_up_average_;
		StdLargeVec<int> is_near_wall_P1_;
		StdLargeVec<int> is_near_wall_P2_;
		StdLargeVec<int> index_nearest_;
		StdLargeVec<Vecd> e_nearest_tau_, e_nearest_normal_;

		LevelSetShape* level_set_shape_;
		StdLargeVec<Vecd> & pos_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		Real fluid_particle_spacing_, wall_particle_spacing_;
		StdVec < StdLargeVec<Vecd>*>  contact_n_;
		int dimension_;
	};
//=================================================================================================//
	class StandardWallFunctionCorrection : 
		public LocalDynamics, public FSIContactData, public WallFunction
	{
	public:
		StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation,Real y_p_constant);
		virtual ~StandardWallFunctionCorrection() {};
		inline void interaction(size_t index_i, Real dt = 0.0);

	protected:
		//Real offset_dist_;
		StdLargeVec<Real> y_p_;
		StdLargeVec<Real> wall_Y_plus_, wall_Y_star_;
		StdLargeVec<Real> velo_tan_;
		StdLargeVec<Vecd> velo_friction_;

		StdLargeVec<Real>& distance_to_dummy_interface_levelset_;
		StdLargeVec<Real>& distance_to_dummy_interface_;
		StdLargeVec<Real>& distance_to_dummy_interface_up_average_;
		StdLargeVec<Real>& turbu_k_;
		StdLargeVec<Real>& turbu_epsilon_;
		StdLargeVec<Real>& turbu_mu_;
		Real molecular_viscosity_;
		StdLargeVec<int>& is_near_wall_P1_;
		StdLargeVec<int>& is_near_wall_P2_;
		StdLargeVec<int>& index_nearest;
		StdLargeVec<Vecd>& vel_;
		StdLargeVec<Real>& rho_;
		StdLargeVec<Matd>& velocity_gradient_;
		StdLargeVec<Real>& k_production_;
		StdLargeVec<Vecd>& e_nearest_tau_;
		StdLargeVec<Vecd>& e_nearest_normal_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec < StdLargeVec<Vecd>*>  contact_n_;
	};
//=================================================================================================//
//*********************TESTING MODULES*********************
//=================================================================================================//
	/** Note this is a temporary treatment *
	* @class BaseGetTimeAverageData
	* @brief  BaseGetTimeAverageData
	*/
	//template <class DataDelegationType>
	class BaseGetTimeAverageData : public BaseTurtbulentModel<Base, FluidDataInner>
	{
	public:
		explicit BaseGetTimeAverageData(BaseInnerRelation& inner_relation, int num_observer_points);
		virtual ~BaseGetTimeAverageData() {};

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
		explicit GetTimeAverageCrossSectionData(BaseInnerRelation& inner_relation, int num_observer_points, const StdVec<Real>& bound_x, Real offset_dist_y = 0.0);
		virtual ~GetTimeAverageCrossSectionData() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		Real x_min_, x_max_;
		Real offset_dist_y_;
		StdVec<Real> monitor_cellcenter_y;
	};
	/** Note this is a temporary treatment *
	* @class GetTimeAverageCrossSectionData_Y
	* @brief  GetTimeAverageCrossSectionData_Y
	*/
	class GetTimeAverageCrossSectionData_Y : public BaseGetTimeAverageData
	{
	public:
		explicit GetTimeAverageCrossSectionData_Y(BaseInnerRelation& inner_relation, int num_observer_points, Real observe_x_ratio,
			const StdVec<Real>& bound_y, const StdVec<Real>& bound_x_f, const StdVec<Real>& bound_x_b);
		virtual ~GetTimeAverageCrossSectionData_Y() {};

		void update(size_t index_i, Real dt = 0.0);
		void output_monitor_x_coordinate();
	protected:
		Real observe_x_ratio_, observe_x_spacing_;
		StdVec<Real> bound_x_f_, bound_x_b_, bound_y_;
	};

//=================================================================================================//
	/**
	 * @class ClearYPositionForTest
	 * @brief  Test
	 */
	class ClearYPositionForTest : public LocalDynamics,
		public FluidDataSimple, public BaseTurbuClosureCoeff
	{
	public:
		explicit ClearYPositionForTest(SPHBody& sph_body);
		virtual ~ClearYPositionForTest() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd>& pos_;
		StdLargeVec<Vecd>& vel_;
	};
//=================================================================================================//
}
}
#endif // K_EPSILON_TURBULENT_MODEL_H