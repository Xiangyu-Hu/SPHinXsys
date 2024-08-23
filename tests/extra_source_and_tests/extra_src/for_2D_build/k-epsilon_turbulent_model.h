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
		Real turbulent_length_ratio_for_epsilon_inlet_;

		//** Start time for laminar law *
		Real start_time_laminar_;
		Real y_star_threshold_laminar_;
	};
//=================================================================================================//
	class WallFunction : public BaseTurbuClosureCoeff
	{
	public:
		explicit WallFunction() {};
		virtual ~WallFunction() {};

		Real get_dimensionless_velocity(Real y_star);
		Real get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity);
		Real get_distance_from_P_to_wall(Real y_p_constant);

		Real log_law_wall_functon(Real y_star);
		Real laminar_law_wall_functon(Real y_star);
		Real log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law);
		Real laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity);
	};
//=================================================================================================//
	template <typename... InteractionTypes>
	class GetVelocityGradient;

	template <class DataDelegationType>
	class GetVelocityGradient<DataDelegationType>
		: public LocalDynamics, public DataDelegationType
	{
	public:
		template <class BaseRelationType>
		explicit GetVelocityGradient(BaseRelationType& base_relation);
		virtual ~GetVelocityGradient() {};
	protected:
		StdLargeVec<Real> &Vol_;
		StdLargeVec<Vecd> &vel_,  &pos_;
		StdLargeVec<int> &is_near_wall_P1_; //** This is used to specially treat near wall region  *
		StdLargeVec<int> &is_near_wall_P2_;

		StdLargeVec<Matd> &velocity_gradient_;
		//**For test*
		StdLargeVec<Matd> &velocity_gradient_wall;
	};
	//** Inner part *
	template <>
	class GetVelocityGradient<Inner<>> : public GetVelocityGradient<DataDelegateInner>
	{
	public:
		explicit GetVelocityGradient(BaseInnerRelation& inner_relation);
		virtual ~GetVelocityGradient() {};
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Matd> &velocity_gradient_;
		StdLargeVec<Matd> &B_;
		StdLargeVec<Matd> &turbu_B_;
	};
	using GetVelocityGradientInner = GetVelocityGradient<Inner<>>;
	
	//** Updated Wall part *
	template <>
	class GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<GetVelocityGradient>
	{
	public:
		explicit GetVelocityGradient(BaseContactRelation& contact_relation);
		virtual ~GetVelocityGradient() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Matd> &velocity_gradient_;
	};
	
	//** Interface part *
	using GetVelocityGradientComplex = ComplexInteraction<GetVelocityGradient<Inner<>, Contact<Wall>>>;

	//using GetVelocityGradientComplex = BaseGetVelocityGradientComplex<Inner<>, Contact<>>;
//=================================================================================================//
	class TransferVelocityGradient : public LocalDynamics,
		public DataDelegateSimple
	{
	public:
		explicit TransferVelocityGradient(SPHBody& sph_body);
		virtual ~TransferVelocityGradient() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<int> &is_near_wall_P1_;
		StdLargeVec<Matd> &velocity_gradient_;
		StdLargeVec<Matd> &vel_grad_;
	};
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
		
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<Matd> &turbu_strain_rate_; //** temporary naming to distinguish the regular strain rate * 
		
		Real mu_, smoothing_length_, particle_spacing_min_;
		StdLargeVec<Real> &rho_,&Vol_;
		StdLargeVec<Vecd> &vel_;
		int dimension_;
	};
//=================================================================================================//
	/**
	 * @class K_TurtbulentModelInner
	 * @brief  K_TurtbulentModelInner
	 */
	class K_TurtbulentModelInner : public BaseTurtbulentModel<Base, DataDelegateInner>
	{
	public:
		explicit K_TurtbulentModelInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values, int is_extr_visc_dissipa = 0);
		virtual ~K_TurtbulentModelInner() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
		//void update_prior_turbulent_value();
	protected:
		StdLargeVec<Real> &dk_dt_;
		StdLargeVec<Real> &k_production_;

		StdLargeVec<int> &is_near_wall_P1_; //** This is used to specially treat near wall region  *
		StdLargeVec<Matd> &velocity_gradient_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<Matd> &turbu_strain_rate_;
        StdLargeVec<int> &is_extra_viscous_dissipation_;

		//** for test */
		StdLargeVec<int>  &turbu_indicator_;
		StdLargeVec<Real>  &k_diffusion_, &vel_x_;
	};
//=================================================================================================//
	/**
	 * @class E_TurtbulentModelInner
	 * @brief  E_TurtbulentModelInner
	 */
	class E_TurtbulentModelInner : public BaseTurtbulentModel<Base, DataDelegateInner>
	{
	public:
		explicit E_TurtbulentModelInner(BaseInnerRelation& inner_relation);
		virtual ~E_TurtbulentModelInner() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Real> &dE_dt_, &ep_production, &ep_dissipation_, &ep_diffusion_;
		
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		StdLargeVec<Real> &k_production_;
		StdLargeVec<int> &is_near_wall_P1_;
	};
//=================================================================================================//
	/**
	 * @class kOmegaSST_TurtbulentModelInner
	 * @brief  kOmegaSST_TurtbulentModelInner
	 */
	class kOmegaSST_kTransportEquationInner : public BaseTurtbulentModel<Base, DataDelegateInner>
	{
	public:
		explicit kOmegaSST_kTransportEquationInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values);
		virtual ~kOmegaSST_kTransportEquationInner() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
		void update_prior_turbulent_value();
	protected:

	};
//=================================================================================================//
	/**
	 * @class kOmegaSST_TurtbulentModelInner
	 * @brief  kOmegaSST_TurtbulentModelInner
	 */
	class kOmegaSST_omegaTransportEquationInner : public BaseTurtbulentModel<Base, DataDelegateInner>
	{
	public:
		explicit kOmegaSST_omegaTransportEquationInner(BaseInnerRelation& inner_relation);
		virtual ~kOmegaSST_omegaTransportEquationInner() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
		void update_prior_turbulent_value();
	protected:

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
		StdLargeVec<Vecd> &force_;
		StdLargeVec<Real> &mass_;
		StdLargeVec<int> &indicator_;
		StdLargeVec<Vecd> &pos_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Vecd> &test_k_grad_rslt_;

		//StdLargeVec<Vecd> &tke_acc_inner_, &tke_acc_wall_;
	};
	//** Inner part *
	template <>
	class TKEnergyForce<Inner<>> : public TKEnergyForce<Base, DataDelegateInner>
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
	class TKEnergyForce<Contact<>> : public TKEnergyForce<Base, DataDelegateContact>
	{
	public:
		explicit TKEnergyForce(BaseContactRelation& contact_relation);
			//: TKEnergyForce<Base, DataDelegateContact>(contact_relation) {};
		virtual ~TKEnergyForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd> &test_k_grad_rslt_;
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
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<Real> &wall_Y_plus_;
		StdLargeVec<Real> &wall_Y_star_;
		StdLargeVec<Vecd> &velo_friction_;
		StdLargeVec<Real> &y_p_;
		StdLargeVec<int> &is_near_wall_P2_;
		Real molecular_viscosity_;
		Real c0_;

		//** For test *
		//StdLargeVec<Matd> visc_direction_matrix_;
		//StdLargeVec<Vecd> &visc_acc_inner_, &visc_acc_wall_;

	};

	//** Inner part *
	template <>
	class TurbuViscousForce<Inner<>> : public TurbuViscousForce<DataDelegateInner>, public ForcePrior
	{
	public:
		explicit TurbuViscousForce(BaseInnerRelation& inner_relation);
		virtual ~TurbuViscousForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<int>  &turbu_indicator_;
		StdLargeVec<int> &is_extra_viscous_dissipation_;
		StdLargeVec<Matd> &B_;
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
		StdLargeVec<Matd> &B_;
	};

	//** Interface part *
	using TurbulentViscousForceWithWall = ComplexInteraction<TurbuViscousForce<Inner<>, Contact<Wall>>>;
//=================================================================================================//
	/**
	 * @class TurbuViscousAccInner
	 * @brief  the turbulent viscosity force induced acceleration
	 */
	class TurbulentEddyViscosity : public LocalDynamics,
		public DataDelegateSimple, public BaseTurbuClosureCoeff
	{
	public:
		explicit TurbulentEddyViscosity(SPHBody& sph_body);
		virtual ~TurbulentEddyViscosity() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Real> &rho_;
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		StdLargeVec<Real> &wall_Y_plus_,  &wall_Y_star_;
		Real mu_;
	};
//=================================================================================================//
	/**
	 * @class TurbulentAdvectionTimeStepSize
	 * @brief Computing the turbulent advection time step size
	 */
	class TurbulentAdvectionTimeStepSize : public LocalDynamicsReduce<ReduceMax>,
		public DataDelegateSimple
	{
	public:
		explicit TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL = 0.25);
		virtual ~TurbulentAdvectionTimeStepSize() {};
		Real reduce(size_t index_i, Real dt = 0.0);
		virtual Real outputResult(Real reduced_value) override;
	protected:
		StdLargeVec<Vecd> &vel_;
		Real smoothing_length_min_;
		Real speed_ref_turbu_, advectionCFL_;
		StdLargeVec<Real> &turbu_mu_;
		Fluid &fluid_;
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
		Real CharacteristicLength_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		Real TurbulentLength_;
		
		virtual Real getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k);
		virtual Real getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E);
	};
//=================================================================================================//
	class JudgeIsNearWall : public LocalDynamics, public DataDelegateContact,
		public BaseTurbuClosureCoeff
	{
	public:
		JudgeIsNearWall(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation);
		virtual ~JudgeIsNearWall() {};
		inline void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Real> &distance_to_dummy_interface_;
		StdLargeVec<Real> &distance_to_dummy_interface_up_average_;
		StdLargeVec<int> &is_near_wall_P1_;
		StdLargeVec<int> &is_near_wall_P2_;
		StdLargeVec<int> &index_nearest_;
		StdLargeVec<Vecd> &e_nearest_tau_, &e_nearest_normal_;

		StdLargeVec<Vecd> &pos_;
		int dimension_;
		Real fluid_particle_spacing_, wall_particle_spacing_;
		StdLargeVec<Vecd> &distance_from_wall_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec < StdLargeVec<Vecd>*>  contact_n_;
		
	};
//=================================================================================================//
	class StandardWallFunctionCorrection : 
		public LocalDynamics, public DataDelegateContact, public WallFunction
	{
	public:
		StandardWallFunctionCorrection(BaseInnerRelation& inner_relation,
			BaseContactRelation& contact_relation,Real y_p_constant);
		virtual ~StandardWallFunctionCorrection() {};
		inline void interaction(size_t index_i, Real dt = 0.0);

	protected:
		//Real offset_dist_;
		StdLargeVec<Real> &y_p_;
		StdLargeVec<Real> &wall_Y_plus_, &wall_Y_star_;
		StdLargeVec<Real> &velo_tan_;
		StdLargeVec<Vecd> &velo_friction_;

		StdLargeVec<Vecd> &vel_;
		StdLargeVec<Real> &rho_;
		Real molecular_viscosity_;
		StdLargeVec<Real> &turbu_k_;
		StdLargeVec<Real> &turbu_epsilon_;
		StdLargeVec<Real> &turbu_mu_;
		StdLargeVec<int> &is_near_wall_P1_;
		StdLargeVec<int> &is_near_wall_P2_;
		StdLargeVec<Matd> &velocity_gradient_;
		StdLargeVec<Real> &k_production_;
		StdLargeVec<Real> &distance_to_dummy_interface_;
		StdLargeVec<Real> &distance_to_dummy_interface_up_average_;
		StdLargeVec<int> &index_nearest;
		StdLargeVec<Vecd> &e_nearest_tau_;
		StdLargeVec<Vecd> &e_nearest_normal_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec < StdLargeVec<Vecd>*>  contact_n_;
	};
//=================================================================================================//
	class ConstrainNormalVelocityInRegionP : public LocalDynamics,
		public DataDelegateSimple
	{
	public:
		explicit ConstrainNormalVelocityInRegionP(SPHBody &sph_body);
		virtual ~ConstrainNormalVelocityInRegionP() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd> &vel_;
		StdLargeVec<int> &is_near_wall_P1_;
		StdLargeVec<Vecd> &e_nearest_normal_;
	};
//=================================================================================================//
	template <typename... InteractionTypes>
	class ExtraTransportForce;

	template <class DataDelegationType, class ParticleScope>
	class ExtraTransportForce<Base, DataDelegationType, ParticleScope>
		: public LocalDynamics, public DataDelegationType
	{
	public:
		template <class BaseRelationType>
		explicit ExtraTransportForce(BaseRelationType &base_relation);
		virtual ~ExtraTransportForce() {};
	protected:
	    StdLargeVec<Real> &rho_;
		StdLargeVec<Vecd> &vel_;
		StdLargeVec<Vecd> &zero_gradient_residue_;
		StdLargeVec<Matd> extra_transport_stress_;
		StdLargeVec<Vecd> extra_transport_vel_;
		ParticleScope within_scope_;
	};
	//** Inner part *
	template <class LimiterType, typename... CommonControlTypes>
	class ExtraTransportForce<Inner<LimiterType>, CommonControlTypes...> 
        : public ExtraTransportForce<Base, DataDelegateInner, CommonControlTypes...>
	{
	public:
		explicit ExtraTransportForce(BaseInnerRelation& inner_relation);
		virtual ~ExtraTransportForce() {};
		void initialization(size_t index_i, Real dt = 0.0);
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
        StdLargeVec<Real> &Vol_;
		const Real h_ref_;
		StdLargeVec<Matd> &extra_transport_stress_;
		StdLargeVec<Vecd> &extra_transport_vel_;
		LimiterType limiter_;
	};

	template <class LimiterType, class ParticleScope>
	using ExtraTransportForceInner =ExtraTransportForce<Inner<LimiterType>, ParticleScope>;

	//** Wall part *
	template <typename... CommonControlTypes>
	class ExtraTransportForce<Contact<Boundary>, CommonControlTypes...> 
	    : public ExtraTransportForce<Base, DataDelegateContact, CommonControlTypes...>
	{
	public:
		explicit ExtraTransportForce(BaseContactRelation& contact_relation);
		virtual ~ExtraTransportForce() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
        StdLargeVec<Matd> &extra_transport_stress_;
		StdLargeVec<Vecd> &extra_transport_vel_;
	    StdVec<StdLargeVec<Real> *> wall_Vol_;
	};

	//** Interface part *
    template <class LimiterType, typename... CommonControlTypes>
	using BaseExtraTransportForceComplex = ComplexInteraction<ExtraTransportForce<Inner<LimiterType>, Contact<Boundary>>, CommonControlTypes...>;
    
	template <class ParticleScope>
	using ExtraTransportForceComplex = BaseExtraTransportForceComplex<NoLimiter, ParticleScope>;

	template <class ParticleScope>
	using ExtraTransportForceLimitedComplex = BaseExtraTransportForceComplex<TruncatedLinear, ParticleScope>;
//=================================================================================================//
	class ConstrainVelocityAt_Y_Direction : public LocalDynamics,
		public DataDelegateSimple
	{
	public:
		explicit ConstrainVelocityAt_Y_Direction(SPHBody& sph_body, Real Length_channel);
		virtual ~ConstrainVelocityAt_Y_Direction() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd> &vel_;
		StdLargeVec<Vecd> &pos_;
		Real length_channel_;
	};
//=================================================================================================//
	class UpdateTurbulentPlugFlowIndicator : public LocalDynamics,
		public DataDelegateSimple
	{
	public:
		explicit UpdateTurbulentPlugFlowIndicator(SPHBody& sph_body, Real DH);
		virtual ~UpdateTurbulentPlugFlowIndicator() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<int> &turbu_plug_flow_indicator_;
		StdLargeVec<Vecd> &pos_;
		Real channel_width_;
	};
//=================================================================================================//
template <int INDICATOR>
class TurbulentIndicatedParticles : public WithinScope
{
    StdLargeVec<int> &indicator_;
    StdLargeVec<int> &turbu_plug_flow_indicator_;
  public:
    explicit TurbulentIndicatedParticles(BaseParticles *base_particles)
        : WithinScope(),
          indicator_(*base_particles->getVariableByName<int>("Indicator")),
		  turbu_plug_flow_indicator_(*base_particles->getVariableByName<int>("TurbulentPlugFlowIndicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] == INDICATOR && turbu_plug_flow_indicator_[index_i] == INDICATOR;
    };
};

using TurbulentPlugFlowParticles = TurbulentIndicatedParticles<0>;
//=================================================================================================//
template <typename... InteractionTypes>
class TurbulentLinearGradientCorrectionMatrix;

template <class DataDelegationType>
class TurbulentLinearGradientCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};

  protected:
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Matd> &turbu_B_;
	StdLargeVec<Matd> &B_;
};

template <>
class TurbulentLinearGradientCorrectionMatrix<Inner<>>
    : public TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>
{
    Real turbu_alpha_;

  public:
    explicit TurbulentLinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation), turbu_alpha_(alpha){};
    template <typename BodyRelationType, typename FirstArg>
    explicit TurbulentLinearGradientCorrectionMatrix(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : TurbulentLinearGradientCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using TurbulentLinearGradientCorrectionMatrixInner = TurbulentLinearGradientCorrectionMatrix<Inner<>>;

//=================================================================================================//
	class GetLimiterOfTransportVelocityCorrection : public LocalDynamics, public DataDelegateSimple
	{
	public:
		explicit GetLimiterOfTransportVelocityCorrection(SPHBody& sph_body, Real slope);
		virtual ~GetLimiterOfTransportVelocityCorrection() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		const Real h_ref_;
		StdLargeVec<Vecd> &zero_gradient_residue_;
		Real slope_;
		StdLargeVec<Real> &limiter_tvc_;
	};
//=================================================================================================//
template <class ParticleScope>
using TVC_Limited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, LinearGradientCorrection, ParticleScope>;
template <class ParticleScope>
using TVC_NoLimiter_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, LinearGradientCorrection, ParticleScope>;
//=================================================================================================//
class ModifiedTruncatedLinear : public Limiter
{
    Real ref_, slope_;

  public:
    ModifiedTruncatedLinear(Real ref, Real slope = 1000.0)
        : Limiter(), ref_(ref), slope_(slope){};
    Real operator()(Real measure)
    {
        Real measure_scale = measure * ref_;
        return SMIN(slope_ * measure_scale, Real(1));
    };
};
template <class ParticleScope>
using TVC_ModifiedLimited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, LinearGradientCorrection, ParticleScope>;
template <class ParticleScope>
using TVC_ModifiedLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, LinearGradientCorrectionWithBulkScope, ParticleScope>;
template <class ParticleScope>
using TVC_ModifiedLimited_withoutLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SingleResolution, ModifiedTruncatedLinear, NoKernelCorrection, ParticleScope>;
//=================================================================================================//
//*********************TESTING MODULES*********************
//=================================================================================================//
	/** Note this is a temporary treatment *
	* @class BaseGetTimeAverageData
	* @brief  BaseGetTimeAverageData
	*/
	//template <class DataDelegationType>
	class BaseGetTimeAverageData : public BaseTurtbulentModel<Base, DataDelegateInner>
	{
	public:
		explicit BaseGetTimeAverageData(BaseInnerRelation& inner_relation, int num_observer_points);
		virtual ~BaseGetTimeAverageData() {};

		void output_time_history_data(Real cutoff_time);
		void get_time_average_data(Real cutoff_time);
	protected:
		PltEngine plt_engine_;

		StdLargeVec<Vecd> &pos_;
		StdLargeVec<int> num_in_cell_;
		StdLargeVec<Real> &turbu_k_,  &turbu_mu_,   &turbu_epsilon_;
		//std::vector<std::vector<Real>>  data_sto_;
		StdLargeVec<std::vector<Real>> data_sto_, data_loaded_;
		StdLargeVec<Real>  data_time_aver_sto_;
		//ConcurrentVec<ConcurrentVec<Real>> data_sto_;
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
		public DataDelegateSimple, public BaseTurbuClosureCoeff
	{
	public:
		explicit ClearYPositionForTest(SPHBody& sph_body);
		virtual ~ClearYPositionForTest() {};

		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd> &pos_;
		StdLargeVec<Vecd> &vel_;
	};
//=================================================================================================//
}
}
#endif // K_EPSILON_TURBULENT_MODEL_H