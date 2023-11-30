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
		StdLargeVec<Real> turbu_mu_;
		StdLargeVec<Real> turbu_k_;
		StdLargeVec<Real> turbu_epsilon_;
		Real smoothing_length_;
		Real particle_spacing_min_;
		Real mu_;
		StdLargeVec<Real>& rho_;
		StdLargeVec<Vecd>& vel_;
		int dimension_;
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
		StdLargeVec<Matd>& velocity_gradient_;
		StdLargeVec<int>& is_near_wall_P1_; //** This is used to specially treat near wall region  *

		//**For test*
		StdLargeVec<Matd> velocity_gradient_wall;
	};
	//** Inner part *
	template <>
	class GetVelocityGradient<Inner<>> : public GetVelocityGradient<Base, FluidDataInner>
	{
	public:
		explicit GetVelocityGradient(BaseInnerRelation& inner_relation)
			: GetVelocityGradient<Base, FluidDataInner>(inner_relation) {};
		virtual ~GetVelocityGradient() {};
		void interaction(size_t index_i, Real dt = 0.0);
	};
	using GetVelocityGradientInner = GetVelocityGradient<Inner<>>;
	//** Interface part *
	//template <class InnerInteractionType>
	//using BaseGetVelocityGradientInner = GetVelocityGradient<InnerInteractionType>;
	
	//using GetVelocityGradientInner = BaseGetVelocityGradientInner<Inner<>>;
//=================================================================================================//
	template <typename... InteractionTypes>
	class TKEnergyAcc;

	template <class DataDelegationType>
	class TKEnergyAcc<Base, DataDelegationType>
		: public BaseTurtbulentModel<Base, DataDelegationType>
	{
	public:
		template <class BaseRelationType>
		explicit TKEnergyAcc(BaseRelationType& base_relation);
		virtual ~TKEnergyAcc() {};
	protected:
		StdLargeVec<Real>& turbu_k_;
		StdLargeVec<Vecd>& acc_prior_;
		StdLargeVec<Vecd>& pos_;
		StdLargeVec<int>& indicator_;
		StdLargeVec<Vecd> tke_acc_inner_, tke_acc_wall_;
		StdLargeVec<Vecd> test_k_grad_rslt_;
	};
	//** Inner part *
	template <>
	class TKEnergyAcc<Inner<>> : public TKEnergyAcc<Base, FluidDataInner>
	{
	public:
		explicit TKEnergyAcc(BaseInnerRelation& inner_relation);
		virtual ~TKEnergyAcc() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd>  &test_k_grad_rslt_;
	};
	//** Wall part *
	template <>
	class TKEnergyAcc<Contact<>> : public TKEnergyAcc<Base, FluidContactData>
	{
	public:
		explicit TKEnergyAcc(BaseContactRelation& contact_relation);
			//: TKEnergyAcc<Base, FluidContactData>(contact_relation) {};
		virtual ~TKEnergyAcc() {};
		void interaction(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Vecd>& test_k_grad_rslt_;
	};

	//** Interface part *
	template <class InnerInteractionType, class... ContactInteractionTypes>
	using BaseTKEnergyAccComplex = ComplexInteraction<TKEnergyAcc<InnerInteractionType, ContactInteractionTypes...>>;

	using TKEnergyAccComplex = BaseTKEnergyAccComplex<Inner<>, Contact<>>;
//=================================================================================================//
	template <typename... InteractionTypes>
	class TurbuViscousAcceleration;

	template <class DataDelegationType>
	class TurbuViscousAcceleration<DataDelegationType>: public ViscousAcceleration<DataDelegationType>
	{
	public:
		template <class BaseRelationType>
		explicit TurbuViscousAcceleration(BaseRelationType& base_relation);
		virtual ~TurbuViscousAcceleration() {};
	protected:
		StdLargeVec<Real>& turbu_mu_;
		StdLargeVec<Real>& wall_Y_plus_;
		StdLargeVec<Vecd>& velo_friction_;
		StdLargeVec<Vecd> visc_acc_inner_, visc_acc_wall_;
		StdLargeVec<Real>& y_p_;
	};

	//** Inner part *
	template <>
	class TurbuViscousAcceleration<Inner<>> : public TurbuViscousAcceleration<FluidDataInner>
	{
	public:
		explicit TurbuViscousAcceleration(BaseInnerRelation& inner_relation);
			//: TurbuViscousAcceleration<FluidDataInner>(inner_relation) {};
		virtual ~TurbuViscousAcceleration() {};
		void interaction(size_t index_i, Real dt = 0.0);
	};

	//** Wall part *
	using BaseTurbuViscousAccelerationWithWall = InteractionWithWall<TurbuViscousAcceleration>;
	template <>
	class TurbuViscousAcceleration<ContactWall<>> : public BaseTurbuViscousAccelerationWithWall
	{
	public:
		explicit TurbuViscousAcceleration(BaseContactRelation& wall_contact_relation);
			//: BaseTurbuViscousAccelerationWithWall(wall_contact_relation) {};
		virtual ~TurbuViscousAcceleration() {};
		void interaction(size_t index_i, Real dt = 0.0);
	};

	//** Interface part *
	using TurbulentViscousAccelerationWithWall = ComplexInteraction<TurbuViscousAcceleration<Inner<>, ContactWall<>>>;

//=================================================================================================//
}
}
#endif // K_EPSILON_TURBULENT_MODEL_H