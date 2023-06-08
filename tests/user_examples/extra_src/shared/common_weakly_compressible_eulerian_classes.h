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
  * @file 	common_weakly_compressible_eulerian_classes.h
  * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
  * @author	Zhentong Wang and Xiangyu Hu
  */
#ifndef COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H
#define COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H

#include "fluid_body.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "compressible_fluid.h"
#include "riemann_solver.h"
#include "fluid_dynamics_complex.h"

namespace SPH
{
	/**
	* @class EulerianWCTimeStepInitialization
	* @brief initialize a time step for a body.
	* including initialize particle acceleration
	* induced by viscous, gravity and other forces,
	* set the number of ghost particles into zero.
	*/
	class EulerianWCTimeStepInitialization : public TimeStepInitialization
	{
	public:
		EulerianWCTimeStepInitialization(SPHBody& sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
		virtual ~EulerianWCTimeStepInitialization() {};
		void update(size_t index_i, Real dt = 0.0);
	protected:
		StdLargeVec<Real>& rho_;
		StdLargeVec<Vecd>& pos_, & vel_;
		StdLargeVec<Vecd>& dmom_dt_prior_;
	};

	/**
	* @class EulerianWCAcousticTimeStepSize
	* @brief Computing the acoustic time step size
	*/
	class EulerianWCAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
	{
	public:
		explicit EulerianWCAcousticTimeStepSize(SPHBody& sph_body, Real CFL = 0.6) : AcousticTimeStepSize(sph_body, CFL) {};
		virtual ~EulerianWCAcousticTimeStepSize() {};
		virtual Real outputResult(Real reduced_value) override;
	};

	//----------------------------------------------------------------------
	//	Remann Solver classes.
	//----------------------------------------------------------------------
	/**
    * @struct AcousticRiemannSolverInEulerianMethod
    * @brief  Acoustic RiemannSolver for weakly-compressible flow in Eulerian method.
    */
	class AcousticRiemannSolverInEulerianMethod
	{
		Fluid& fluid_i_, & fluid_j_;

	public:
		AcousticRiemannSolverInEulerianMethod(Fluid& compressible_fluid_i, Fluid& compressible_fluid_j)
			: fluid_i_(compressible_fluid_i), fluid_j_(compressible_fluid_j) {};
		FluidStarState getInterfaceState(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
		{
			Real ul = -e_ij.dot(state_i.vel_);
			Real ur = -e_ij.dot(state_j.vel_);
			Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
			Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
			Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);

			Real p_star = (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ + rhol_cl * rhor_cr * (ul - ur)
				* SMIN(15.0 * SMAX((ul - ur) / clr, 0.0), 1.0)) / (rhol_cl + rhor_cr);
            Real u_star = (rhol_cl * ul + rhor_cr * ur + (state_i.p_ - state_j.p_) * pow(SMIN(15.0 * SMAX((ul - ur) / clr, 0.0), 1.0),2)) / (rhol_cl + rhor_cr);
			Vecd vel_star = (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_)
				- e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));

			FluidStarState interface_state(vel_star, p_star);
			interface_state.vel_ = vel_star;
			interface_state.p_ = p_star;

			return interface_state;
		};
	};

	//----------------------------------------------------------------------
	//	Viscous force inner and between two bodies
	//----------------------------------------------------------------------
	/**
	* @class EulerianViscousAccelerationInner
	* @brief  the viscosity force induced acceleration
	*/
	class EulerianViscousAccelerationInner : public fluid_dynamics::BaseViscousAccelerationInner
	{
	public:
		explicit EulerianViscousAccelerationInner(BaseInnerRelation& inner_relation);
		virtual ~EulerianViscousAccelerationInner() {};
		void interaction(size_t index_i, Real dt = 0.0);
		StdLargeVec<Vecd>& dmom_dt_prior_;
	};

	/**
	* @class InteractionWithWall
	* @brief  template class viscous acceleration with wall boundary
	*/
	template <class BaseIntegrationType>
	class InteractionWithWall : public BaseIntegrationType, public fluid_dynamics::FluidWallData
	{
	public:
		template <class BaseBodyRelationType>
		InteractionWithWall(BaseBodyRelationType& base_body_relation, BaseContactRelation& wall_contact_relation)
			: BaseIntegrationType(base_body_relation), fluid_dynamics::FluidWallData(wall_contact_relation)
		{
			if (&base_body_relation.getSPHBody() != &wall_contact_relation.getSPHBody())
			{
				std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (size_t k = 0; k != fluid_dynamics::FluidWallData::contact_particles_.size(); ++k)
			{
				Real rho_0_k = fluid_dynamics::FluidWallData::contact_bodies_[k]->base_material_->ReferenceDensity();
				wall_inv_rho0_.push_back(1.0 / rho_0_k);
				wall_vel_ave_.push_back(fluid_dynamics::FluidWallData::contact_particles_[k]->AverageVelocity());
				wall_acc_ave_.push_back(fluid_dynamics::FluidWallData::contact_particles_[k]->AverageAcceleration());
				wall_n_.push_back(&(fluid_dynamics::FluidWallData::contact_particles_[k]->n_));
			}
		};
		virtual ~InteractionWithWall() {};

	protected:
		StdVec<Real> wall_inv_rho0_;
		StdVec<StdLargeVec<Vecd>*> wall_vel_ave_, wall_acc_ave_, wall_n_;
	};

	/**
	* @class ViscousWithWall
	* @brief  template class viscous acceleration with wall boundary
	*/
	template <class BaseViscousAccelerationType>
	class ViscousWithWall : public InteractionWithWall<BaseViscousAccelerationType>
	{
	public:
		// template for different combination of constructing body relations
		template <class BaseBodyRelationType>
		ViscousWithWall(BaseBodyRelationType& base_body_relation, BaseContactRelation& wall_contact_relation)
			: InteractionWithWall<BaseViscousAccelerationType>(base_body_relation, wall_contact_relation) {};
		virtual ~ViscousWithWall() {};
		void interaction(size_t index_i, Real dt = 0.0)
		{
			BaseViscousAccelerationType::interaction(index_i, dt);

			Real rho_i = this->rho_[index_i];
			const Vecd& vel_i = this->vel_[index_i];

			Vecd acceleration = Vecd::Zero();
			Vecd vel_derivative = Vecd::Zero();
			for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_ave_k = *(this->wall_vel_ave_[k]);
				Neighborhood& contact_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];

					vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
					acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
				}
			}

			this->dmom_dt_prior_[index_i] += acceleration * rho_i;
		};;
	};

	/** template interface class for different pressure relaxation with wall schemes */
	template <class BaseViscousAccelerationType>
	class BaseViscousAccelerationWithWall : public BaseViscousAccelerationType
	{
	public:
		explicit BaseViscousAccelerationWithWall(ComplexRelation& fluid_wall_relation)
			: BaseViscousAccelerationType(fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()) {};
		BaseViscousAccelerationWithWall(BaseInnerRelation& fluid_inner_relation,
			BaseContactRelation& wall_contact_relation)
			: BaseViscousAccelerationType(fluid_inner_relation, wall_contact_relation) {};
		BaseViscousAccelerationWithWall(ComplexRelation& fluid_complex_relation,
			BaseContactRelation& wall_contact_relation)
			: BaseViscousAccelerationType(fluid_complex_relation, wall_contact_relation) {};
	};
	using ViscousAccelerationWithWall = BaseViscousAccelerationWithWall<ViscousWithWall<EulerianViscousAccelerationInner>>;

	/**
	* @class EulerianBaseIntegration
	* @brief Pure abstract base class for all Eulerian fluid relaxation schemes
	*/
	class EulerianBaseIntegration : public fluid_dynamics::BaseIntegration
	{
	public:
		explicit EulerianBaseIntegration(BaseInnerRelation& inner_relation);
		virtual ~EulerianBaseIntegration() {};

	protected:
		StdLargeVec<Real>& Vol_;
		StdLargeVec<Vecd>& mom_, & dmom_dt_, & dmom_dt_prior_;
	};

	/**
	* @class BaseIntegration1stHalf
	* @brief Template class for pressure relaxation scheme with the Riemann solver
	* as template variable
	*/
	template <class RiemannSolverType>
	class BaseIntegration1stHalf : public EulerianBaseIntegration
	{
	public:
		explicit BaseIntegration1stHalf(BaseInnerRelation& inner_relation)
			: EulerianBaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {};
		virtual ~BaseIntegration1stHalf() {};
		RiemannSolverType riemann_solver_;
		void initialization(size_t index_i, Real dt = 0.0)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = fluid_.getPressure(rho_[index_i]);
		};
		void interaction(size_t index_i, Real dt = 0.0)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Vecd momentum_change_rate = dmom_dt_prior_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

				momentum_change_rate -= 2.0 *
					((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
			}
			dmom_dt_[index_i] = momentum_change_rate;
		};
		void update(size_t index_i, Real dt = 0.0)
		{
			mom_[index_i] += dmom_dt_[index_i] * dt;
			vel_[index_i] = mom_[index_i] / rho_[index_i];
		};
	};
	/** define the mostly used pressure relaxation scheme using Riemann solver */
	using Integration1stHalfAcousticRiemann = BaseIntegration1stHalf<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseIntegration1stHalfWithWall
	* @brief  template class pressure relaxation scheme with wall boundary
	*/
	template <class BaseIntegration1stHalfType>
	class BaseIntegration1stHalfWithWall : public InteractionWithWall<BaseIntegration1stHalfType>
	{
	public:
		// template for different combination of constructing body relations
		template <class BaseBodyRelationType>
		BaseIntegration1stHalfWithWall(BaseBodyRelationType& base_body_relation, BaseContactRelation& wall_contact_relation)
			: InteractionWithWall<BaseIntegration1stHalfType>(base_body_relation, wall_contact_relation) {};
		explicit BaseIntegration1stHalfWithWall(ComplexRelation& fluid_wall_relation)
			: BaseIntegration1stHalfWithWall(fluid_wall_relation.getInnerRelation(),
				fluid_wall_relation.getContactRelation()) {};
		virtual ~BaseIntegration1stHalfWithWall() {};
		void interaction(size_t index_i, Real dt = 0.0)
		{
			BaseIntegration1stHalfType::interaction(index_i, dt);

			FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);

			Vecd momentum_change_rate = Vecd::Zero();
			for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
				Neighborhood& wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd& e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;
					FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
					FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

					momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
				}
			}
			this->dmom_dt_[index_i] += momentum_change_rate;
		};
	};
	using Integration1stHalfAcousticRiemannWithWall = BaseIntegration1stHalfWithWall<Integration1stHalfAcousticRiemann>;

	/**
	* @class BaseIntegration2ndHalf
	* @brief  Template density relaxation scheme with different Riemann solver
	*/
	template <class RiemannSolverType>
	class BaseIntegration2ndHalf : public EulerianBaseIntegration
	{
	public:
		explicit BaseIntegration2ndHalf(BaseInnerRelation& inner_relation)
			: EulerianBaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {};;
		virtual ~BaseIntegration2ndHalf() {};
		RiemannSolverType riemann_solver_;
		void interaction(size_t index_i, Real dt = 0.0)
		{
			FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
			Real density_change_rate = 0.0;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

				FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
				FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

				Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
				density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
			}
			drho_dt_[index_i] = density_change_rate;
		};
		void update(size_t index_i, Real dt = 0.0)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
		};
	};
	using Integration2ndHalfAcousticRiemann = BaseIntegration2ndHalf<AcousticRiemannSolverInEulerianMethod>;

	/**
	* @class BaseIntegration2ndHalfWithWall
	* @brief template density relaxation scheme with using  Riemann solver.
	*/
	template <class BaseIntegration2ndHalfType>
	class BaseIntegration2ndHalfWithWall : public InteractionWithWall<BaseIntegration2ndHalfType>
	{
	public:
		// template for different combination of constructing body relations
		template <class BaseBodyRelationType>
		BaseIntegration2ndHalfWithWall(BaseBodyRelationType& base_body_relation, BaseContactRelation& wall_contact_relation)
			: InteractionWithWall<BaseIntegration2ndHalfType>(base_body_relation, wall_contact_relation) {};
		explicit BaseIntegration2ndHalfWithWall(ComplexRelation& fluid_wall_relation)
			: BaseIntegration2ndHalfWithWall(fluid_wall_relation.getInnerRelation(),
				fluid_wall_relation.getContactRelation()) {};
		virtual ~BaseIntegration2ndHalfWithWall() {};
		void interaction(size_t index_i, Real dt = 0.0)
		{
			BaseIntegration2ndHalfType::interaction(index_i, dt);

			FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
			Real density_change_rate = 0.0;
			for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& n_k = *(this->wall_n_[k]);
				Neighborhood& wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
				for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
				{
					size_t index_j = wall_neighborhood.j_[n];
					Vecd& e_ij = wall_neighborhood.e_ij_[n];
					Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

					Vecd vel_in_wall = -state_i.vel_;
					Real p_in_wall = state_i.p_;
					Real rho_in_wall = state_i.rho_;

					FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
					FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
					Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

					density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
				}
			}
			this->drho_dt_[index_i] += density_change_rate;
		};
	};
	using Integration2ndHalfAcousticRiemannWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalfAcousticRiemann>;

	//----------------------------------------------------------------------
	//	Non-Reflective Boundary
	//----------------------------------------------------------------------
	class NonReflectiveBoundaryVariableCorrection : public LocalDynamics, public DataDelegateInner<FluidParticles>
	{
	public:
		NonReflectiveBoundaryVariableCorrection(BaseInnerRelation& inner_relation);
		virtual ~NonReflectiveBoundaryVariableCorrection() {};
		void initialization(size_t index_i, Real dt = 0.0);
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	protected:
		Fluid& fluid_;
		Real rho_farfield_, sound_speed_;
		Vecd vel_farfield_;
		StdLargeVec<Real>& rho_, & p_, & Vol_;
		StdLargeVec<Vecd>& vel_, & mom_, & pos_;
		StdLargeVec<Vecd> n_;
		StdLargeVec<Real> inner_weight_summation_, rho_average_, vel_normal_average_;
		StdLargeVec<Vecd> vel_tangential_average_, vel_average_;
		StdLargeVec<int>& surface_indicator_;
		StdLargeVec<int> surface_inner_particle_indicator_;
	};
}
#endif // COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H