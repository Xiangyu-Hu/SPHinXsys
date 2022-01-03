/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file fluid_dynamics_multi_phase.h
* @brief Here, we define the algorithm classes for the dynamics involving multiple fluids.   
* @author	Chi ZHang and Xiangyu Hu
*/


#ifndef FLUID_DYNAMICS_MULTI_PHASE_H
#define FLUID_DYNAMICS_MULTI_PHASE_H

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
	namespace fluid_dynamics
	{
		typedef DataDelegateContact<FluidBody, FluidParticles, Fluid,
			FluidBody, FluidParticles, Fluid, DataDelegateEmptyBase> MultiPhaseContactData;
		typedef DataDelegateContact<FluidBody, FluidParticles, Fluid, FluidBody, FluidParticles, Fluid> MultiPhaseData;
		/**
		 * @class ViscousAccelerationMultiPhase
		 * @brief  the viscosity force induced acceleration
		 */
		class ViscousAccelerationMultiPhase : public ViscousAccelerationInner, public MultiPhaseContactData
		{
		public:
			ViscousAccelerationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation);
			explicit ViscousAccelerationMultiPhase(ComplexBodyRelation &complex_relation);
			virtual ~ViscousAccelerationMultiPhase() {};
		protected:
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_n_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using ViscousAccelerationMultiPhaseWithWall = 
			BaseViscousAccelerationWithWall<ViscousWithWall<ViscousAccelerationMultiPhase>>;

		/**
		* @class ViscousAccelerationMultiPhase
		* @brief Abstract base class for general multiphase fluid dynamics
		*/
		template<class RelaxationInnerType>
		class RelaxationMultiPhase : public RelaxationInnerType, public MultiPhaseContactData
		{
		public:
			RelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation);
			virtual ~RelaxationMultiPhase() {};
		protected:
			StdVec<StdLargeVec<Real>*> contact_Vol_, contact_p_, contact_rho_n_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_n_;
		};

		/**
		 * @class BasePressureRelaxationMultiPhase
		 * @brief  template class for multiphase pressure relaxation scheme
		 */
		template<class PressureRelaxationInnerType>
		class BasePressureRelaxationMultiPhase : public RelaxationMultiPhase<PressureRelaxationInnerType>
		{
		public:
			BasePressureRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation);
			explicit BasePressureRelaxationMultiPhase(ComplexBodyRelation &complex_relation);
			virtual ~BasePressureRelaxationMultiPhase() {};
		protected:
			using CurrentRiemannSolver = decltype(PressureRelaxationInnerType::riemann_solver_);
			StdVec<CurrentRiemannSolver> riemann_solvers_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};
		using MultiPhasePressureRelaxation = BasePressureRelaxationMultiPhase<PressureRelaxationInner>;
		using MultiPhasePressureRelaxationRiemann = BasePressureRelaxationMultiPhase<PressureRelaxationRiemannInner>;

		using MultiPhasePressureRelaxationWithWall = 
			BasePressureRelaxationWithWall<PressureRelaxation<MultiPhasePressureRelaxation>>;
		using MultiPhasePressureRelaxationRiemannWithWall = 
			BasePressureRelaxationWithWall<PressureRelaxation<MultiPhasePressureRelaxationRiemann>>;
		using ExtendMultiPhasePressureRelaxationRiemannWithWall = 
			ExtendPressureRelaxationWithWall<ExtendPressureRelaxation<MultiPhasePressureRelaxationRiemann>>;


		/**
		 * @class BaseDensityRelaxationMultiPhase
		 * @brief  template class pressure relaxation scheme with wall boundary
		 */
		template<class DensityRelaxationInnerType>
		class BaseDensityRelaxationMultiPhase : public RelaxationMultiPhase<DensityRelaxationInnerType>
		{
		public:
			BaseDensityRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation);
			explicit BaseDensityRelaxationMultiPhase(ComplexBodyRelation &complex_relation);
			virtual ~BaseDensityRelaxationMultiPhase() {};
		protected:
			using CurrentRiemannSolver = decltype(DensityRelaxationInnerType::riemann_solver_);
			StdVec<CurrentRiemannSolver> riemann_solvers_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		using MultiPhaseDensityRelaxation = BaseDensityRelaxationMultiPhase<DensityRelaxationInner>;
		using MultiPhaseDensityRelaxationRiemann = BaseDensityRelaxationMultiPhase<DensityRelaxationRiemannInner>;
		using MultiPhaseDensityRelaxationWithWall = BaseDensityRelaxationMultiPhase<MultiPhaseDensityRelaxation>;
		using MultiPhaseDensityRelaxationRiemannWithWall = BaseDensityRelaxationWithWall<MultiPhaseDensityRelaxationRiemann>;

		/**
		 * @class MultiPhaseColorFunctionGradient
		 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
		 */
		class MultiPhaseColorFunctionGradient : public InteractionDynamics, public MultiPhaseData
		{
		public:
			explicit MultiPhaseColorFunctionGradient(BaseBodyRelationContact &contact_relation);
			virtual ~MultiPhaseColorFunctionGradient() {};
		protected:
			Real rho0_;
			StdVec<Real> contact_rho0_;
			StdLargeVec<Real>& Vol_, & pos_div_;
			StdLargeVec<int>& surface_indicator_;
			StdLargeVec<Vecd> color_grad_, surface_norm_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //FLUID_DYNAMICS_MULTI_PHASE_H