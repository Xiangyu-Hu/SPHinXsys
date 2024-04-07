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
 * @file 	particle_dynamics_diffusion_reaction.h
 * @brief 	This is the particle dynamics applicable for all type bodies
 * 			TODO: there is an issue on applying corrected configuration for contact bodies.. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef EXTRA_TWO_PHASE_HEAT_TRANSFER_H
#define EXTRA_TWO_PHASE_HEAT_TRANSFER_H

#include "all_particle_dynamics.h"
#include "two_phase_heat_transfer_particles.h"
#include "diffusion_reaction.h"

namespace SPH
{
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	using TwoPhaseDiffusionReactionSimpleData =
		DataDelegateSimple<TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>>;

	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	using TwoPhaseDiffusionReactionInnerData =
		DataDelegateInner<TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>>;

	template <class BaseParticlesType, class BaseMaterialType,
			  class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES = 1>
	using TwoPhaseDiffusionReactionContactData =
		DataDelegateContact<TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>,
							TwoPhaseDiffusionReactionParticles<ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>,
							DataDelegateEmptyBase>;
	/**
	 * @class  DiffusionReactionInitialCondition
	 * @brief pure abstract class for initial conditions
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class TwoPhaseDiffusionReactionInitialCondition
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	public:
		explicit TwoPhaseDiffusionReactionInitialCondition(SPHBody &sph_body);
		virtual ~TwoPhaseDiffusionReactionInitialCondition(){};

	protected:
		StdLargeVec<Vecd> &pos_;
		StdVec<StdLargeVec<Real>> &species_n_;
		StdLargeVec<Real>& thermal_conductivity_;
	};

	/**
	 * @class GetDiffusionTimeStepSize
	 * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class TwoPhaseGetDiffusionTimeStepSize
		: public BaseDynamics<Real>,
		  public TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	public:
		explicit TwoPhaseGetDiffusionTimeStepSize(SPHBody &sph_body);
		virtual ~TwoPhaseGetDiffusionTimeStepSize(){};

		virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };
		virtual Real parallel_exec(Real dt = 0.0) override { return exec(dt); };

	protected:
		Real diff_time_step_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesInner
	 * @brief Compute the diffusion relaxation process of all species
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesInner
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionInnerData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
	protected:
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion *> species_diffusion_;
		StdVec<StdLargeVec<Real>> &species_n_;
		StdVec<StdLargeVec<Real>> &diffusion_dt_;

		void initializeDiffusionChangeRate(size_t particle_i);
		void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij);
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

	public:
		static const int number_of_species_ = NUM_SPECIES;
		typedef BaseParticlesType InnerBaseParticlesType;
		typedef BaseMaterialType InnerBaseMaterialType;
		typedef BaseInnerRelation BodyRelationType;
		DiffusionReaction<BaseMaterialType, NUM_SPECIES> &diffusion_reaction_material_;

		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesInner(BaseInnerRelation &inner_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesInner(){};
		void interaction(size_t index_i, Real dt = 0.0);
		void update(size_t index_i, Real dt = 0.0);
	};


	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES = 1>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesComplex
		: public TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>,
		public TwoPhaseDiffusionReactionContactData<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>
	{
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_;
		StdLargeVec<Real>& thermal_conductivity_;
		StdVec<StdLargeVec<Real>>& diffusion_dt_;
		StdLargeVec<Real>& external_diffusion_dt_;
		StdLargeVec<Real>& external_diffusion_dt_sum_;
		StdVec<StdVec<StdLargeVec<Real>>*> contact_species_n_;
		StdVec<StdLargeVec<Real>*> contact_thermal_conductivity_;
	protected:
		void initializeExternalDiffusionChangeRate(size_t particle_i);
		void getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd& e_ij,
			Real surface_area_ij, const StdVec<StdLargeVec<Real>>& species_n_k, const StdLargeVec<Real>& thermal_conductivity_k);
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		typedef ComplexRelation BodyRelationType;
		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesComplex(ComplexRelation& complex_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesComplex() {};
		void interaction(size_t index_i, Real dt = 0.0);

	};

	/**
	 * @class InitializationRK
	 * @brief initialization of a runge-kutta integration scheme
	 */
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES = 1>
	class TwoPhaseInitializationRK
		: public LocalDynamics,
		  public TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>
	{
		StdVec<BaseDiffusion *> species_diffusion_;
		StdVec<StdLargeVec<Real>> &species_n_, &species_s_;

		void initializeIntermediateValue(size_t particle_i);

	public:
		TwoPhaseInitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s);
		virtual ~TwoPhaseInitializationRK(){};

		void update(size_t index_i, Real dt = 0.0);
	};

	
	/**
	 * @class SecondStageRK2
	 * @brief the second stage of the 2nd-order Runge-Kutta scheme
	 */
	template <class FirstStageType>
	class TwoPhaseSecondStageRK2 : public FirstStageType
	{
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_;
		StdVec<StdLargeVec<Real>>& diffusion_dt_;

	protected:
		StdVec<StdLargeVec<Real>>& species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;
	public:
		TwoPhaseSecondStageRK2(typename FirstStageType::BodyRelationType& body_relation,
			StdVec<StdLargeVec<Real>>& species_s);
		virtual ~TwoPhaseSecondStageRK2() {};
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	template <class FirstStageType>
	class TwoPhaseRelaxationOfAllDiffusionSpeciesRK2 : public BaseDynamics<void>
	{
	protected:
		StdVec<BaseDiffusion*> species_diffusion_;
		/** Intermediate Value */
		StdVec<StdLargeVec<Real>> species_s_;

		SimpleDynamics<TwoPhaseInitializationRK<typename FirstStageType::InnerBaseParticlesType,
			typename FirstStageType::InnerBaseMaterialType,
			FirstStageType::number_of_species_>>
			rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<TwoPhaseSecondStageRK2<FirstStageType>> rk2_2nd_stage_;

	public:
		explicit TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType& body_relation);
		virtual ~TwoPhaseRelaxationOfAllDiffusionSpeciesRK2() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_H