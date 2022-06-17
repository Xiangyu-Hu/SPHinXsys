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
 * @file 	electro_physiology.h
 * @brief 	In is file, we declaim the dynamics relevant to electrophysiology,
 * including diffusion, reaction and muscle activation.
 * @author 	Chi Zhang and Xiangyu Hu
 */

#ifndef ELECTRO_PHYSIOLOGY_H
#define ELECTRO_PHYSIOLOGY_H

#include "particle_dynamics_diffusion_reaction.h"
#include "solid_particles.h"

namespace SPH
{
	class ElectroPhysiologyReaction : public BaseReactionModel
	{
	protected:
		Real k_a_;
		size_t voltage_;
		size_t gate_variable_;
		size_t active_contraction_stress_;

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) = 0;
		virtual Real getProductionActiveContractionStress(StdVec<StdLargeVec<Real>> &species, size_t particle_i);
		virtual Real getLossRateActiveContractionStress(StdVec<StdLargeVec<Real>> &species, size_t particle_i);

	public:
		explicit ElectroPhysiologyReaction(Real k_a)
			: BaseReactionModel({"Voltage", "GateVariable", "ActiveContractionStress"}),
			  k_a_(k_a), voltage_(species_indexes_map_["Voltage"]),
			  gate_variable_(species_indexes_map_["GateVariable"]),
			  active_contraction_stress_(species_indexes_map_["ActiveContractionStress"])
		{
			reaction_model_ = "ElectroPhysiologyReaction";
			initializeElectroPhysiologyReaction();
		};
		virtual ~ElectroPhysiologyReaction(){};

		void initializeElectroPhysiologyReaction();
	};

	/**
	 * @class AlievPanfilowModel
	 * @brief The simplest Electrophysiology Reaction model,
	 * which reduces the complex of array of ion currents to two variables that
	 * describe excitation and recovery.
	 */
	class AlievPanfilowModel : public ElectroPhysiologyReaction
	{
	protected:
		/** Parameters for two variable cell model. */
		Real k_, a_, b_, mu_1_, mu_2_, epsilon_, c_m_;

		virtual Real getProductionRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getLossRateIonicCurrent(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getProductionRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;
		virtual Real getLossRateGateVariable(StdVec<StdLargeVec<Real>> &species, size_t particle_i) override;

	public:
		explicit AlievPanfilowModel(Real k_a, Real c_m, Real k, Real a, Real b, Real mu_1, Real mu_2, Real epsilon)
			: ElectroPhysiologyReaction(k_a), k_(k), a_(a), b_(b), mu_1_(mu_1), mu_2_(mu_2),
			  epsilon_(epsilon), c_m_(c_m)
		{
			reaction_model_ = "AlievPanfilowModel";
		};
		virtual ~AlievPanfilowModel(){};
	};

	/**
	 * @class MonoFieldElectroPhysiology
	 * @brief material class for electro_physiology.
	 */
	class MonoFieldElectroPhysiology : public DiffusionReaction<Solid>
	{
	public:
		explicit MonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
											Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
		virtual ~MonoFieldElectroPhysiology(){};
	};

	/**
	 * @class LocalMonoFieldElectroPhysiology
	 * @brief material class for electro_physiology with locally oriented fibers.
	 */
	class LocalMonoFieldElectroPhysiology : public DiffusionReaction<Solid>
	{
	public:
		explicit LocalMonoFieldElectroPhysiology(ElectroPhysiologyReaction &electro_physiology_reaction,
												 Real diff_cf, Real bias_diff_cf, Vecd bias_direction);
		virtual ~LocalMonoFieldElectroPhysiology(){};

		virtual void readFromXmlForLocalParameters(const std::string &filefullpath) override;
	};

	/**
	 * @class ElectroPhysiologyParticles
	 * @brief A group of particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyParticles
		: public DiffusionReactionParticles<SolidParticles>
	{
	public:
		ElectroPhysiologyParticles(SPHBody &sph_body, DiffusionReaction<Solid> *diffusion_reaction_material);
		virtual ~ElectroPhysiologyParticles(){};
		virtual ElectroPhysiologyParticles *ThisObjectPtr() override { return this; };
	};
	/**
	 * @class ElectroPhysiologyReducedParticles
	 * @brief A group of reduced particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyReducedParticles
		: public DiffusionReactionParticles<SolidParticles>
	{
	public:
		/** Constructor. */
		ElectroPhysiologyReducedParticles(SPHBody &sph_body, DiffusionReaction<Solid> *diffusion_reaction_material);
		/** Destructor. */
		virtual ~ElectroPhysiologyReducedParticles(){};
		virtual ElectroPhysiologyReducedParticles *ThisObjectPtr() override { return this; };

		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd &e_ij) override
		{
			return dW_ij * e_ij;
		};
	};

	namespace electro_physiology
	{
		typedef DiffusionReactionSimpleData<RealBody, SolidParticles, Solid> ElectroPhysiologyDataDelegateSimple;
		typedef DiffusionReactionInnerData<RealBody, SolidParticles, Solid> ElectroPhysiologyDataDelegateInner;
		/**
		 * @class ElectroPhysiologyInitialCondition
		 * @brief  set initial condition for a muscle body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElectroPhysiologyInitialCondition : public ParticleDynamicsSimple,
												  public ElectroPhysiologyDataDelegateSimple
		{
		public:
			explicit ElectroPhysiologyInitialCondition(RealBody &real_body);
			virtual ~ElectroPhysiologyInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_;
			StdVec<StdLargeVec<Real>> &species_n_;
		};
		/**
		 * @class GetElectroPhysiologyTimeStepSize
		 * @brief Computing the time step size from diffusion criteria
		 */
		class GetElectroPhysiologyTimeStepSize
			: public GetDiffusionTimeStepSize<RealBody, SolidParticles, Solid>
		{
		public:
			explicit GetElectroPhysiologyTimeStepSize(RealBody &real_body)
				: GetDiffusionTimeStepSize<RealBody, SolidParticles, Solid>(real_body){};
			virtual ~GetElectroPhysiologyTimeStepSize(){};
		};
		/**
		 * @class ElectroPhysiologyDiffusionRelaxationInner
		 * @brief Compute the diffusion relaxation process
		 */
		class ElectroPhysiologyDiffusionRelaxationInner
			: public RelaxationOfAllDiffusionSpeciesRK2<
				  RelaxationOfAllDiffussionSpeciesInner<RealBody, SolidParticles, Solid>>
		{
		public:
			explicit ElectroPhysiologyDiffusionRelaxationInner(BaseBodyRelationInner &inner_relation)
				: RelaxationOfAllDiffusionSpeciesRK2(inner_relation){};
			virtual ~ElectroPhysiologyDiffusionRelaxationInner(){};
		};
		/**
		 * @class ElectroPhysiologyDiffusionRelaxationComplex
		 * @brief Compute the diffusion relaxation process
		 */
		class ElectroPhysiologyDiffusionRelaxationComplex
			: public RelaxationOfAllDiffusionSpeciesRK2<
				  RelaxationOfAllDiffussionSpeciesComplex<RealBody, SolidParticles, Solid, RealBody, SolidParticles, Solid>>
		{
		public:
			explicit ElectroPhysiologyDiffusionRelaxationComplex(ComplexBodyRelation &complex_relation)
				: RelaxationOfAllDiffusionSpeciesRK2(complex_relation){};
			virtual ~ElectroPhysiologyDiffusionRelaxationComplex(){};
		};
		/**
		 * @class ElectroPhysiologyReactionRelaxationForward
		 * @brief Solve the reaction ODE equation of trans-membrane potential
		 * using forward sweeping
		 */
		class ElectroPhysiologyReactionRelaxationForward
			: public RelaxationOfAllReactionsForward<RealBody, SolidParticles, Solid>
		{
		public:
			explicit ElectroPhysiologyReactionRelaxationForward(RealBody &real_body)
				: RelaxationOfAllReactionsForward<RealBody, SolidParticles, Solid>(real_body){};
			virtual ~ElectroPhysiologyReactionRelaxationForward(){};
		};
		/**
		 * @class ElectroPhysiologyReactionRelaxationForward
		 * @brief Solve the reaction ODE equation of trans-membrane potential
		 * using backward sweeping
		 */
		class ElectroPhysiologyReactionRelaxationBackward
			: public RelaxationOfAllReactionsBackward<RealBody, SolidParticles, Solid>
		{
		public:
			explicit ElectroPhysiologyReactionRelaxationBackward(RealBody &real_body)
				: RelaxationOfAllReactionsBackward<RealBody, SolidParticles, Solid>(real_body){};
			virtual ~ElectroPhysiologyReactionRelaxationBackward(){};
		};
		/**
		 * @class ApplyStimulusCurrents
		 * @brief Apply specific stimulus currents
		 * This is a abstract class to be override for case specific implementations.
		 */
		class ApplyStimulusCurrents : public ParticleDynamicsSimple,
									  public ElectroPhysiologyDataDelegateSimple
		{
		public:
			explicit ApplyStimulusCurrents(RealBody &real_body)
				: ParticleDynamicsSimple(real_body),
				  ElectroPhysiologyDataDelegateSimple(real_body) {}
			virtual ~ApplyStimulusCurrents(){};
		};
	}
}
#endif // ELECTRO_PHYSIOLOGY_H