/**
 * @file 	eletro_physiology.h
 * @brief 	In is file, we declear the diffusion dynamics and the corresponing function. 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
#pragma once

#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "base_kernel.h"

namespace SPH
{
	class ElectrophysiologyReaction;
	
	namespace electro_physiology
	{
        typedef ParticleDynamicsSimple<SolidBody, MuscleParticles, Muscle> ElectroPhysiologySimple;
        typedef ParticleDynamicsReduce<Real, ReduceMin, SolidBody, MuscleParticles, Muscle> ElectroPhysiologyMinimum;
        typedef ParticleDynamicsInner1Level<SolidBody, MuscleParticles, Muscle> ElectroPhysiologyInnver1Level;
		/**
		 * @class  OffsetInitialParticlePosition
		 * @brief  set initial condition for a cardiac muscle body
		*/
		class OffsetInitialParticlePosition : public ElectroPhysiologySimple
		{
			Vecd offset_;
		protected:
			//default for set all particle at rest
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			OffsetInitialParticlePosition(SolidBody *body, Vecd offset)
				: ElectroPhysiologySimple(body), offset_(offset){};
			virtual ~OffsetInitialParticlePosition() {};
		};
        /**
		 * @class ElectroPhysiologyInitialCondition
		 * @brief  set initial condition for a muscle body
		*/
		class ElectroPhysiologyInitialCondition : public ElectroPhysiologySimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ElectroPhysiologyInitialCondition(SolidBody *body)
				: ElectroPhysiologySimple(body) {};
			virtual ~ElectroPhysiologyInitialCondition() {};
		};
        /**
		* @class getDiffusionTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class getDiffusionTimeStepSize : public ElectroPhysiologyMinimum
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit getDiffusionTimeStepSize(SolidBody* body);
			virtual ~getDiffusionTimeStepSize() {};
		};
        /**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration : public ElectroPhysiologyInnver1Level
		{
		protected:
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CorrectConfiguration(SolidBody *body) : ElectroPhysiologyInnver1Level(body) {};
			virtual ~CorrectConfiguration() {};
		};
        /**
		* @class DiffusionRelaxation
		* @brief Compute the diffusion relaxation process
		*/
		class DiffusionRelaxation : public ElectroPhysiologyInnver1Level
		{
		protected:
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			DiffusionRelaxation(SolidBody *body): ElectroPhysiologyInnver1Level(body) {};
			virtual ~DiffusionRelaxation() {};
		};
        /**
		 * @class TransmembranePotentialReaction
		 * @brief Solv the reaction ODE equation of trans-membrane potential
		*/
		class TransmembranePotentialReaction : public ElectroPhysiologySimple
		{
		protected:
            ElectrophysiologyReaction* reaction_model_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			TransmembranePotentialReaction(SolidBody *body);
			virtual ~TransmembranePotentialReaction() {};
		};
        /**
		 * @class GateVariableReaction
		 * @brief Solvin the reaction ODE equation of gate variable.
		*/
		class GateVariableReaction : public ElectroPhysiologySimple
		{
		protected:
		 	ElectrophysiologyReaction* reaction_model_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			GateVariableReaction(SolidBody *body);
			virtual ~GateVariableReaction() {};
		};    
        /**
		 * @class ApplyStimulusCurrents
		 * @brief Apply specific stimulus currents
		*/
		class ApplyStimulusCurrents : public ElectroPhysiologySimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ApplyStimulusCurrents(SolidBody *body);
			virtual ~ApplyStimulusCurrents() {};
		};
    }
}
