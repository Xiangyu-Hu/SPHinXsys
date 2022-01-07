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
 * @file 	eletro_physiology.h
 * @brief 	In is file, we declaim the dynamics relevant to electrophysiology,
 * including diffusion, reaction and muscle activation. 
 * @author 	Chi Zhang and Xiangyu Hu
 */

#ifndef ELECTRO_PHYSIOLOGY_H
#define ELECTRO_PHYSIOLOGY_H



#include "all_particle_dynamics.h"
#include "diffusion_reaction_particles.h"
#include "particle_dynamics_diffusion_reaction.h"
#include "diffusion_reaction.h"
#include "elastic_solid.h"
#include "base_kernel.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
	namespace electro_physiology
	{
		typedef DiffusionReactionSimpleData<SolidBody, SolidParticles, Solid> ElectroPhysiologyDataDelegateSimple;
		typedef DiffusionReactionInnerData<SolidBody, SolidParticles, Solid> ElectroPhysiologyDataDelegateInner;
		/**
		 * @class ElectroPhysiologyInitialCondition
		 * @brief  set initial condition for a muscle body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElectroPhysiologyInitialCondition : 
			public ParticleDynamicsSimple,
			public ElectroPhysiologyDataDelegateSimple
		{
		public:
			explicit ElectroPhysiologyInitialCondition(SolidBody &solid_body);
			virtual ~ElectroPhysiologyInitialCondition() {};
		protected:
			StdLargeVec<Vecd>& pos_n_;
			StdVec<StdLargeVec<Real>>& species_n_;
		};
        /**
		* @class GetElectroPhysiologyTimeStepSize
		* @brief Computing the time step size from diffusion criteria
		*/
		class GetElectroPhysiologyTimeStepSize 
			: public GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid>
		{
		public:
			explicit GetElectroPhysiologyTimeStepSize(SolidBody &solid_body)
				: GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid>(solid_body) {};
			virtual ~GetElectroPhysiologyTimeStepSize() {};
		};
        /**
		* @class ElectroPhysiologyDiffusionRelaxationInner
		* @brief Compute the diffusion relaxation process
		*/
		class ElectroPhysiologyDiffusionRelaxationInner : 
			public RelaxationOfAllDiffusionSpeciesRK2<SolidBody, SolidParticles, Solid,
			RelaxationOfAllDiffussionSpeciesInner<SolidBody, SolidParticles, Solid>, 
			BaseBodyRelationInner>
		{
		public:
			explicit ElectroPhysiologyDiffusionRelaxationInner(BaseBodyRelationInner &inner_relation)
				: RelaxationOfAllDiffusionSpeciesRK2(inner_relation) {};
			virtual ~ElectroPhysiologyDiffusionRelaxationInner() {};
		};
		/**
		* @class ElectroPhysiologyDiffusionRelaxationComplex
		* @brief Compute the diffusion relaxation process
		*/
		class ElectroPhysiologyDiffusionRelaxationComplex : 
			public RelaxationOfAllDiffusionSpeciesRK2<SolidBody, SolidParticles, Solid,
			RelaxationOfAllDiffussionSpeciesComplex<SolidBody, SolidParticles, Solid, SolidBody, SolidParticles, Solid>, 
			ComplexBodyRelation>
		{
		public:
			explicit ElectroPhysiologyDiffusionRelaxationComplex(ComplexBodyRelation &complex_relation)
				: RelaxationOfAllDiffusionSpeciesRK2(complex_relation) {};
			virtual ~ElectroPhysiologyDiffusionRelaxationComplex() {};
		};
        /**
		 * @class ElectroPhysiologyReactionRelaxationForward
		 * @brief Solve the reaction ODE equation of trans-membrane potential
		 * using forward sweeping
		*/
		class ElectroPhysiologyReactionRelaxationForward
			: public RelaxationOfAllReactionsForward<SolidBody, SolidParticles, Solid>
		{
		public:
			explicit ElectroPhysiologyReactionRelaxationForward(SolidBody &solid_body)
				: RelaxationOfAllReactionsForward<SolidBody, SolidParticles, Solid>(solid_body) {};
			virtual ~ElectroPhysiologyReactionRelaxationForward() {};
		};
		/**
		 * @class ElectroPhysiologyReactionRelaxationForward
		 * @brief Solve the reaction ODE equation of trans-membrane potential
		 * using backward sweeping
		*/
		class ElectroPhysiologyReactionRelaxationBackward
			: public RelaxationOfAllReactionsBackward<SolidBody, SolidParticles, Solid>
		{
		public:
			explicit ElectroPhysiologyReactionRelaxationBackward(SolidBody &solid_body)
				: RelaxationOfAllReactionsBackward<SolidBody, SolidParticles, Solid>(solid_body) {};
			virtual ~ElectroPhysiologyReactionRelaxationBackward() {};
		};
		/**
		 * @class ApplyStimulusCurrents
		 * @brief Apply specific stimulus currents
		 * This is a abstract class to be override for case specific implementations.
		*/
		class ApplyStimulusCurrents :
			public ParticleDynamicsSimple,
			public ElectroPhysiologyDataDelegateSimple
		{
		public:
			explicit ApplyStimulusCurrents(SolidBody &solid_body) : 
				ParticleDynamicsSimple(solid_body), 
				ElectroPhysiologyDataDelegateSimple(solid_body) {}
			virtual ~ApplyStimulusCurrents() {};
		};
    }
}
#endif //ELECTRO_PHYSIOLOGY_H