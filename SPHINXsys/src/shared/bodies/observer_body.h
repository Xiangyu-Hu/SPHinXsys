/**
 * @file 	observer_body.h
 * @brief 	This is the class for bodies recording the state of the flow or solid
 * 			in given locations by observer particles at these locations. 
 *			This body has no inner configuration so that no
 * 			cell linked list is required. However, it has contact configuration to
 * 			the body it is obeserving.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_body.h"

using namespace std;

namespace SPH {
	/**
	 * @brief preclaimed class. 
	 */
	class SPHSystem;
	class ObserverParticles;
	/**
	 * @class ObserverBody
	 * @brief Declaration of observerbody which is derived from FictitiousBody
	 */
	class ObserverBody : public FictitiousBody
	{
	public:
		ObserverBody(SPHSystem &system, string body_name,
			ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op);
		virtual ~ObserverBody() {};
		
		/** particles in this body */
		ObserverParticles &observer_particles_;	
	};

	class ObserverLagrangianBody : public ObserverBody
	{
	public:
		ObserverLagrangianBody(SPHSystem &system, string body_name,
			ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op);
		virtual ~ObserverLagrangianBody() {};

		/** Build contact configuration. */
		virtual void BuildContactConfiguration() override;
	};

	class ObserverEulerianBody : public ObserverBody
	{
	public:
		ObserverEulerianBody(SPHSystem &system, string body_name,
			ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op);
		virtual ~ObserverEulerianBody() {};

		/** Build contact configuration. */
		virtual void BuildContactConfiguration() override;
	};
}