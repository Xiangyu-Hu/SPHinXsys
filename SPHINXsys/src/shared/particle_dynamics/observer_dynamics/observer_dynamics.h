/**
 * @file 	obeserver_dynamics.h
 * @brief 	This is the class for obsevser bodies to record the state of the flow or
 *			solid in given locations. Mostly, this is done by interpolation alogortim.   
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	namespace observer_dynamics
	{
		/**
		 * @class ObserveABody
		 * @brief observer physical field from fluid
		 */	
		template <class InteractingBodytype, class InteractingParticlesType>
		class ObserveABody 
			: public ParticleDynamicsContact<ObserverBody, ObserverParticles, InteractingBodytype, InteractingParticlesType>
		{
		
		protected:
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;

		public:
			ObserveABody(ObserverBody *body, InteractingBodytype *interacting_body)
				: ParticleDynamicsContact(body, { interacting_body }) {};
			virtual ~ObserveABody() {};
		};
		
		/**
		 * @class ObserveFluidPressure
		 * @brief observe fluid pressure
		 */
		class ObserveFluidPressure : public ObserveABody<WeaklyCompressibleFluidBody, WeaklyCompressibleFluidParticles>
		{

		protected:

			StdVec<Real>  pressures_;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			ObserveFluidPressure(ObserverBody *body, WeaklyCompressibleFluidBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_particles_; ++i) pressures_.push_back(Real(0));
			};
			virtual ~ObserveFluidPressure() {};
		};

		/**
		 * @class ObserveElasticDisplacement
		 * @brief observe elastic displacement
		 */
		class ObserveElasticDisplacement : public ObserveABody<ElasticBody, ElasticBodyParticles>
		{

		protected:
			StdVec<Vecd>  displacements_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;

		public:
			ObserveElasticDisplacement(ObserverBody *body, ElasticBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_particles_; ++i) displacements_.push_back(Vecd(0));
			};
			virtual ~ObserveElasticDisplacement() {};
		};
	}
}