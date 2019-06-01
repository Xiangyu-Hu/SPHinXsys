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
		 * @class ObserveAFluidQuantity
		 * @brief observe a fluid quantity
		 */
		template <typename FluidQuantityType>
		class ObserveAFluidQuantity : public ObserveABody<WeaklyCompressibleFluidBody, WeaklyCompressibleFluidParticles>
		{

		protected:
			StdVec<FluidQuantityType>  fluid_quantities_;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual FluidQuantityType GetAFluidQuantity(size_t index_particle_j, WeaklyCompressibleFluidParticles &particles) = 0;

		public:
			ObserveAFluidQuantity(ObserverBody *body, WeaklyCompressibleFluidBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_particles_; ++i) fluid_quantities_.push_back(FluidQuantityType(0));
			};
			virtual ~ObserveAFluidQuantity() {};
		};
		
		/**
		 * @class ObserveFluidPressure
		 * @brief observe fluid pressure
		 */
		class ObserveFluidPressure : public ObserveAFluidQuantity<Real>
		{
		protected:
			virtual Real GetAFluidQuantity(size_t index_particle_j, WeaklyCompressibleFluidParticles &particles) override
			{ return particles.fluid_data_[index_particle_j].p_; };
		public:
			ObserveFluidPressure(ObserverBody *body, WeaklyCompressibleFluidBody *interacting_body)
				: ObserveAFluidQuantity(body, interacting_body) {};
			virtual ~ObserveFluidPressure() {};
		};


		/**
		 * @class ObserveFluidPressure
		 * @brief observe fluid pressure
		 */
		class ObserveFluidVelocity : public ObserveAFluidQuantity<Vecd>
		{
		protected:
			virtual Vecd GetAFluidQuantity(size_t index_particle_j, WeaklyCompressibleFluidParticles &particles) override
			{
				return particles.base_particle_data_[index_particle_j].vel_n_;
			};
		public:
			ObserveFluidVelocity(ObserverBody *body, WeaklyCompressibleFluidBody *interacting_body)
				: ObserveAFluidQuantity(body, interacting_body) {};
			virtual ~ObserveFluidVelocity() {};
		};

		/**
		 * @class ObserveElasticDisplacement
		 * @brief observe elastic displacement
		 */
		template <typename ElasticBodyQuantityType>
		class ObserveAnElasticBodyQuantity : public ObserveABody<ElasticBody, ElasticBodyParticles>
		{

		protected:
			StdVec<ElasticBodyQuantityType>  elastic_body_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual ElasticBodyQuantityType GetAnElasticBodyQuantity(size_t index_particle_j, ElasticBodyParticles &particles) = 0;

		public:
			ObserveAnElasticBodyQuantity(ObserverBody *body, ElasticBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_particles_; ++i) elastic_body_quantities_.push_back(ElasticBodyQuantityType(0));
			};
			virtual ~ObserveAnElasticBodyQuantity() {};
		};

		/**
		 * @class ObserveElasticDisplacement
		 * @brief observe elastic displacement
		 */
		class ObserveElasticDisplacement : public ObserveAnElasticBodyQuantity<Vecd>
		{

		protected:
			virtual Vecd GetAnElasticBodyQuantity(size_t index_particle_j, ElasticBodyParticles &particles) override
			{
				return particles.base_particle_data_[index_particle_j].pos_n_;
			};

		public:
			ObserveElasticDisplacement(ObserverBody *body, ElasticBody *interacting_body)
				: ObserveAnElasticBodyQuantity(body, interacting_body) {};
			virtual ~ObserveElasticDisplacement() {};
		};
	}
}