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
		template <class InteractingBodyType, class InteractingParticlesType>
		using ObserverDynamicsContact = ParticleDynamicsContact<ObserverBody, ObserverParticles, Material,
			InteractingBodyType, InteractingParticlesType>;
		
		/**
		 * @class ObserveABody
		 * @brief observer physical field from fluid
		 */	
		template <class InteractingBodyType, class InteractingParticlesType>
		class ObserveABody : public ObserverDynamicsContact<InteractingBodyType, InteractingParticlesType>
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		public:
			explicit ObserveABody(ObserverBody *body, InteractingBodyType *interacting_body)
				: ObserverDynamicsContact<InteractingBodyType, InteractingParticlesType>(body, { interacting_body }) {};
			virtual ~ObserveABody() {};
		};
		
		/**
		 * @class ObserveAFluidQuantity
		 * @brief observe a fluid quantity
		 */
		template <typename FluidQuantityType>
		class ObserveAFluidQuantity : public ObserveABody<FluidBody, FluidParticles>
		{
		protected:
			StdVec<FluidQuantityType>  fluid_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual FluidQuantityType GetAFluidQuantity(size_t index_particle_j, FluidParticles &particles) = 0;
		public:
			explicit ObserveAFluidQuantity(ObserverBody *body, FluidBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_real_particles_; ++i) fluid_quantities_.push_back(FluidQuantityType(0));
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
			virtual Real GetAFluidQuantity(size_t index_particle_j, FluidParticles &particles) override { 
				return particles.fluid_particle_data_[index_particle_j].p_; 
			};
		public:
			explicit ObserveFluidPressure(ObserverBody *body, FluidBody *interacting_body)
				: ObserveAFluidQuantity(body, interacting_body) {};
			virtual ~ObserveFluidPressure() {};
		};

		/**
		 * @class ObserveFluidVelocity
		 * @brief observe fluid velocity
		 */
		class ObserveFluidVelocity : public ObserveAFluidQuantity<Vecd>
		{
		protected:
			virtual Vecd GetAFluidQuantity(size_t index_particle_j, FluidParticles &particles) override {
				return particles.base_particle_data_[index_particle_j].vel_n_;
			};
		public:
			explicit ObserveFluidVelocity(ObserverBody *body, FluidBody *interacting_body)
				: ObserveAFluidQuantity(body, interacting_body) {};
			virtual ~ObserveFluidVelocity() {};
		};

		/**
		 * @class ObserveAnElasticSolidQuantity
		 * @brief observe an elastic solid quantity
		 */
		template <typename ElasticSolidQuantityType>
		class ObserveAnElasticSolidQuantity : public ObserveABody<SolidBody, ElasticSolidParticles>
		{

		protected:
			StdVec<ElasticSolidQuantityType>  elastic_body_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual ElasticSolidQuantityType GetAnElasticSolidQuantity(size_t index_particle_j, ElasticSolidParticles &particles) = 0;

		public:
			explicit ObserveAnElasticSolidQuantity(ObserverBody *body, SolidBody *interacting_body)
				: ObserveABody(body, interacting_body) {
				for (size_t i = 0; i < body->number_of_real_particles_; ++i) 
					elastic_body_quantities_.push_back(ElasticSolidQuantityType(0));
			};
			virtual ~ObserveAnElasticSolidQuantity() {};
		};

		/**
		 * @class ObserveElasticDisplacement
		 * @brief observe elastic displacement
		 */
		class ObserveElasticDisplacement : public ObserveAnElasticSolidQuantity<Vecd>
		{

		protected:
			virtual Vecd GetAnElasticSolidQuantity(size_t index_particle_j, ElasticSolidParticles &particles) override {
				return particles.base_particle_data_[index_particle_j].pos_n_;
			};

		public:
			ObserveElasticDisplacement(ObserverBody *body, SolidBody *interacting_body)
				: ObserveAnElasticSolidQuantity(body, interacting_body) {};
			virtual ~ObserveElasticDisplacement() {};
		};
	}
}