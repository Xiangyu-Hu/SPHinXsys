/**
 * @file 	obeserver_dynamics.h
 * @brief 	There are the classes for obsevser bodies to record the state of the flow or
 *			solid in given locations. Mostly, this is done by an interpolation alogortim.   
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	namespace observer_dynamics
	{
		template <class BodyType, class ParticlesType, class MaterialType, class InteractingBodyType, class InteractingParticlesType>
		using ContactInterpolation = ParticleDynamicsContact<BodyType, ParticlesType, MaterialType, 
			InteractingBodyType, InteractingParticlesType>;
		
		/**
		 * @class InterpolationFromABody
		 * @brief Interpolate any variable from general body
		 */	
		template <class BodyType, class ParticlesType, class MaterialType, class InteractingBodyType, class InteractingParticlesType>
		class InterpolationFromABody : public ContactInterpolation<BodyType, ParticlesType, MaterialType, InteractingBodyType, InteractingParticlesType>
		{
		protected:
			/** Abstract method for observing. */
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		public:
			explicit InterpolationFromABody(BodyType *body, InteractingBodyType *interacting_body)
				: ContactInterpolation<BodyType, ParticlesType, MaterialType, InteractingBodyType, InteractingParticlesType>(body, { interacting_body }) {};
			virtual ~InterpolationFromABody() {};
		};

		/**
		 * @class ObserveABody
		 * @brief Observering general body
		 */	
		template <class InteractingBodyType, class InteractingParticlesType>
		class ObserveABody : public ContactInterpolation<ObserverBody, ObserverParticles, BaseMaterial, InteractingBodyType, InteractingParticlesType>
		{
		protected:
			/** Abstract method for observing. */
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		public:
			explicit ObserveABody(ObserverBody *body, InteractingBodyType *interacting_body)
				: ContactInterpolation<ObserverBody, ObserverParticles, BaseMaterial, InteractingBodyType, InteractingParticlesType>(body, { interacting_body }) {};
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
			/** Observed quantities saved here. */
			StdVec<FluidQuantityType>  fluid_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			/** Abstarct method to define the quantity to be observed. */
			virtual FluidQuantityType GetAFluidQuantity(size_t index_particle_j, FluidParticles &particles) = 0;
		public:
			explicit ObserveAFluidQuantity(ObserverBody *body, FluidBody *interacting_body)
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
			/** Define to observe fluid pressure. */
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
			/** Define to observe fluid velocity. */
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
		 * @brief Observe an elastic solid quantity.
		 * This class is the couterpart to the class
		 * ObserveAFluidQuantity
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
				for (size_t i = 0; i < body->number_of_particles_; ++i) 
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
			/** Define to observe the solid dispalacement. */
			virtual Vecd GetAnElasticSolidQuantity(size_t index_particle_j, ElasticSolidParticles &particles) override {
				return particles.base_particle_data_[index_particle_j].pos_n_;
			};

		public:
			ObserveElasticDisplacement(ObserverBody *body, SolidBody *interacting_body)
				: ObserveAnElasticSolidQuantity(body, interacting_body) {};
			virtual ~ObserveElasticDisplacement() {};
		};
		 /**
		 * @class ObserveAnElasticSolidQuantity
		 * @brief Observe an muscle quantity.
		 * This class is the couterpart to the class
		 * ObserveAFluidQuantity and ObserveAnElasticSolidQuantity
		 */
		template <typename ElectroPhysiologyQuantityType>
		class ObserveAElectroPhysiologyQuantity : public ObserveABody<SolidBody, ElectroPhysiologyParticles>
		{

		protected:
			ElectroPhysiologyParticles* electro_physiology_particles_;
			map<string, size_t> species_indexes_map_;

			StdVec<ElectroPhysiologyQuantityType> electro_physiology_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override;
			virtual ElectroPhysiologyQuantityType GetAMuscleQuantity(size_t index_particle_j, ElectroPhysiologyParticles&particles) = 0;

		public:
			explicit ObserveAElectroPhysiologyQuantity(ObserverBody *body, SolidBody *interacting_body)
				: ObserveABody(body, interacting_body) 
				{
				ElectroPhysiologyParticles* electro_physiology_particles_
					= interacting_particles_[0];
				species_indexes_map_ = electro_physiology_particles_->getSpeciesIndexMap();
					for (size_t i = 0; i < body->number_of_particles_; ++i) 
						electro_physiology_quantities_.push_back(ElectroPhysiologyQuantityType(0));
				};
			virtual ~ObserveAElectroPhysiologyQuantity() {};
		};
		/**
		 * @class ObserveElectroPhysiologyVoltage
		 * @brief observe elastic displacement
		 */
		class ObserveElectroPhysiologyVoltage : public ObserveAElectroPhysiologyQuantity<Real>
		{

		protected:
			/** Index of voltage. */
			size_t voltage_;

			/** Define to observe the solid dispalacement. */
			virtual Real GetAMuscleQuantity(size_t index_particle_j, ElectroPhysiologyParticles& particles) override
			{
				return particles.diffusion_reaction_data_[index_particle_j].species_n_[voltage_];
			};

		public:
			ObserveElectroPhysiologyVoltage(ObserverBody *body, SolidBody *interacting_body)
				: ObserveAElectroPhysiologyQuantity(body, interacting_body) 
			{
				voltage_ = species_indexes_map_["Voltage"];
			};
			virtual ~ObserveElectroPhysiologyVoltage() {};
		};

		/**
		 * @class Interpolate
		 * @brief Observe an muscle quantity.
		 * This class is the couterpart to the class
		 * ObserveAFluidQuantity and ObserveAnElasticSolidQuantity
		 */
		typedef InterpolationFromABody<SolidBody, ActiveMuscleParticles, ActiveMuscle, SolidBody, ElectroPhysiologyParticles> ElectroPhysiologyInterpolation;
		
		template <typename ElectroPhysiologyQuantityType>
		class ElectroPhysiologyQuantityInterpolation : public ElectroPhysiologyInterpolation
		{

		protected:
			ElectroPhysiologyParticles* electro_physiology_particles_;
			map<string, size_t> species_indexes_map_;

			virtual ElectroPhysiologyQuantityType GetAMuscleQuantity(size_t index_particle_j, ElectroPhysiologyParticles&particles) = 0;
			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override ;

		public:
			explicit ElectroPhysiologyQuantityInterpolation(SolidBody *body, SolidBody *interacting_body)
				: ElectroPhysiologyInterpolation(body, interacting_body) 
				{
					ElectroPhysiologyParticles* electro_physiology_particles_ = interacting_particles_[0];
					species_indexes_map_ = electro_physiology_particles_->getSpeciesIndexMap();
				};
			virtual ~ElectroPhysiologyQuantityInterpolation() {};
		};

		/**
		 * @class ObserveMuscleVoltage
		 * @brief observe elastic displacement
		 */
		class ActiveContractStressInterpolation : public ElectroPhysiologyQuantityInterpolation<Real>
		{

		protected:
			/** Index of active constract stress. */
			size_t active_contract_stress_;

			/** Define to observe the solid dispalacement. */
			virtual Real GetAMuscleQuantity(size_t index_particle_j, ElectroPhysiologyParticles& particles) override
			{
				return particles.diffusion_reaction_data_[index_particle_j].species_n_[active_contract_stress_];
			};

		public:
			ActiveContractStressInterpolation(SolidBody *body, SolidBody *interacting_body)
				: ElectroPhysiologyQuantityInterpolation(body, interacting_body) 
			{
				active_contract_stress_ = species_indexes_map_["ActiveContractionStress"];
			};
			virtual ~ActiveContractStressInterpolation() {};
		};
	}
}