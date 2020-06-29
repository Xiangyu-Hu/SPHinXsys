/**
 * @file 	obeserver_dynamics.h
 * @brief 	There are the classes for obsevser bodies to record the state of the flow or
 *			solid in given locations. Mostly, this is done by an interpolation alogortim.   
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 * @note 	c++ knowledge : 
 * 			https://stackoverflow.com/questions/1120833/derived-template-class-access-to-base-class-member-data
 * 			Chi ZHANG
 */

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	namespace observer_dynamics
	{
		template <class ObserverParticlesType, class TargetParticlesType>
		using ContactInterpolation = ParticleDynamicsContact<SPHBody, ObserverParticlesType, 
			BaseMaterial, SPHBody, TargetParticlesType>;

		template <class ObserverParticlesType, class TargetParticlesType>
		using ComplexInterpolation = ParticleDynamicsComplex<SPHBody, ObserverParticlesType, 
			BaseMaterial, SPHBody, TargetParticlesType>;

		template <class BaseParticlesType, class TargetParticlesType>
		using ObservingDyanmics = ParticleDynamicsContact<SPHBody, BaseParticlesType, BaseMaterial,SPHBody, TargetParticlesType, BaseMaterial>;

		/**
		 * @class ObservingAQuantityFromABody
		 * @brief Observering a given member data in the particles of a general body
		 */
		template <class DataType, class TargetParticlesType, class TargetDataType,
			StdLargeVec<TargetDataType> TargetParticlesType:: * TrgtDataMemPtr, DataType TargetDataType:: * TrgtMemPtr>
		class ObservingAQuantityFromABody : public ObservingDyanmics<BaseParticles, TargetParticlesType>
		{
		protected:
			/** Observed quantities saved here. */
			StdLargeVec<DataType>  observed_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override
			{
				TargetParticlesType* target_particles = this->interacting_particles_[interacting_body_index];
				StdLargeVec<TargetDataType>& target_data = target_particles->*TrgtDataMemPtr;

				DataType observed_quantity(0);
				Real ttl_weight(0);
				Neighborhood& contact_neighborhood 
					= (*this->current_interacting_configuration_[interacting_body_index])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*this->interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
					Real Vol_j = base_particle_data_j.Vol_;
					TargetDataType& target_data_j = target_data[index_particle_j];

					Real weight_j = neighboring_particle->W_ij_ * Vol_j;
					observed_quantity += weight_j * target_data_j.* TrgtMemPtr;
					ttl_weight += weight_j;
				}
				observed_quantities_[index_particle_i] = observed_quantity / ttl_weight;
			};
		public:
			explicit ObservingAQuantityFromABody(SPHBody* observer, SPHBody* target)
				: ObservingDyanmics<BaseParticles, TargetParticlesType>(observer, { target }) {
				for (size_t i = 0; i < observer->number_of_particles_; ++i) observed_quantities_.push_back(DataType(0));
			};
			virtual ~ObservingAQuantityFromABody() {};
		};

		/**
		 * @class ObservingADiffusionReactionQuantityFromABody
		 * @brief Observing a diffusion-reaction quantity from a body
		 */
		template <class DiffusionReactionParticlesType>
		class ObservingADiffusionReactionQuantityFromABody 
				: public ObservingDyanmics<BaseParticles, DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** Observed quantities saved here. */
			StdLargeVec<Real>  observed_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) override
			{
				DiffusionReactionParticlesType* target_particles = this->interacting_particles_[interacting_body_index];
				StdLargeVec<DiffusionReactionData>& target_data = target_particles->diffusion_reaction_data_;

				Real observed_quantity(0);
				Real ttl_weight(0);
				Neighborhood& contact_neighborhood 
					= (*this->current_interacting_configuration_[interacting_body_index])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*this->interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
					Real Vol_j = base_particle_data_j.Vol_;

					Real weight_j = neighboring_particle->W_ij_ * Vol_j;
					observed_quantity += weight_j * target_data[index_particle_j].species_n_[species_index_];
					ttl_weight += weight_j;
				}
				observed_quantities_[index_particle_i] = observed_quantity / ttl_weight;
			};
		public:
			explicit ObservingADiffusionReactionQuantityFromABody(string species_name, SPHBody* observer, SPHBody* target)
				: ObservingDyanmics<BaseParticles, DiffusionReactionParticlesType>(observer, { target }) {
				species_indexes_map_ = this->interacting_particles_[0]->getSpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];
				for (size_t i = 0; i < observer->number_of_particles_; ++i) observed_quantities_.push_back(0.0);
			};
			virtual ~ObservingADiffusionReactionQuantityFromABody() {};
		};

		/**
		 * @class InterpolatingAQuantity
		 * @brief interpolate a given member data in the particles of a general body from another body
		 */
		template <class DataType, class ObserverParticlesType, class ObserverDataType, class TargetParticlesType, class TargetDataType,
			StdLargeVec<ObserverDataType> ObserverParticlesType:: * ObrsvrDataMemPtr, StdLargeVec<TargetDataType> TargetParticlesType:: * TrgtDataMemPtr,
			DataType ObserverDataType:: * ObrsvrMemPtr, DataType TargetDataType:: * TrgtMemPtr>
			class InterpolatingAQuantity : public ComplexInterpolation<ObserverParticlesType, TargetParticlesType>
		{
		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				ObserverParticlesType* particles = this->particles_;
				StdLargeVec<ObserverDataType>& observer_data = this->particles_->*ObrsvrDataMemPtr;
				ObserverDataType& observer_data_i = observer_data[index_particle_i];

				DataType observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->current_interacting_configuration_.size(); ++k)
				{
					TargetParticlesType* target_particles = this->interacting_particles_[k];
					StdLargeVec<BaseParticleData>& target_base_particle_data = target_particles->base_particle_data_;
					StdLargeVec<TargetDataType>& target_data = target_particles->*TrgtDataMemPtr;

					Neighborhood& contact_neighborhood = (*this->current_interacting_configuration_[k])[index_particle_i];
					NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
					for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
					{
						BaseNeighborRelation* neighboring_particle = contact_neighors[n];
						size_t index_particle_j = neighboring_particle->j_;
						BaseParticleData& base_particle_data_j = target_base_particle_data[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;
						TargetDataType& target_data_j = target_data[index_particle_j];

						Real weight_j = neighboring_particle->W_ij_ * Vol_j;
						observed_quantity += weight_j * target_data_j.*TrgtMemPtr;
						ttl_weight += weight_j;
					}
				}
				observer_data_i.* ObrsvrMemPtr = observed_quantity / ttl_weight;
			};
		public:
			explicit InterpolatingAQuantity(SPHBody* observer, StdVec<SPHBody*> target_bodies)
				: ComplexInterpolation<ObserverParticlesType, TargetParticlesType>(observer, target_bodies) {};
			virtual ~InterpolatingAQuantity() {};
		};

		/**
		 * @class InterpolatingADiffusionReactionQuantity
		 * @brief Observering general body
		 */
		template <class ObserverParticlesType, class ObserverDataType, class DiffusionReactionParticlesType,
			StdLargeVec<ObserverDataType> ObserverParticlesType:: * ObrsvrDataMemPtr, Real ObserverDataType:: * ObrsvrMemPtr>
			class InterpolatingADiffusionReactionQuantity
			: public ComplexInterpolation<ObserverParticlesType, DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** Observed quantities saved here. */
			StdLargeVec<Real>  observed_quantities_;

			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				ObserverParticlesType* particles = this->particles_;
				StdLargeVec<ObserverDataType>& observer_data = this->particles_->*ObrsvrDataMemPtr;
				ObserverDataType& observer_data_i = observer_data[index_particle_i];

				Real observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->current_interacting_configuration_.size(); ++k)
				{
					DiffusionReactionParticlesType* target_particles = this->interacting_particles_[k];
					StdLargeVec<BaseParticleData>& target_base_particle_data = target_particles->base_particle_data_;
					StdLargeVec<DiffusionReactionData>& target_diffusion_reaction_data = target_particles->diffusion_reaction_data_;

					Neighborhood& contact_neighborhood = (*this->current_interacting_configuration_[k])[index_particle_i];
					NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
					for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
					{
						BaseNeighborRelation* neighboring_particle = contact_neighors[n];
						size_t index_particle_j = neighboring_particle->j_;
						BaseParticleData& base_particle_data_j = target_base_particle_data[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;
						DiffusionReactionData& target_data_j = target_diffusion_reaction_data[index_particle_j];

						Real weight_j = neighboring_particle->W_ij_ * Vol_j;
						observed_quantity += weight_j * target_data_j.species_n_[species_index_];
						ttl_weight += weight_j;
					}
				}
				observer_data_i.*ObrsvrMemPtr = observed_quantity / ttl_weight;
			};
		public:
			explicit InterpolatingADiffusionReactionQuantity(string species_name, SPHBody* observer, StdVec<SPHBody*> target_bodies)
				: ComplexInterpolation<ObserverParticlesType, DiffusionReactionParticlesType>(observer, target_bodies) 
			{
				species_indexes_map_ = this->interacting_particles_[0]->getSpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];
			};
			virtual ~InterpolatingADiffusionReactionQuantity() {};
		};
		
		/**
		* @class CorrectKenelWeightsForInterpolation
		* @brief  correct kenel weights for interpolation
		*/
		class CorrectKenelWeightsForInterpolation : 
			public ParticleDynamicsComplex<SPHBody, BaseParticles, BaseMaterial, SPHBody, BaseParticles, BaseMaterial>
		{
		protected:
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CorrectKenelWeightsForInterpolation(SPHBody *body, StdVec<SPHBody*> interacting_bodies)
				: ParticleDynamicsComplex<SPHBody, BaseParticles, BaseMaterial, SPHBody, BaseParticles, BaseMaterial>(body, interacting_bodies) {};
			virtual ~CorrectKenelWeightsForInterpolation() {};
		};
	}
}
