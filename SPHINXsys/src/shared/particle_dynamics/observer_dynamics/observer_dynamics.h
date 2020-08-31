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
 * @file 	observer_dynamics.h
 * @brief 	There are the classes for observser bodies to record the state of the flow or
 *			solid in given locations. Mostly, this is done by an interpolation algorithm.   
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
		using InterpolationDynamics 
			= ParticleDynamicsContact<SPHBody, ObserverParticlesType, BaseMaterial, SPHBody, TargetParticlesType>;

		template <class TargetParticlesType>
		using ObservationDynamics 
			= ParticleDynamicsContact<SPHBody, BaseParticles, BaseMaterial,SPHBody, TargetParticlesType, BaseMaterial>;

		/**
		 * @class ObservingAQuantity
		 * @brief Observering a given member data in the particles of a general body
		 */
		template <class DataType, class TargetParticlesType, class TargetDataType,
			StdLargeVec<TargetDataType> TargetParticlesType:: * TrgtDataMemPtr, DataType TargetDataType:: * TrgtMemPtr>
		class ObservingAQuantity : public ObservationDynamics<TargetParticlesType>
		{
		protected:
			/** Observed quantities saved here. */
			StdLargeVec<DataType>  observed_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				DataType observed_quantity(0);
				Real ttl_weight(0);

				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					TargetParticlesType* target_particles = this->contact_particles_[k];
					StdLargeVec<TargetDataType>& target_data = target_particles->*TrgtDataMemPtr;

					Neighborhood& contact_neighborhood 
						= this->contact_configuration_[k][index_particle_i];
					KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
					CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						CommonRelation& common_relation = contact_common_relations[n];
						size_t index_particle_j = common_relation.j_;
						BaseParticleData& base_particle_data_j
							= this->contact_particles_[k]->base_particle_data_[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;
						TargetDataType& target_data_j = target_data[index_particle_j];

						Real weight_j = kernel_value_list[n] * Vol_j;
						observed_quantity += weight_j * target_data_j.*TrgtMemPtr;
						ttl_weight += weight_j;
					}
				}
				observed_quantities_[index_particle_i] = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit ObservingAQuantity(SPHBodyContactRelation* body_contact_relation)
				: ObservationDynamics<TargetParticlesType>(body_contact_relation) {
				for (size_t i = 0; i < this->body_->number_of_particles_; ++i) observed_quantities_.push_back(DataType(0));
			};
			virtual ~ObservingAQuantity() {};
		};

		/**
		 * @class ObservingADiffusionReactionQuantity
		 * @brief Observing a diffusion-reaction quantity from a body
		 */
		template <class DiffusionReactionParticlesType>
		class ObservingADiffusionReactionQuantity : public ObservationDynamics<DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** Observed quantities saved here. */
			StdLargeVec<Real>  observed_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				Real observed_quantity(0);
				Real ttl_weight(0);
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					DiffusionReactionParticlesType* target_particles = this->contact_particles_[k];
					StdLargeVec<DiffusionReactionData>& target_data = target_particles->diffusion_reaction_data_;

					Neighborhood& contact_neighborhood
						= this->contact_configuration_[k][index_particle_i];
					KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
					CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						CommonRelation& common_relation = contact_common_relations[n];
						size_t index_particle_j = common_relation.j_;
						BaseParticleData& base_particle_data_j
							= this->contact_particles_[k]->base_particle_data_[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;

						Real weight_j = kernel_value_list[n] * Vol_j;
						observed_quantity += weight_j * target_data[index_particle_j].species_n_[species_index_];
						ttl_weight += weight_j;
					}
					observed_quantities_[index_particle_i] = observed_quantity / (ttl_weight + TinyReal);
				}
			};
		public:
			explicit ObservingADiffusionReactionQuantity(string species_name, SPHBodyContactRelation* body_contact_relation)
				: ObservationDynamics<DiffusionReactionParticlesType>(body_contact_relation) {
				species_indexes_map_ = this->contact_particles_[0]->getSpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];
				for (size_t i = 0; i < this->body_->number_of_particles_; ++i) observed_quantities_.push_back(0.0);
			};
			virtual ~ObservingADiffusionReactionQuantity() {};
		};

		/**
		 * @class InterpolatingAQuantity
		 * @brief interpolate a given member data in the particles of a general body from another body
		 */
		template <class DataType, class ObserverParticlesType, class ObserverDataType, class TargetParticlesType, class TargetDataType,
			StdLargeVec<ObserverDataType> ObserverParticlesType:: * ObrsvrDataMemPtr, StdLargeVec<TargetDataType> TargetParticlesType:: * TrgtDataMemPtr,
			DataType ObserverDataType:: * ObrsvrMemPtr, DataType TargetDataType:: * TrgtMemPtr>
			class InterpolatingAQuantity : public InterpolationDynamics<ObserverParticlesType, TargetParticlesType>
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				ObserverParticlesType* particles = this->particles_;
				StdLargeVec<ObserverDataType>& observer_data = this->particles_->*ObrsvrDataMemPtr;
				ObserverDataType& observer_data_i = observer_data[index_particle_i];

				DataType observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					TargetParticlesType* target_particles = this->contact_particles_[k];
					StdLargeVec<BaseParticleData>& target_base_particle_data = target_particles->base_particle_data_;
					StdLargeVec<TargetDataType>& target_data = target_particles->*TrgtDataMemPtr;

					Neighborhood& contact_neighborhood = this->contact_configuration_[k][index_particle_i];
					KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
					CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						CommonRelation& common_relation = contact_common_relations[n];
						size_t index_particle_j = common_relation.j_;
						BaseParticleData& base_particle_data_j = target_base_particle_data[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;
						TargetDataType& target_data_j = target_data[index_particle_j];

						Real weight_j = kernel_value_list[n] * Vol_j;
						observed_quantity += weight_j * target_data_j.*TrgtMemPtr;
						ttl_weight += weight_j;
					}
				}
				observer_data_i.* ObrsvrMemPtr = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit InterpolatingAQuantity(SPHBodyContactRelation* body_contact_relation)
				: InterpolationDynamics<ObserverParticlesType, TargetParticlesType>(body_contact_relation) {};
			virtual ~InterpolatingAQuantity() {};
		};

		/**
		 * @class InterpolatingADiffusionReactionQuantity
		 * @brief interpolate a diffusion-reaction member data in the particles of a general body from another body
		 */
		template <class ObserverParticlesType, class ObserverDataType, class DiffusionReactionParticlesType,
			StdLargeVec<ObserverDataType> ObserverParticlesType:: * ObrsvrDataMemPtr, Real ObserverDataType:: * ObrsvrMemPtr>
			class InterpolatingADiffusionReactionQuantity
			: public InterpolationDynamics<ObserverParticlesType, DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** Observed quantities saved here. */
			StdLargeVec<Real>  observed_quantities_;

			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override
			{
				ObserverParticlesType* particles = this->particles_;
				StdLargeVec<ObserverDataType>& observer_data = this->particles_->*ObrsvrDataMemPtr;
				ObserverDataType& observer_data_i = observer_data[index_particle_i];

				Real observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					DiffusionReactionParticlesType* target_particles = this->contact_particles_[k];
					StdLargeVec<BaseParticleData>& target_base_particle_data = target_particles->base_particle_data_;
					StdLargeVec<DiffusionReactionData>& target_diffusion_reaction_data = target_particles->diffusion_reaction_data_;

					Neighborhood& contact_neighborhood = this->contact_configuration_[k][index_particle_i];
					KernelValueList& kernel_value_list = contact_neighborhood.kernel_value_list_;
					CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						CommonRelation& common_relation = contact_common_relations[n];
						size_t index_particle_j = common_relation.j_;
						BaseParticleData& base_particle_data_j = target_base_particle_data[index_particle_j];
						Real Vol_j = base_particle_data_j.Vol_;
						DiffusionReactionData& target_data_j = target_diffusion_reaction_data[index_particle_j];

						Real weight_j = kernel_value_list[n] * Vol_j;
						observed_quantity += weight_j * target_data_j.species_n_[species_index_];
						ttl_weight += weight_j;
					}
				}
				observer_data_i.*ObrsvrMemPtr = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit InterpolatingADiffusionReactionQuantity(string species_name, SPHBodyContactRelation* body_contact_relation)
				: InterpolationDynamics<ObserverParticlesType, DiffusionReactionParticlesType>(body_contact_relation)
			{
				species_indexes_map_ = this->contact_particles_[0]->getSpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];
			};
			virtual ~InterpolatingADiffusionReactionQuantity() {};
		};
		
		/**
		* @class CorrectInterpolationKernelWeights
		* @brief  correct kernel weights for interpolation between general bodies
		*/
		class CorrectInterpolationKernelWeights : 
			public ParticleDynamicsContact<SPHBody, BaseParticles, BaseMaterial, SPHBody, BaseParticles, BaseMaterial>
		{
		protected:
			virtual void ContactInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CorrectInterpolationKernelWeights(SPHBodyContactRelation* body_contact_relation)
				: ParticleDynamicsContact<SPHBody, BaseParticles, 
				BaseMaterial, SPHBody, BaseParticles, BaseMaterial>(body_contact_relation) {};
			virtual ~CorrectInterpolationKernelWeights() {};
		};
	}
}
