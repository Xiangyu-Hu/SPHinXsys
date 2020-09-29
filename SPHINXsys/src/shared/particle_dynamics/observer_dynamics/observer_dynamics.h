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
		using InterpolationDataDelegate
			= DataDelegateContact<SPHBody, ObserverParticlesType, BaseMaterial, SPHBody, TargetParticlesType>;

		template <class TargetParticlesType>
		using ObservationDataDelegate
			= DataDelegateContact<SPHBody, BaseParticles, BaseMaterial, SPHBody, TargetParticlesType, BaseMaterial>;

		/**
		 * @class ObservingAQuantity
		 * @brief Observering a given member data in the particles of a general body
		 */
		template <class DataType, class TargetParticlesType, StdLargeVec<DataType> TargetParticlesType:: * TrgtMemPtr>
		class ObservingAQuantity : 
			public ParticleDynamicsContact, public ObservationDataDelegate<TargetParticlesType>
		{
		protected:
			/** Observed quantities saved here. */
			StdLargeVec<DataType>  observed_quantities_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			StdVec<StdLargeVec<DataType>*> contact_data_;

			virtual void ContactInteraction(size_t index_i, Real dt = 0.0) override
			{
				DataType observed_quantity(0);
				Real ttl_weight(0);

				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<DataType>& data_k = *(contact_data_[k]);
					Neighborhood& contact_neighborhood
						= this->contact_configuration_[k][index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

						observed_quantity += weight_j * data_k[index_j];
						ttl_weight += weight_j;
					}
				}
				observed_quantities_[index_i] = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit ObservingAQuantity(SPHBodyContactRelation* body_contact_relation)
				: ParticleDynamicsContact(body_contact_relation),
				ObservationDataDelegate<TargetParticlesType>(body_contact_relation)
			{
				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
					contact_data_.push_back(&(this->contact_particles_[k]->*TrgtMemPtr));
				}
				for (size_t i = 0; i < this->body_->number_of_particles_; ++i) 
					observed_quantities_.push_back(DataType(0));
			};
			virtual ~ObservingAQuantity() {};
		};

		/**
		 * @class ObservingADiffusionReactionQuantity
		 * @brief Observing a diffusion-reaction quantity from a body
		 */
		template <class DiffusionReactionParticlesType>
		class ObservingADiffusionReactionQuantity :
			public ParticleDynamicsContact, public ObservationDataDelegate<DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** Observed quantities saved here. */
			StdLargeVec<Real>  observed_quantities_;
			StdVec<StdLargeVec<Real>*> contact_Vol_, contact_data_;

			virtual void ContactInteraction(size_t index_i, Real dt = 0.0) override
			{
				Real observed_quantity(0);
				Real ttl_weight(0);
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Real>& data_k = *(contact_data_[k]);
					Neighborhood& contact_neighborhood
						= this->contact_configuration_[k][index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

						observed_quantity += weight_j * data_k[index_j];
						ttl_weight += weight_j;
					}
					observed_quantities_[index_i] = observed_quantity / (ttl_weight + TinyReal);
				}
			};
		public:
			explicit ObservingADiffusionReactionQuantity(string species_name, SPHBodyContactRelation* body_contact_relation)
				: ParticleDynamicsContact(body_contact_relation),
				ObservationDataDelegate<DiffusionReactionParticlesType>(body_contact_relation)
			{
				species_indexes_map_ = this->contact_particles_[0]->SpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];
				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
					contact_data_.push_back(&(this->contact_particles_[k]->species_n_[species_index_]));
				}

				for (size_t i = 0; i < this->body_->number_of_particles_; ++i) 
					observed_quantities_.push_back(0.0);
			};
			virtual ~ObservingADiffusionReactionQuantity() {};
		};

		/**
		 * @class InterpolatingAQuantity
		 * @brief interpolate a given member data in the particles of a general body from another body
		 */
		template <class DataType, class ObserverParticlesType, class TargetParticlesType, 
			StdLargeVec<DataType> ObserverParticlesType:: * ObrsvrMemPtr, StdLargeVec<DataType> TargetParticlesType:: * TrgtMemPtr>
		class InterpolatingAQuantity : 
			public ParticleDynamicsContact,
			public InterpolationDataDelegate<ObserverParticlesType, TargetParticlesType>
		{
		protected:
			/** interpolated quantities saved here. */
			StdLargeVec<DataType>  &interpolated_data_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			StdVec<StdLargeVec<DataType>*> contact_data_;

			virtual void ContactInteraction(size_t index_i, Real dt = 0.0) override
			{
				DataType observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<DataType>& data_k = *(contact_data_[k]);
					Neighborhood& contact_neighborhood = this->contact_configuration_[k][index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

						observed_quantity += weight_j * data_k[index_j];
						ttl_weight += weight_j;
					}
				}

				interpolated_data_[index_i] = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit InterpolatingAQuantity(SPHBodyContactRelation* body_contact_relation)
				: ParticleDynamicsContact(body_contact_relation),
				InterpolationDataDelegate<ObserverParticlesType, TargetParticlesType>(body_contact_relation),
				interpolated_data_(this->particles_->*ObrsvrMemPtr)
			{
				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
					contact_data_.push_back(&(this->contact_particles_[k]->*TrgtMemPtr));
				}
			};
			virtual ~InterpolatingAQuantity() {};
		};

		/**
		 * @class InterpolatingADiffusionReactionQuantity
		 * @brief interpolate a diffusion-reaction member data in the particles of a general body from another body
		 */
		template <class ObserverParticlesType, class DiffusionReactionParticlesType,
			StdLargeVec<Real> ObserverParticlesType:: * ObrsvrDataMemPtr>
		class InterpolatingADiffusionReactionQuantity : 
			public ParticleDynamicsContact,
			public InterpolationDataDelegate<ObserverParticlesType, DiffusionReactionParticlesType>
		{
		protected:
			/** Index of voltage. */
			size_t species_index_;
			map<string, size_t> species_indexes_map_;
			/** interpolated quantities saved here. */
			StdLargeVec<Real>& interpolated_data_;
			StdVec<StdLargeVec<Real>*> contact_Vol_, contact_data_;

			virtual void ContactInteraction(size_t index_i, Real dt = 0.0) override
			{

				Real observed_quantity(0);
				Real ttl_weight(0);
				/** Compute the first order consistent kernel weights */
				for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
				{
					StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
					StdLargeVec<Real>& data_k = *(contact_data_[k]);
					Neighborhood& contact_neighborhood = this->contact_configuration_[k][index_i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						Real weight_j = contact_neighborhood.W_ij_[n] * Vol_k[index_j];

						observed_quantity += weight_j * data_k[index_j];
						ttl_weight += weight_j;
					}
				}
				interpolated_data_[index_i] = observed_quantity / (ttl_weight + TinyReal);
			};
		public:
			explicit InterpolatingADiffusionReactionQuantity(string species_name, SPHBodyContactRelation* body_contact_relation)
				: ParticleDynamicsContact(body_contact_relation),
				InterpolationDataDelegate<ObserverParticlesType, DiffusionReactionParticlesType>(body_contact_relation),
				interpolated_data_(this->particles_->*ObrsvrDataMemPtr)
			{
				species_indexes_map_ = this->contact_particles_[0]->SpeciesIndexMap();
				species_index_ = species_indexes_map_[species_name];

				for (size_t k = 0; k != this->contact_particles_.size(); ++k)
				{
					contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
					contact_data_.push_back(&(this->contact_particles_[k]->species_n_[species_index_]));
				}
			};
			virtual ~InterpolatingADiffusionReactionQuantity() {};
		};
		
		/**
		* @class CorrectInterpolationKernelWeights
		* @brief  correct kernel weights for interpolation between general bodies
		*/
		class CorrectInterpolationKernelWeights : 
			public ParticleDynamicsContact,
			public DataDelegateContact<SPHBody, BaseParticles, BaseMaterial, SPHBody, BaseParticles, BaseMaterial>
		{
		public:
			CorrectInterpolationKernelWeights(SPHBodyContactRelation* body_contact_relation);
			virtual ~CorrectInterpolationKernelWeights() {};
		protected:
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			virtual void ContactInteraction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
