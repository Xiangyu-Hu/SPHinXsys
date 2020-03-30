/**
 * @file 	observer_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "observer_dynamics.h"
//=================================================================================================//
using namespace SimTK;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace observer_dynamics
	{
//=================================================================================================//
		template <typename FluidQuantityType>
		void ObserveAFluidQuantity<FluidQuantityType>
			::ContactInteraction(size_t index_particle_i,	size_t interacting_body_index, Real dt)
		{
			FluidQuantityType observed_quantity(0);
			Real weight(0);
			NeighborList& contact_neighors
				= getNeighborList(current_interacting_configuration_[interacting_body_index], index_particle_i);
			for (size_t n = 0; n < contact_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = contact_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j 
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAFluidQuantity(index_particle_j, *interacting_particles_[interacting_body_index]) 
					* neighboring_particle->W_ij_ * base_particle_data_j.Vol_ ;
				weight += neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
			}

			fluid_quantities_[index_particle_i] = observed_quantity / weight;
		}
		//template definitions should be instantiated here
		template class ObserveAFluidQuantity<Real>;
		template class ObserveAFluidQuantity<Vecd>;
//=================================================================================================//
		template <typename ElasticSolidQuantityType>
		void ObserveAnElasticSolidQuantity<ElasticSolidQuantityType>
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			ElasticSolidQuantityType observed_quantity(0);
			Real total_weight(0);
			NeighborList& contact_neighors
				= getNeighborList(current_interacting_configuration_[interacting_body_index], index_particle_i);
			for (size_t n = 0; n < contact_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = contact_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAnElasticSolidQuantity(index_particle_j, *interacting_particles_[interacting_body_index])
					* neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
			}

			elastic_body_quantities_[index_particle_i] = observed_quantity / total_weight;
		}
		//template definitions should be instantiated here
		template class ObserveAnElasticSolidQuantity<Real>;
		template class ObserveAnElasticSolidQuantity<Vecd>;
//=================================================================================================//
		template <typename ElectroPhysiologyQuantityType>
		void ObserveAElectroPhysiologyQuantity<ElectroPhysiologyQuantityType>
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			ElectroPhysiologyQuantityType observed_quantity(0);
			Real total_weight(0);
			NeighborList& contact_neighors
				= getNeighborList(current_interacting_configuration_[interacting_body_index], index_particle_i);
			for (size_t n = 0; n < contact_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = contact_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAMuscleQuantity(index_particle_j, *interacting_particles_[interacting_body_index])
					* neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
			}

			electro_physiology_quantities_[index_particle_i] = observed_quantity / total_weight;
		}
		//template definitions should be instantiated here
		template class ObserveAElectroPhysiologyQuantity<Real>;
//=================================================================================================//
		template <typename ElectroPhysiologyQuantityType>
		void ElectroPhysiologyQuantityInterpolation<ElectroPhysiologyQuantityType>::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			BaseParticleData &base_particle_data_i 		= particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i 			= particles_->solid_body_data_[index_particle_i];
			ElasticSolidParticleData &elastic_data_i 	= particles_->elastic_body_data_[index_particle_i];
			ActiveMuscleData &active_muscle_data_i 		= particles_->active_muscle_data_[index_particle_i];

			ElectroPhysiologyQuantityType observed_quantity(0);
			Real total_weight(0);
			NeighborList& contact_neighors
				= getNeighborList(current_interacting_configuration_[interacting_body_index], index_particle_i);
			for (size_t n = 0; n < contact_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = contact_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData &base_particle_data_j = (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAMuscleQuantity(index_particle_j, *interacting_particles_[interacting_body_index])
					* neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle->W_ij_ * base_particle_data_j.Vol_;
			}

			active_muscle_data_i.active_contraction_stress_ = observed_quantity / total_weight;
		}
		//template definitions should be instantiated here
		template class ElectroPhysiologyQuantityInterpolation<Real>;
//=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//


