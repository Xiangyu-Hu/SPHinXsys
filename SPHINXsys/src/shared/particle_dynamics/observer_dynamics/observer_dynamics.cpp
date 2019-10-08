#include "observer_dynamics.h"
#include "all_particles.h"
#include "all_types_of_bodies.h"
#include "all_materials.h"
#include "neighboring_particle.h"
#include "mesh_cell_linked_list.h"
#include "elastic_solid.h"

using namespace SimTK;

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
			StdVec<NeighboringParticle>  &neighors 
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j 
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAFluidQuantity(index_particle_j, *interacting_particles_[interacting_body_index]) 
					* neighboring_particle.W_ij_ * base_particle_data_j.Vol_ ;
				weight += neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
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
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAnElasticSolidQuantity(index_particle_j, *interacting_particles_[interacting_body_index])
					* neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
			}

			elastic_body_quantities_[index_particle_i] = observed_quantity / total_weight;
		}
		//template definitions should be instantiated here
		template class ObserveAnElasticSolidQuantity<Real>;
		template class ObserveAnElasticSolidQuantity<Vecd>;
//=================================================================================================//
		template <typename MuscleQuantityType>
		void ObserveAMuscleQuantity<MuscleQuantityType>
			::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			MuscleQuantityType observed_quantity(0);
			Real total_weight(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				observed_quantity += GetAMuscleQuantity(index_particle_j, *interacting_particles_[interacting_body_index])
					* neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
			}

			muscle_quantities_[index_particle_i] = observed_quantity / total_weight;
		}
		//template definitions should be instantiated here
		template class ObserveAMuscleQuantity<Real>;
//=================================================================================================//
	}
//=================================================================================================//
}



