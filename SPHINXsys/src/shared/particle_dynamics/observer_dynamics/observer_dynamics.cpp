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
	namespace observer_dynamics
	{
		//===============================================================//
		void ObserveFluidPressure::ContactInteraction(size_t index_particle_i,	size_t interacting_body_index, Real dt)
		{
			Real pressure(0);
			Real weight(0);
			StdVec<NeighboringParticle>  &neighors 
				= (*current_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j 
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];
				WeaklyCompressibleFluidParticleData &fluid_data_j 
					= (*interacting_particles_[interacting_body_index]).fluid_data_[index_particle_j];

				pressure += fluid_data_j.p_ * neighboring_particle.W_ij_ * base_particle_data_j.Vol_ ;
				weight += neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
			}

			pressures_[index_particle_i] = pressure / weight;
		}
		//===============================================================//
		void ObserveElasticDisplacement::ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt)
		{
			Vecd displacement(0);
			Real total_weight(0);
			StdVec<ReferenceNeighboringParticle>  &neighors
				= (*reference_interacting_configuration_[interacting_body_index])[index_particle_i];
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				ReferenceNeighboringParticle &neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle.j_;
				BaseParticleData &base_particle_data_j
					= (*interacting_particles_[interacting_body_index]).base_particle_data_[index_particle_j];

				displacement += base_particle_data_j.pos_n_ * neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
				total_weight += neighboring_particle.W_ij_ * base_particle_data_j.Vol_;
			}

			displacements_[index_particle_i] = displacement / total_weight;
		}
	}
}



