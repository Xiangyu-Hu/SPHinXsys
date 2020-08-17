/**
 * @file 	body_relation.h
 * @brief 	The topological relations between bodies are described here.
 * @author	Xiangyu Hu
 */

#pragma once
#include "neighbor_relation.h"
#include "base_body.h"

namespace SPH
{
	/**
	 * @class SPHBodyBaseRelation
	 * @brief The relation within a SPH body or with its contact SPH bodies
	 */
	class SPHBodyBaseRelation
	{
	public:
		SPHBody* body_;
		SplitCellLists& split_cell_lists_;
		BaseParticles* base_particles_;
		BaseMeshCellLinkedList* base_mesh_cell_linked_list_;

		SPHBodyBaseRelation(SPHBody* body);
		virtual ~SPHBodyBaseRelation() {};

		void subscribe_to_body() { body_->body_relations_.push_back(this); };
		virtual void updateConfigurationMemories() = 0;
		virtual void updateConfiguration() = 0;
	protected:
		virtual void createNeighborRelation(Neighborhood& neighborhood,
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);
		virtual void initializeNeighborRelation(Neighborhood& neighborhood, size_t current_count_of_neighbors, 
			Kernel& kernel, Vecd& vec_r_ij, size_t i_index, size_t j_index);

	};

	/**
	 * @class SPHBodyInnerRelation
	 * @brief The relation within a SPH body
	 */
	class SPHBodyInnerRelation : public SPHBodyBaseRelation
	{
	public:
		/** inner configuration for the neighbor relations. */
		ParticleConfiguration inner_configuration_;

		SPHBodyInnerRelation(SPHBody* body);
		virtual ~SPHBodyInnerRelation() {};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SPHBodyContactRelation
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class SPHBodyContactRelation : public SPHBodyBaseRelation
	{
	protected:
		StdVec<BaseMeshCellLinkedList*> target_mesh_cell_linked_lists_;
		virtual bool checkNeighbor(Real particle_distance, Real cutoff_radius,
			BaseParticleData& base_particle_data_i, BaseParticleData& base_particle_data_j);

	public:
		SPHBodyVector relation_bodies_;

		/** Configurations for particle interaction between bodies. */
		ContatcParticleConfiguration contact_configuration_;

		SPHBodyContactRelation(SPHBody* body, SPHBodyVector relation_bodies);
		virtual ~SPHBodyContactRelation() {};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};

	/**
	 * @class SPHBodyCollisionRelation
	 * @brief The relation between a SPH body and its contact SPH bodies
	 */
	class SPHBodyCollisionRelation : public SPHBodyContactRelation
	{
	protected:
		virtual bool checkNeighbor(Real particle_distance, Real cutoff_radius,
			BaseParticleData& base_particle_data_i, BaseParticleData& base_particle_data_j) override;

	public:
		SPHBodyCollisionRelation(SPHBody* body,	SPHBodyVector relation_bodies)
			: SPHBodyContactRelation(body, relation_bodies) {};
		virtual ~SPHBodyCollisionRelation() {};
	};

	/**
	 * @class SPHBodyComplexRelation
	 * @brief The relation within a SPH body and with its contact SPH bodies.
	 * The interaction is in a inner-boundary-condition fashion. Here inner interaction is
	 * different from conact interaction
	 */
	class SPHBodyComplexRelation : public SPHBodyBaseRelation
	{
	protected:
		SPHBodyInnerRelation* inner_relation_;
		SPHBodyContactRelation* contact_relation_;
	public:
		SPHBodyVector relation_bodies_;

		/** inner configuration for the neighbor relations. */
		ParticleConfiguration& inner_configuration_;
		/** Configurations for updated Lagrangian formulation. **/
		ContatcParticleConfiguration& contact_configuration_;

		SPHBodyComplexRelation(SPHBody* body, SPHBodyVector contact_bodies);
		SPHBodyComplexRelation(SPHBodyInnerRelation* body_inner_relation, SPHBodyVector contact_bodies);
		virtual ~SPHBodyComplexRelation() {
			delete inner_relation_;
			delete contact_relation_;
		};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration()  override;
	};
}
