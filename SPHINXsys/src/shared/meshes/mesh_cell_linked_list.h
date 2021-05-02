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
* @file mesh_cell_linked_list.h
* @brief Here gives the classes for managing cell linked lists. This is the basic class 
* for building the particle configurations.
* @details  The cell linked list saves for each body a list of particles
* located within the cell.
* @author	Yongchuan Yu, Chi ZHang and Xiangyu Hu
*/

#pragma once

#include "base_mesh.h"
#include "neighbor_relation.h"

namespace SPH {

	class SPHSystem;
	class SPHBody;
	class BaseParticles;
	class Kernel;

	/**
	 * @class CellList
	 * @brief The linked list for one cell
	 */
	class CellList
	{
	public:
		/** using concurrent vectors due to writting conflicts when building the list */
		ConcurrentIndexVector concurrent_particle_indexes_;
		/** non-concurrent cell linked list rewritten for building neighbor list */
		CellListDataVector cell_list_data_;
		/** the index vector for real particles. */
		IndexVector real_particle_indexes_;

		CellList();
		~CellList() {};
	};

	/**
	 * @class BaseMeshCellLinkedList
	 * @brief Abstract class for mesh cell linked list.
	 */
	class BaseMeshCellLinkedList : public Mesh
	{
	protected:
		SPHBody& sph_body_;
		Kernel& kernel_;
		BaseParticles* base_particles_;

		/** clear split cell lists in this mesh*/
		virtual void clearSplitCellLists(SplitCellLists& split_cell_lists);
		/** update split particle list in this mesh */
		virtual void updateSplitCellLists(SplitCellLists& split_cell_lists) = 0;
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		BaseMeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation, 
			BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width = 2);
		virtual ~BaseMeshCellLinkedList() {};

		/** Assign base particles to the mesh cell linked list,
		 * and is important because particles are not defined in the constructor.  */
		virtual void assignBaseParticles(BaseParticles* base_particles) = 0;

		/** update the cell lists */
		virtual void UpdateCellLists() = 0;
		/** Insert a cell-linked_list entry to the concurrent index list. */
		virtual void insertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) = 0;
		/** Insert a cell-linked_list entry of the index and particle position pair. */
		virtual void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) = 0;
		/** find the nearest list data entry */
		virtual ListData findNearestListDataEntry(Vecd& position) = 0;
		/** computing the sequence which indicate the order of sorted particle data */
		virtual void computingSequence(StdLargeVec<size_t>& sequence) = 0;
		/** Tag body part by cell, call by body part */
		virtual void tagBodyPartByCell(CellLists& cell_lists, std::function<bool(Vecd, Real)>& check_included) = 0;
		/** Tag domain bounding cells in an axis direction, called by domain boudning classes */
		virtual void tagBodyDomainBoundingCells(StdVec<CellLists>& cell_lists, BoundingBox& body_domain_bounds, int axis) = 0;
		/** Tag mirror bounding cells, called by mirror boundary condition */
		virtual void tagMirrorBoundingCells(CellLists& cell_lists, BoundingBox& body_domain_bounds, int axis, bool positive) = 0;
	};

	/**
	 * @class MeshCellLinkedList
	 * @brief Defining a mesh cell linked list for a body.
	 * The meshes for all bodies share the same global coordinates.
	 */
	class MeshCellLinkedList : public BaseMeshCellLinkedList
	{
	protected:
		/** The array for of mesh cells, i.e. mesh data.
		 * Within each cell, a list is saved with the indexes of particles.*/
		matrix_cell cell_linked_lists_;

		virtual void updateSplitCellLists(SplitCellLists& split_cell_lists) override;
	public:
		MeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation, 
			BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width = 2);
		virtual ~MeshCellLinkedList() { deleteMeshDataMatrix(); };

		virtual void allocateMeshDataMatrix() override;
		virtual void deleteMeshDataMatrix() override;
		virtual void assignBaseParticles(BaseParticles* base_particles) override;

		void clearCellLists();
		void UpdateCellListData();
		virtual void UpdateCellLists() override;
		void insertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override;
		virtual ListData findNearestListDataEntry(Vecd& position) override;
		virtual void computingSequence(StdLargeVec<size_t>& sequence) override;
		virtual void tagBodyPartByCell(CellLists& cell_lists, std::function<bool(Vecd, Real)>& check_included) override;
		virtual void tagBodyDomainBoundingCells(StdVec<CellLists>& cell_lists, BoundingBox& body_domain_bounds, int axis) override;
		virtual void tagMirrorBoundingCells(CellLists& cell_lists, BoundingBox& body_domain_bounds, int axis, bool positive) override;
		virtual void writeMeshToPltFile(std::ofstream& output_file) override;

		/** generalized particle search algorithm */
		template<typename GetParticleIndex, typename GetSearchRange, typename GetNeighborRelation>
		void searchNeighborsByParticles(size_t total_real_particles, BaseParticles& source_particles, 
			ParticleConfiguration& particle_configuration, GetParticleIndex& get_particle_index,
			GetSearchRange& get_search_range, GetNeighborRelation& get_neighbor_relation);
	};

	/**
	  * @class MultilevelMeshCellLinkedList
	  * @brief Defining a multilevel mesh cell linked list for a body
	  * for multiresolution particle configuration.
	  */
	class MultilevelMeshCellLinkedList : 
		public MultilevelMesh<SPHBody, BaseMeshCellLinkedList, MeshCellLinkedList>
	{
	protected:
		StdLargeVec<Real>& h_ratio_;	
		virtual void updateSplitCellLists(SplitCellLists& split_cell_lists) override {};
		/** determine mesh level from particle cutoff radius */
		inline size_t getMeshLevel(Real particle_cutoff_radius);
	public:
		MultilevelMeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation, 
			BoundingBox tentative_bounds, Real reference_grid_spacing, 
			size_t total_levels, Real maximum_spacing_ratio);
		virtual ~MultilevelMeshCellLinkedList() {};

		virtual void assignBaseParticles(BaseParticles* base_particles) override;
		virtual void UpdateCellLists() override;
		void insertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override;
		virtual ListData findNearestListDataEntry(Vecd& position) override { return ListData(0, Vecd(0)); };
		virtual void computingSequence(StdLargeVec<size_t>& sequence) override {};
		virtual void tagBodyPartByCell(CellLists& cell_lists, std::function<bool(Vecd, Real)>& check_included) override;
		virtual void tagBodyDomainBoundingCells(StdVec<CellLists>& cell_lists, BoundingBox& body_domain_bounds, int axis) override {};
		virtual void tagMirrorBoundingCells(CellLists& cell_lists, BoundingBox& body_domain_bounds, int axis, bool positive) override {};
	};
}
