/**
* @file mesh_cell_linked_list.h
* @brief Here gives the classes for managing cell linked lists. This is the basic class 
* for building the particle configurations.
* @details  The cell linked list saves for each body a list of particles
* located within the cell.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_mesh.h"

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
		/** the index vector for iterate particles in a split scheme. */
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
		SPHBody* body_;
		BaseParticles* base_particles_;
		Kernel* kernel_;

		/** clear the cell lists */
		void ClearCellLists(Vecu& number_of_cells, matrix_cell cell_linked_lists);
		/** clear split cell lists in this mesh*/
		void ClearSplitCellLists(SplitCellLists& split_cell_lists);
		/** update split particle list in this mesh */
		void UpdateSplitCellLists(SplitCellLists& split_cell_lists,
			Vecu& number_of_cells, matrix_cell cell_linked_lists);
		/** update cell linked list data in this mesh */
		void UpdateCellListData(matrix_cell cell_linked_lists);
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		BaseMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound, 
			Real cell_spacing, size_t buffer_size = 2);
		/** Constructor with the direct information of the mesh. */
		BaseMeshCellLinkedList(SPHBody* body, 
			Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~BaseMeshCellLinkedList() {};

		/** computing search range for building contact configuration */
		int ComputingSearchRage(int origin_refinement_level,
			int target_refinement_level);
		/** choose a kernel for building up inter refinement level configuration */
		Kernel& ChoosingKernel(Kernel* original_kernel, Kernel* target_kernel);
		/** get the address of cell list */
		virtual CellList* CellListFormIndex(Vecu cell_index) = 0;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() = 0;

		/** Assign base particles to the mesh cell linked list. */
		void assignParticles(BaseParticles* base_particles);
		/** Assign kernel to the mesh cell linked list. */
		void reassignKernel(Kernel* kernel);
		/** allcate memories for mesh data */
		virtual void allocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() = 0;

		/** update the cell lists */
		virtual void UpdateCellLists() = 0;

		/** Insert a cell-linked_list entry. */
		virtual void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) = 0;
		virtual void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) = 0;
	};

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
	 * @class MeshCellLinkedList
	 * @brief Defining a mesh cell linked list for a body.
	 * The meshes for all bodies share the same global coordinates.
	 */
	class MeshCellLinkedList : public BaseMeshCellLinkedList
	{
	protected:
		/** cut_off radius */
		Real cutoff_radius_;
		/** The array for of mesh cells, i.e. mesh data.
		 * Within each cell, a list is saved with the indexes of particles.*/
		matrix_cell cell_linked_lists_;
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		MeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real cell_spacing, size_t buffer_size = 2);
		/** direct construct with mesh information. */
		MeshCellLinkedList(SPHBody* body, Vecd mesh_lower_bound,
			Vecu number_of_cells, Real cell_spacing);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~MeshCellLinkedList() { deleteMeshDataMatrix(); };

		/** access protected members */
		virtual CellList* CellListFormIndex(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() override { return cell_linked_lists_; };

		/** allcate memories for mesh data */
		virtual void allocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;

		/** output mesh data for visualization */
		virtual void writeMeshToVtuFile(ofstream &output_file) override {};
		virtual void writeMeshToPltFile(ofstream &output_file) override {};

		/** Insert a cell-linked_list entry. */
		void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override;
	};

	/**
	  * @class MultilevelMeshCellLinkedList
	  * @brief Defining a multimesh cell linked list for a body
	  * for multiresolution particle configuration.
	  */
	class MultilevelMeshCellLinkedList : public BaseMeshCellLinkedList
	{
	protected:
		/** total levels of the storage pyramid.
		  * it is an odd value so that the reference level is the middle level.*/
		size_t total_levels_;
		/**cell spacing for the mesh levels*/
		StdVec<Real> cell_spacing_levels_;
		/** number of cells by dimension */
		StdVec<Vecu> number_of_cells_levels_;
		/** split cell list for building configuration.*/
		StdVec<SplitCellLists> split_cell_lists_levels_;
		/** Projected cell linked lists for building configuration.*/
		StdVec<MeshDataMatrix<CellList>> cell_linked_lists_levels_;
		/** point to every mesh level. */
		StdVec<MeshCellLinkedList*> mesh_cell_linked_list_levels_;

		/** determine mesh level of a particle. */
		size_t getLevelFromCutOffRadius(Real smoothing_length);
		/** find cell indexes from point position in a level */
		Vecu getLevelCellIndexesFromPosition(Vecd& position, size_t level);
		/** Insert a cell-linked_list entry to the projected particle list. */
		void InsertACellLinkedListEntryAtALevel(size_t particle_index, Vecd& position, Vecu& cell_index, size_t level);
	public:
		/** Constructor to achieve alignment of all mesh levels. */
		MultilevelMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real reference_cell_spacing, size_t total_levels = 1, size_t buffer_size = 2);
		/**In the destructor, the dynamically located memory is released.*/
		virtual ~MultilevelMeshCellLinkedList() {};

		/** access protected members */
		virtual CellList* CellListFormIndex(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell CellLinkedLists() override { return cell_linked_lists_levels_[0]; };

		/** allcate memories for mesh data */
		virtual void allocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void deleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;

		/** Insert a cell-linked_list entry to the projected particle list. */
		void InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position) override;
		void InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position) override {};
	};
}
