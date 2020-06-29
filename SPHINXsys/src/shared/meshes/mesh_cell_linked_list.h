/**
* @file mesh_cell_linked_list.h
* @brief Here gives the classes for managing cell linked lists. This is the baic class 
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
		/** the saved lists currently using concurrent vectors
		 * due to writting conflicts when building the lists */
		ConcurrentListDataVector particle_data_lists_;
		/** the index vector for itreate particles in a split scheme. */
		IndexVector real_particle_indexes_;
		Vecu cell_location_;
		size_t real_particle_count_;

		CellList();
		~CellList() {};

		/**Initialize cell information. */
		void setCellInformation(Vecu cell_location) {
			cell_location_ = cell_location;
		};
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
		SPHBodyContactMap contact_map_;
		Kernel* kernel_;

		/** clear the cell lists */
		void ClearCellLists(Vecu& number_of_cells, matrix_cell cell_linked_lists);
		/** clear split cell lists in this mesh*/
		void ClearSplitCellLists(SplitCellLists& split_cell_lists);
		/** update split particle list in this mesh */
		void UpdateSplitCellLists(SplitCellLists& split_cell_lists,
			Vecu& number_of_cells, matrix_cell cell_linked_lists);
	public:
		/** The buffer size 2 used to expand computational domian for particle searching. */
		BaseMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound, 
			Real cell_spacing, size_t buffer_size = 2);
		/** Constructor with the direct information of the mesh. */
		BaseMeshCellLinkedList(SPHBody* body, 
			Vecd mesh_lower_bound, Vecu number_of_cells, Real cell_spacing);
		/**In the destructor, the dynamically located memeory is released.*/
		virtual ~BaseMeshCellLinkedList() {};

		/** computing search range for building contact configuration */
		int ComputingSearchRage(int orign_refinement_level,
			int target_refinement_level);
		/** choose a kernel for building up inter refinement level configuration */
		Kernel& ChoosingKernel(Kernel* orignal_kernel, Kernel* target_kernel);
		/** get the address of cell list */
		virtual CellList* getCellList(Vecu cell_index) = 0;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell getCellLinkedLists() = 0;

		/** Assign base particles to the mesh cell linked list. */
		void assignParticles(BaseParticles* base_particles);
		/** Assign contact map to the mesh cell linked list. */
		void assignContactMap(SPHBodyContactMap contact_map);
		/** Assign kernel to the mesh cell linked list. */
		void reassignKernel(Kernel* kernel);
		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() = 0;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() = 0;

		/** update the cell lists */
		virtual void UpdateCellLists() = 0;

		/** build reference inner configurtion */
		virtual void BuildInnerConfiguration(ParticleConfiguration& inner_configuration);
		/** build reference contact configuration */
		virtual void BuildContactConfiguration();
		/** update contact configuration */
		virtual void UpdateContactConfiguration();

		/** update inner configuration */
		virtual void UpdateInnerConfiguration(ParticleConfiguration& inner_configuration) = 0;
		/** update interaction configuration */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) = 0;

		/** Insert a cell-linked_list entry. */
		virtual void InsertACellLinkedListEntry(size_t particle_index, Vecd particle_position) = 0;
	};

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
		/**In the destructor, the dynamically located memeory is released.*/
		virtual ~MeshCellLinkedList() { DeleteMeshDataMatrix(); };

		/** access protected members */
		virtual CellList* getCellList(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell getCellLinkedLists() override { return cell_linked_lists_; };

		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;

		/** update inner configuration */
		virtual void UpdateInnerConfiguration(ParticleConfiguration& inner_configuration) override;
		/** update interaction configuration */
		virtual void UpdateInteractionConfiguration(SPHBodyVector interacting_bodies) override;

		/** output mesh data for visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) override {};
		virtual void WriteMeshToPltFile(ofstream &output_file) override {};

		/** Insert a cell-linked_list entry. */
		void InsertACellLinkedListEntry(size_t particle_index, Vecd particle_position) override;
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
		/** split cell list for building confifuration.*/
		StdVec<SplitCellLists> split_cell_lists_levels_;
		/** Projected cell linked lists for building configuration.*/
		StdVec<MeshDataMatrix<CellList>> cell_linked_lists_levels_;
		/** point to every mesh level. */
		StdVec<MeshCellLinkedList*> mesh_cell_linked_list_levels_;

		/** determine mesh level of a particle. */
		size_t getLevelFromCutOffRadius(Real smoothing_length);
		/** find cell indexes from point poistion in a level */
		Vecu getLevelCellIndexesFromPosition(Vecd& position, size_t level);
		/** Insert a cell-linked_list entry to the preojected particle list. */
		void InsertACellLinkedListEntryAtALevel(size_t particle_index, Vecd& position, Vecu& cell_index, size_t level);
	public:
		/** Constructor to achieve alignment of all mesh levels. */
		MultilevelMeshCellLinkedList(SPHBody* body, Vecd lower_bound, Vecd upper_bound,
			Real reference_cell_spacing, size_t total_levels = 1, size_t buffer_size = 2);
		/**In the destructor, the dynamically located memeory is released.*/
		virtual ~MultilevelMeshCellLinkedList() {};

		/** access protected members */
		virtual CellList* getCellList(Vecu cell_index) override;
		/** Get the array for of mesh cell linked lists.*/
		virtual matrix_cell getCellLinkedLists() override { return cell_linked_lists_levels_[0]; };

		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() override;

		/** update the cell lists */
		virtual void UpdateCellLists() override;
		/** update inner configuration */
		virtual void UpdateInnerConfiguration(ParticleConfiguration& inner_configuration) override;

		/** Insert a cell-linked_list entry to the preojected particle list. */
		void InsertACellLinkedListEntry(size_t particle_index, Vecd particle_position) override;
	};
}
