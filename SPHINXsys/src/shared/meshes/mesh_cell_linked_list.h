/**
* @file mesh_cell_linked_list.h
* @brief This is the classes for managing cell linked lists. This is the baic class 
* for building the particle configurations.
* @details  The cell linked list saves for each body a list of particles
* located within the cell
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "sph_data_conainers.h"
#include "base_mesh.h"

namespace SPH {

	class SPHSystem;
	class SPHBody;
	class Kernel;

	/**
	 * @class ListData
	 * @brief Save particle index and other fast access information 
	 */
	class ListData
	{
	public:
		size_t particle_index_;
		Vecd particle_position_;

		ListData(size_t particle_index, Vecd &particle_position) { 
			particle_index_ = particle_index;
			particle_position_ = particle_position;
		};
		~ListData() {};
	};

	/**
	 * @class CellList
	 * @brief The linked list for one cell
	 */
	class CellList
	{

	public:
		CellList();
		~CellList() {};

		/** the saved lists currently using concurrent vectors
		 * due to writting conflicts when building the lists */
		ListDataVector particle_data_lists_;
	};

	/**
	 * @class MeshCellLinkedList
	 * @brief defining a mesh cell linked list for a body
	 * the meshes for all bodies share the same global coordiantes
	 */
	class MeshCellLinkedList : public Mesh
	{
	protected:
		/** cut_off radius */
		Real cutoff_radius_;
		/** cell spacing */
		Real cell_spacing_;
		/** number cells by dimension */
		Vecu number_of_cells_;

		/** clear the cell lists */
		void ClearCellLists();

		/** computing search range for building contact configuration */
		int ComputingSearchRage(int orign_refinement_level, 
			int target_refinement_level);
		/** choose a kernel for building up inter refinement level configuration */
		Kernel& ChoosingKernel(Kernel *orignal_kernel, Kernel *target_kernel);

	public:
		MeshCellLinkedList(Vecd lower_bound, Vecd upper_bound, 
			Real cell_spacing, size_t buffer_size = 0);
		virtual ~MeshCellLinkedList() {};

		/** access protected members */
		Real GetCutoffRadius() { return cutoff_radius_; };
		Real GetCellSpacing() { return cell_spacing_; };
		Vecu GetNumberOfCells() { return number_of_cells_; };

		/** number of lists in each cell,
		 * which is the the number of bodies
		 * all the lists numer of cell times number of bodies
		 * a list saved with the index of partilces within the cell 
		 */
		matrix_cell cell_linked_lists_;

		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		virtual void DeleteMeshDataMatrix() override;

		/** update cell lists */
		void UpdateCellLists(SPHBody &body);
		/** only update particle locations in cells */
		void UpdateParticleCellLocation(SPHBody &body);
		/** update by cell particle list for splitting algorithm */
		void UpdateByCellParticleLists(SPHBody &body);

		/** build reference inner configurtion */
		void BuildReferenceInnerConfiguration(SPHBody &body);
		/** build reference contact configuration */
		void BuildReferenceContactConfiguration(SPHBody &body);

		/** update inner configuration */
		void UpdateInnerConfiguration(SPHBody &body);
		/** update contact configuration */
		void UpdateContactConfiguration(SPHBody &body);
		/** update interaction configuration */
		void UpdateInteractionConfiguration(SPHBody &body,
			SPHBodyVector interacting_bodies);

		/** output mesh data for visuallization */
		virtual void WriteMeshToVtuFile(ofstream &output_file) override {};
		virtual void WriteMeshToPltFile(ofstream &output_file) override {};

		/** check if wether the two points are in order repsected to the coordinates */
		bool IsInSequence(Vecd pos_1, Vecd pos_2);
	};
}
