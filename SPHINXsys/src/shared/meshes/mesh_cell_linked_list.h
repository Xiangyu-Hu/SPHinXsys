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

#include "sph_data_conainers.h"
#include "base_mesh.h"

namespace SPH {

	class SPHSystem;
	class SPHBody;
	class Kernel;

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
	 * @brief Defining a mesh cell linked list for a body.
	 * The meshes for all bodies share the same global coordiantes.
	 */
	class MeshCellLinkedList : public Mesh
	{
	protected:
		/** cut_off radius */
		Real cutoff_radius_;
		/** cell spacing */
		Real cell_spacing_;
		/** number of cells by dimension */
		Vecu number_of_cells_;

		/** clear the cell lists */
		void ClearCellLists();
		/** clear by cell particle lists for a body*/
		void ClearByCellParticleLists(SPHBody &body);

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
		Real GetCellSpacing() { return cell_spacing_; };
		Vecu GetNumberOfCells() { return number_of_cells_; };

		/** The array for of mesh cells, i.e. mesh data.
		 * Within each cell, a list is saved with the index of partilces.*/
		matrix_cell cell_linked_lists_;

		/** allcate memories for mesh data */
		virtual void AllocateMeshDataMatrix() override;
		/** delete memories for mesh data */
		virtual void DeleteMeshDataMatrix() override;

		/** update the cell lists */
		void UpdateCellLists(SPHBody &body);
		/** Only update particles' cell locations. 
		 *  This is for fiticious body in which no cell linked list 
		 *  is required because no inner particle configuration to be constructed.
		 */
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
