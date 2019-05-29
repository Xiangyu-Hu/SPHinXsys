#include "weakly_compressible_fluid_body.h"
#include "weakly_compressible_fluid.h"
#include "mesh_cell_linked_list.h"
#include "sph_system.h"


namespace SPH 
{
	//===============================================================//
	void FluidBodyPart::TagBodyPartCells()
	{
		MeshCellLinkedList *mesh_cell_linked_list
			= fluid_body_->mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->GetNumberOfCells();

		for (int i = 0; i < number_of_cells[0]; ++i)
			for (int j = 0; j < number_of_cells[1]; ++j) {
				bool is_contained = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++l)
					{
						Vecd cell_position = mesh_cell_linked_list
							->CellPositionFromIndexes(Vecu(k, l));
						if (fluid_body_part_region_.contain(cell_position))
							is_contained = true;
					}
				if (is_contained == true) body_part_cells_.push_back(Vecu(i, j));
			}
	}
	//===============================================================//
}
