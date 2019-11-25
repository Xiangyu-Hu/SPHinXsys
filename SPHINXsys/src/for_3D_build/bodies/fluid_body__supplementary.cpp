#include "fluid_body.h"
#include "mesh_cell_linked_list.h"

namespace SPH 
{
	//===============================================================//
	void FluidBodyPart::TagBodyPartCells()
	{
		MeshCellLinkedList *mesh_cell_linked_list
			= fluid_body_->mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->GetNumberOfCells();

		for (int i = 0; i < number_of_cells[0]; ++i)
			for (int j = 0; j < number_of_cells[1]; ++j)
				for (int k = 0; k < number_of_cells[2]; ++k)
				{
					bool is_contained = false;
					for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++m)
							for (int n = SMAX(k - 1, 0); n <= SMIN(k + 1, int(number_of_cells[2]) - 1); ++n)
							{
								Vecd cell_position = mesh_cell_linked_list
									->CellPositionFromIndexes(Vecu(l, m, n));
								if (body_part_region_.contain(cell_position))
									is_contained = true;
							}
					if (is_contained == true) body_part_cells_.push_back(Vecu(i, j, k));
				}
	}
	//===============================================================//
}
