/**
 * @file 	base_body_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_body.h"
#include "base_particles.h"
#include "mesh_cell_linked_list.h"
#include "level_set.h"

namespace SPH 
{
	//=================================================================================================//
	Real SPHBody::computeReferenceNumberDensity()
	{
		Real sigma(0);
		Real cutoff_radius = kernel_->GetCutOffRadius();
		Real particle_spacing = particle_spacing_;
		int search_range = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_range; j <= search_range; ++j)
			for (int i = -search_range; i <= search_range; ++i)
			{
				Point particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				if (particle_location.norm() < cutoff_radius)
					sigma +=  kernel_->W(particle_location);
			}
		return sigma;
	}
	//=================================================================================================//
	void BodyPartByCell::tagBodyPart()
	{
		BaseMeshCellLinkedList *mesh_cell_linked_list
			= body_->mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->NumberOfCells();

		for (int i = 0; i < (int)number_of_cells[0]; ++i)
			for (int j = 0; j < (int)number_of_cells[1]; ++j) {
				bool is_contained = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++l)
					{
						Vecd cell_position = mesh_cell_linked_list
							->CellPositionFromIndexes(Vecu(k, l));
						if (body_part_shape_->checkNotFar(cell_position, mesh_cell_linked_list->GridSpacing()))
						{
							if (body_part_shape_->checkContain(cell_position))
								is_contained = true;
						}
					}
				if (is_contained == true) 
					body_part_cells_.push_back(mesh_cell_linked_list->CellListFromIndex(Vecu(i, j)));
			}
	}
	//=================================================================================================//
	void NearBodySurface::tagBodyPart()
	{
		BaseMeshCellLinkedList* mesh_cell_linked_list
			= body_->mesh_cell_linked_list_;
		Vecu number_of_cells = mesh_cell_linked_list->NumberOfCells();

		for (int i = 0; i < (int)number_of_cells[0]; ++i)
			for (int j = 0; j < (int)number_of_cells[1]; ++j) {

				bool is_near = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++l)
					{
						Vecd cell_position = mesh_cell_linked_list->CellPositionFromIndexes(Vecu(k, l));
						if (body_->body_shape_->checkNotFar(cell_position, mesh_cell_linked_list->GridSpacing()))
						{
							Real phii = body_->body_shape_->findSignedDistance(cell_position);
							if (fabs(phii) <= mesh_cell_linked_list->GridSpacing()) is_near = true;
						}
					}
				if (is_near == true)
					body_part_cells_.push_back(mesh_cell_linked_list->CellListFromIndex(Vecu(i, j)));
			}
	}
	//=================================================================================================//
}
