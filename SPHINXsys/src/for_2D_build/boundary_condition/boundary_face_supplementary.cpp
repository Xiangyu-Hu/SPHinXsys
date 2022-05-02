
#include "boundary_face.h"

namespace SPH
{
    void BodyRegionByCellsWithFace::tagBodyDomainBoundingCells(CellLists& bound_cells)
    {
        bound_cells.clear();
        BaseCellLinkedList* base_cell_linked_list = real_body_.cell_linked_list_;
        CellLinkedList* cell_linked_list = dynamic_cast<CellLinkedList*>(base_cell_linked_list);
        if (cell_linked_list == nullptr)
        {
            std::cout << "\n FAILURE: convert cell linked list failed " << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        Vecu lower = cell_linked_list->CellIndexFromPosition(region_bounds.first);
        Vecu upper = cell_linked_list->CellIndexFromPosition(region_bounds.second);
        Vecu cell_num = cell_linked_list->NumberOfCells();

        for (size_t j = SMAX(lower[1] - 1, size_t(0)); j < SMIN(upper[1] + 1, cell_num[1] - 1); ++j)
            for (size_t i = SMAX(lower[0] - 1, size_t(0)); i < SMIN(upper[0] + 1, cell_num[0] - 1); ++i)
            {
                bound_cells.emplace_back(&cell_linked_list->getCellLists()[i][j]);
            }
    }

} // namespace SPH