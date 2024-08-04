#ifndef UPDATE_CELL_LINKED_LIST_SYCL_HPP
#define UPDATE_CELL_LINKED_LIST_SYCL_HPP

#include "update_cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType, ParallelDevicePolicy>::exec(Real dt = 0.0)
{
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_HPP
