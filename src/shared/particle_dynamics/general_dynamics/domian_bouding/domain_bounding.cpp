#include "domain_bounding.h"

namespace SPH
{
//=================================================================================================//
PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::
    PeriodicCellLinkedList(StdVec<CellLists> &bound_cells_data,
                           RealBody &real_body, BoundingBox bounding_bounds, int axis)
    : PeriodicBounding(bound_cells_data, real_body, bounding_bounds, axis),
      bound_cells_data_(bound_cells_data),
      cell_linked_list_(real_body.getCellLinkedList()),
      cut_off_radius_max_(real_body.sph_adaptation_->getKernel()->CutOffRadius()) {}
//=================================================================================================//
void PeriodicConditionUsingCellLinkedList::
    PeriodicCellLinkedList::checkUpperBound(ListDataVector &cell_list_data, Real dt)
{
    for (size_t num = 0; num < cell_list_data.size(); ++num)
    {
        Vecd particle_position = std::get<1>(cell_list_data[num]);
        if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
            particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
        {
            Vecd translated_position = particle_position - periodic_translation_;
            /** insert ghost particle to cell linked list */
            mutex_cell_list_entry_.lock();
            cell_linked_list_.InsertListDataEntry(std::get<0>(cell_list_data[num]),
                                                  translated_position, std::get<2>(cell_list_data[num]));
            mutex_cell_list_entry_.unlock();
        }
    }
}
//=================================================================================================//
void PeriodicConditionUsingCellLinkedList::
    PeriodicCellLinkedList::checkLowerBound(ListDataVector &cell_list_data, Real dt)
{
    for (size_t num = 0; num < cell_list_data.size(); ++num)
    {
        Vecd particle_position = std::get<1>(cell_list_data[num]);
        if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
            particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
        {
            Vecd translated_position = particle_position + periodic_translation_;
            /** insert ghost particle to cell linked list */
            mutex_cell_list_entry_.lock();
            cell_linked_list_.InsertListDataEntry(std::get<0>(cell_list_data[num]),
                                                  translated_position, std::get<2>(cell_list_data[num]));
            mutex_cell_list_entry_.unlock();
        }
    }
}
//=================================================================================================//
void PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::exec(Real dt)
{
    setupDynamics(dt);

    particle_for(execution::ParallelPolicy(), bound_cells_data_[0].second,
                 [&](ListDataVector *cell_ist)
                 { checkLowerBound(*cell_ist, dt); });

    particle_for(execution::ParallelPolicy(), bound_cells_data_[1].second,
                 [&](ListDataVector *cell_ist)
                 { checkUpperBound(*cell_ist, dt); });
}
//=================================================================================================//
} // namespace SPH
