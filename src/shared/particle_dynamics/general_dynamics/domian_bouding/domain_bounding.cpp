#include "domain_bounding.h"

namespace SPH
{
//=================================================================================================//
PeriodicConditionUsingCellLinkedList::
    PeriodicConditionUsingCellLinkedList(RealBody &real_body, PeriodicAlongAxis &periodic_box)
    : BasePeriodicCondition<execution::ParallelPolicy>(real_body, periodic_box),
      bounding_(bound_cells_data_, real_body, periodic_box),
      update_cell_linked_list_(bound_cells_data_, real_body, periodic_box) {}
//=================================================================================================//
PeriodicConditionUsingCellLinkedList::PeriodicCellLinkedList::
    PeriodicCellLinkedList(StdVec<CellLists> &bound_cells_data,
                           RealBody &real_body, PeriodicAlongAxis &periodic_box)
    : PeriodicBounding(bound_cells_data, real_body, periodic_box),
      cell_linked_list_(real_body.getCellLinkedList()) {}
//=================================================================================================//
void PeriodicConditionUsingCellLinkedList::
    PeriodicCellLinkedList::InsertListDataNearUpperBound(ListDataVector &cell_list_data, Real dt)
{
    for (size_t num = 0; num < cell_list_data.size(); ++num)
    {
        Vecd particle_position = std::get<1>(cell_list_data[num]);
        if (particle_position[axis_] < bounding_bounds_.upper_[axis_] &&
            particle_position[axis_] > (bounding_bounds_.upper_[axis_] - cut_off_radius_max_))
        {
            Vecd translated_position = particle_position - periodic_translation_;
            /** insert ghost particle to cell linked list */
            mutex_cell_list_entry_.lock();
            cell_linked_list_.InsertListDataEntry(cell_list_data[num].first, translated_position);
            mutex_cell_list_entry_.unlock();
        }
    }
}
//=================================================================================================//
void PeriodicConditionUsingCellLinkedList::
    PeriodicCellLinkedList::InsertListDataNearLowerBound(ListDataVector &cell_list_data, Real dt)
{
    for (size_t num = 0; num < cell_list_data.size(); ++num)
    {
        Vecd particle_position = std::get<1>(cell_list_data[num]);
        if (particle_position[axis_] > bounding_bounds_.lower_[axis_] &&
            particle_position[axis_] < (bounding_bounds_.lower_[axis_] + cut_off_radius_max_))
        {
            Vecd translated_position = particle_position + periodic_translation_;
            /** insert ghost particle to cell linked list */
            mutex_cell_list_entry_.lock();
            cell_linked_list_.InsertListDataEntry(cell_list_data[num].first, translated_position);
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
                 { InsertListDataNearLowerBound(*cell_ist, dt); });

    particle_for(execution::ParallelPolicy(), bound_cells_data_[1].second,
                 [&](ListDataVector *cell_ist)
                 { InsertListDataNearUpperBound(*cell_ist, dt); });
}
//=================================================================================================//
} // namespace SPH
