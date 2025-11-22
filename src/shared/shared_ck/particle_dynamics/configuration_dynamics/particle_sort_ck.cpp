#include "particle_sort_ck.hpp"

namespace SPH
{
//=================================================================================================//
UpdateSortableVariables::UpdateSortableVariables(UnsignedInt data_size)
    : initialize_temp_variables_()
{
    initialize_temp_variables_(temp_variables_, data_size);
}
//=================================================================================================//
} // namespace SPH