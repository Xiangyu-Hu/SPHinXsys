#include "scalar_functions.h"

namespace SPH
{
//=================================================================================================//
int NextAxis(int first_axis)
{
    return first_axis == 2 ? 0 : first_axis + 1;
}
} // namespace SPH
