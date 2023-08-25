#include "scalar_functions.h"

namespace SPH
{
//=================================================================================================//
int SecondAxis(int first_axis)
{
    return first_axis == 1 ? 0 : 1;
}
} // namespace SPH
