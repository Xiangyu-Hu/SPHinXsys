#include "scalar_functions.h"

namespace SPH {
	//=================================================================================================//
	int SecondAxis(int axis_direction) {
		return axis_direction == 2 ? 0 : axis_direction + 1;
	}
}
