/**
 * @file 	scalar_functions_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "scalar_functions.h"

namespace SPH {
	//=================================================================================================//
	int SeondAxis(int axis_direction) {
		return axis_direction == 1 ? 0 : 1;
	}
}
