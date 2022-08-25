/**
 * @file 	scalar_functions_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "scalar_functions.h"

namespace SPH {
	//=================================================================================================//
	int SecondAxis(int first_axis) {
		return first_axis == 1 ? 0 : 1;
	}
}
