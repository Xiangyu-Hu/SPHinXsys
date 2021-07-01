/**
 * @file 	scalar_functions.cpp
 * @author	Xiangyu Hu
 * @version	0.1
 */

#include "scalar_functions.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	int ThirdAxis(int axis_direction) {
		return SecondAxis(SecondAxis(axis_direction));
	}

}
