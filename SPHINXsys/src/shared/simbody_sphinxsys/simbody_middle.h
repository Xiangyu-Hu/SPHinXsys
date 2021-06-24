/**
 * @file 	simbody_middle.h
 * @brief 	file to include Simbody headers and suppress their warnings
 * @author	Wen-Yang Chu
 */

#ifndef SIMBODY_MIDDLE_H
#define SIMBODY_MIDDLE_H
#pragma GCC system_header //for GCC/CLANG
#pragma warning(push, 0) //for MSVC

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

#pragma warning(pop) //for MSVC
#endif //SIMBODY_MIDDLE_H