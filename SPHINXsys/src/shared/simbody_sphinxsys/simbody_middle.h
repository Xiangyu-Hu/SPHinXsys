/**
 * @file 	simbody_middle.h
 * @brief 	file to include Simbody headers and suppress their warnings
 * @author	Wen-Yang Chu
 */

#ifndef SIMBODY_MIDDLE_H
#define SIMBODY_MIDDLE_H

#ifdef __linux__ 
#pragma GCC system_header //for GCC/CLANG
#elif _WIN32
#pragma warning(push, 0) //for MSVC
#endif

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

#ifdef _WIN32 
#pragma warning(pop) //for MSVC
#endif
#endif //SIMBODY_MIDDLE_H