/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file complex_shape_mesh.h
* @brief x. 
* @details x 
*			x. 
* @author	Yijin Mao
*/

#ifndef COMPLEX_SHAPE_MESH_H
#define COMPLEX_SHAPE_MESH_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "SimTKcommon.h"
#include "SimTKmath.h"

#include <iostream>
#include <string>
#include <fstream>

/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH
{
	class ComplexShapeMesh
	{
	public:
		//constructor for load stl file from out side
		ComplexShapeMesh(std::string name) : name_(name) {};
		virtual ~ComplexShapeMesh() {};
		virtual BoundingBox findBounds()  = 0;
		virtual bool checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED = true) = 0;
		virtual bool checkNotFar(const Vec3d& input_pnt, Real threshold) = 0;
		virtual bool checkNearSurface(const Vec3d& input_pnt, Real threshold) = 0;
		/** Signed distance is negative for point within the complex shape. */
		virtual Real findSignedDistance(const Vec3d& input_pnt) = 0;
		/** Normal direction point toward outside of the complex shape. */
		virtual Vec3d findNormalDirection(const Vec3d& input_pnt) = 0;

		virtual Vec3d findClosestPoint(const Vec3d& input_pnt) = 0;

		std::string getName() { return name_; }

	protected:
		std::string name_;
	};
}

#endif //COMPLEX_SHAPE_MESH_H
