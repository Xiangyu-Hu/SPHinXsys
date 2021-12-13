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
* @file image_mesh_shape.h
* @brief x 
* @details x
*			x
* @author	Yijin Mao
*/

#ifndef IMAGE_MESH_SHAPE_H
#define IMAGE_MESH_SHAPE_H

#ifndef __EMSCRIPTEN__

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"
#include "simbody_middle.h"
#include "geometry.h"
#include "image_mhd.h"
#include "base_data_type.h"

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
//- imageMesh is input 3d image
	class ImageMeshShape
	{
	public:
		//constructor for load mhd/raw file from out side
		ImageMeshShape(std::string file_path_name);
		ImageMeshShape(Real radius, Vec3d spacings, Vec3d center);
		virtual ~ImageMeshShape();

		bool checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED = true);
		Vec3d findClosestPoint(const Vec3d& input_pnt);
		BoundingBox findBounds();

		Real findValueAtPoint(const Vec3d& input_pnt);
		Vec3d findNormalAtPoint(const Vec3d & input_pnt);

	protected:
		//- distance map has to be float type image
		Vec3d translation_;
		Mat3d rotation_;
		std::unique_ptr<ImageMHD<float, 3>> image_;

		Real max_distance_;
		Real min_distance_;

	};
}

#endif // __EMSCRIPTEN__

#endif //IMAGE_MESH_SHAPE_H
