/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
* @file image_shape.h
* @brief x 
* @details x 
*			x 
* @author	Yijin Mao
*/

#ifndef IMAGE_SHAPE_3D_H
#define IMAGE_SHAPE_3D_H

#ifndef __EMSCRIPTEN__

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_geometry.h"
#include "image_mhd.h"

#include <iostream>
#include <string>
#include <fstream>

/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace SPH
{
	class ImageShape : public Shape
	{
	public:
		explicit ImageShape(const std::string &shape_name)
			: Shape(shape_name), translation_(0.0), rotation_(1.0),
			  max_distance_(-INFINITY), min_distance_(INFINITY){};

		virtual bool checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec3d findClosestPoint(const Vec3d &probe_point) override;

	protected:
		//- distance map has to be float type image
		Vec3d translation_;
		Mat3d rotation_;
		std::unique_ptr<ImageMHD<float, 3>> image_;
		Real max_distance_;
		Real min_distance_;
	
		virtual BoundingBox findBounds() override;
	};

	class ImageShapeFromFile : public ImageShape
	{
	public:
		//constructor for load mhd/raw file from out side
		explicit ImageShapeFromFile(const std::string &file_path_name,
									const std::string &shape_name = "ImageShapeFromFile");
		virtual ~ImageShapeFromFile(){};
	};

	class ImageShapeSphere : public ImageShape
	{
	public:
		//constructor for load mhd/raw file from out side
		ImageShapeSphere(Real radius, Vec3d spacings, Vec3d center,
						 const std::string &shape_name = "ImageShapeSphere");
		virtual ~ImageShapeSphere(){};
	};
}

#endif //__EMSCRIPTEN__

#endif //IMAGE_SHAPE_3D_H
