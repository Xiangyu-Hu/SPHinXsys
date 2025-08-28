/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	image_shape.h
 * @brief 	Geometry processing with image shape.
 * @author	Yijin Mao, Chi Zhang and Xiangyu Hu
 */

#ifndef IMAGE_SHAPE_3D_H
#define IMAGE_SHAPE_3D_H

#ifndef __EMSCRIPTEN__

#include "base_geometry.h"
#include "image_mhd.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

namespace SPH
{
class ImageShape : public Shape
{
  public:
    explicit ImageShape(const std::string &shape_name)
        : Shape(shape_name), translation_(Vecd::Zero()), rotation_(Matd::Identity()),
          max_distance_(-INFINITY), min_distance_(INFINITY){};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;

  protected:
    Vecd translation_;
    Matd rotation_;
    std::unique_ptr<ImageMHD<float, 3>> image_;
    Real max_distance_;
    Real min_distance_;

    virtual BoundingBox findBounds() override;
};

class ImageShapeFromFile : public ImageShape
{
  public:
    explicit ImageShapeFromFile(const std::string &file_path_name,
                                const std::string &shape_name = "ImageShapeFromFile");
    virtual ~ImageShapeFromFile(){};
};

class ImageShapeSphere : public ImageShape
{
  public:
    ImageShapeSphere(Real radius, Vecd spacings, Vecd center,
                     const std::string &shape_name = "ImageShapeSphere");
    virtual ~ImageShapeSphere(){};
};
} // namespace SPH

#endif //__EMSCRIPTEN__

#endif // IMAGE_SHAPE_3D_H
