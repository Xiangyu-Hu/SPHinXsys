#ifndef __EMSCRIPTEN__

#include "image_shape.h"

namespace SPH
{
//=================================================================================================//
bool ImageShape::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    Real value = image_->findValueAtPoint(probe_point);
    return BOUNDARY_INCLUDED ? value <= 0.0 : value < 0.0;
}
//=================================================================================================//
Vecd ImageShape::findClosestPoint(const Vecd &probe_point)
{
    return image_->findClosestPoint(probe_point);
}
//=================================================================================================//
BoundingBoxd ImageShape::findBounds()
{
    return image_->findBounds();
}
//=================================================================================================//
ImageShapeFromFile::
    ImageShapeFromFile(const std::string &file_path_name, const std::string &shape_name)
    : ImageShape(shape_name)
{
    image_.reset(new ImageMHD<float, 3>(file_path_name));
}
//=================================================================================================//
ImageShapeSphere::
    ImageShapeSphere(Real radius, Vecd spacings, Vecd center, const std::string &shape_name)
    : ImageShape(shape_name)
{
    Real extend = 1.5;
    int length = int(std::ceil(2.0 * extend * radius));
    Arrayi NxNyNz(length, length, length);
    image_.reset(new ImageMHD<float, 3>(radius, NxNyNz, spacings));
}
//=================================================================================================//
} // namespace SPH

#endif // __EMSCRIPTEN__
