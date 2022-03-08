#ifndef __EMSCRIPTEN__

#include "image_shape.h"

namespace SPH
{
	//=================================================================================================//
	bool ImageShape::checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED)
	{
		Real value = findSignedDistance(input_pnt);
		if (BOUNDARY_INCLUDED == true)
		{
			if (value > 0.0)
				return false;
			else
				return true;
		}
		else
		{
			if (value >= 0.0)
				return false;
			else
				return true;
		}
	}
	//=================================================================================================//
	Vec3d ImageShape::findClosestPoint(const Vec3d &input_pnt)
	{
		return image_->findClosestPoint(input_pnt);
	}
	//=================================================================================================//
	BoundingBox ImageShape::findBounds()
	{
		return image_->findBounds();
	}
	//=================================================================================================//
	Real ImageShape::findSignedDistance(const Vec3d &input_pnt)
	{
		return image_->findValueAtPoint(input_pnt);
	}
	//=================================================================================================//
	Vec3d ImageShape::findNormalDirection(const Vec3d &input_pnt)
	{
		return image_->findNormalAtPoint(input_pnt);
	}
	//=================================================================================================//
	bool ImageShape::checkNotFar(const Vec3d &input_pnt, Real threshold)
	{
		return checkContain(input_pnt) || checkNearSurface(input_pnt, threshold) ? true : false;
	}
	//=================================================================================================//
	bool ImageShape::checkNearSurface(const Vec3d &input_pnt, Real threshold)
	{
		return getMaxAbsoluteElement(input_pnt - findClosestPoint(input_pnt)) < threshold ? true : false;
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
		ImageShapeSphere(Real radius, Vec3d spacings, Vec3d center, const std::string &shape_name)
		: ImageShape(shape_name)
	{
		double extend = 1.5;
		int length = int(std::ceil(2.0 * extend * radius));
		Vec3i NxNyNz(length, length, length);
		image_.reset(new ImageMHD<float, 3>(radius, NxNyNz, spacings));
	}
	//=================================================================================================//
}

#endif // __EMSCRIPTEN__
