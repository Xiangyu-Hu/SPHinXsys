#include "image_mesh_shape.h"

namespace SPH
{
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(std::string file_path_name) :
	translation_(0.0),
	rotation_(1.0),
	image_(nullptr),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{
		if (image_ == nullptr)
			image_.reset(new ImageMHD<float,3>(file_path_name));
	}
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(Real radius, Vec3d spacings, Vec3d center) :
	translation_(0.0),
	rotation_(1.0),
	image_(nullptr),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{
		double extend = 1.5;
		int length = int(std::ceil(2.0*extend*radius));
		Vec3i NxNyNz(length, length, length);
		if(image_ == nullptr)
			image_.reset(new ImageMHD<float,3>(radius, NxNyNz,spacings));
	}
	//=================================================================================================//
	ImageMeshShape::~ImageMeshShape()
	{
	}
	//=================================================================================================//
	bool ImageMeshShape::checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED)
	{
        Real value = findValueAtPoint(input_pnt);
        if (BOUNDARY_INCLUDED == true)
        {
            if (value > 0.0) return false;
            else return true;
        }
        else
        {
            if (value >= 0.0) return false;
            else return true;
        }
	}
	//=================================================================================================//
	Vec3d ImageMeshShape::findClosestPoint(const Vec3d& input_pnt)
	{
		return image_->findClosestPoint(input_pnt);
	}
	//=================================================================================================//
	BoundingBox ImageMeshShape::findBounds()
	{
        return image_->findBounds();
	}
	//=================================================================================================//
	Real ImageMeshShape::findValueAtPoint(const Vec3d& input_pnt)
	{
		return image_->findValueAtPoint(input_pnt);
	}
	//=================================================================================================//
	Vec3d ImageMeshShape::findNormalAtPoint(const Vec3d & input_pnt)
	{
		return image_->findNormalAtPoint(input_pnt);
	}
	//=================================================================================================//

}
