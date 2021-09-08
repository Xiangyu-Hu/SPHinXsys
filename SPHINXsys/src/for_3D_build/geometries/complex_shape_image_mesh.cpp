#include "complex_shape_image_mesh.h"
#include "image_mesh_shape.h"

namespace SPH
{
	//=================================================================================================//
	bool ComplexShapeImageMesh::checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED)
	{
		bool exist = false;
		bool inside = false;
		
		for (auto& each_shape : image_mesh_shapes_)
		{
			ImageMeshShape* sp = each_shape.first;
			ShapeBooleanOps operation_string = each_shape.second;

			switch (operation_string)
			{
			case ShapeBooleanOps::add:
			{
				inside = sp->checkContain(input_pnt);
				exist = exist || inside;
				break;
			}
			case ShapeBooleanOps::sub:
			{
				inside = sp->checkContain(input_pnt);
				exist = exist && (!inside);
				break;
			}
			default:
			{
				throw std::runtime_error("FAILURE: the boolean operation is not applicable!");
				break;
			}
			}
		}
		return exist;
	}
	//=================================================================================================//
	Vec3d ComplexShapeImageMesh::findClosestPoint(const Vec3d& input_pnt)
	{
		//a big positive number
		Real large_number(Infinity);
		Real dist_min = large_number;
		Vec3d pnt_closest(0);
		Vec3d pnt_found(0);

		for (auto& each_shape : image_mesh_shapes_)
		{
			ImageMeshShape* sp = each_shape.first;
			pnt_found  = sp->findClosestPoint(input_pnt);
			Real dist = (input_pnt - pnt_found).norm();

			if(dist <= dist_min)
			{
				dist_min = dist;
				pnt_closest = pnt_found;
			}
		}

		return pnt_closest;
	}
	//=================================================================================================//
	Real ComplexShapeImageMesh::findSignedDistance(const Vec3d& input_pnt)
	{
		return image_mesh_shapes_[0].first->findValueAtPoint(input_pnt);
	}
	//=================================================================================================//
	Vec3d ComplexShapeImageMesh::findNormalDirection(const Vec3d& input_pnt)
	{
		Vecd direction_to_surface = image_mesh_shapes_[0].first->findNormalAtPoint(input_pnt);
		bool is_contain = image_mesh_shapes_[0].first->checkContain(input_pnt);
		return is_contain ? direction_to_surface : -1.0 * direction_to_surface;
	}
	//=================================================================================================//
	void ComplexShapeImageMesh::addImageMeshShape(ImageMeshShape* image_mesh_shape, ShapeBooleanOps op)
	{
		std::pair<ImageMeshShape*, ShapeBooleanOps> shape_and_op(image_mesh_shape, op);
		image_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeImageMesh::addBrick(Vec3d halfsize, int resolution, Vec3d translation, Mat3d rotation, ShapeBooleanOps op)
	{
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		throw std::runtime_error("addBrick is not implemented");
	}
	//=================================================================================================//
	void ComplexShapeImageMesh::addSphere(Real radius, int resolution, Vec3d translation, Mat3d rotation, ShapeBooleanOps op)
	{
		Vec3d spacings(resolution,resolution,resolution);
		Vec3d center(translation);
		ImageMeshShape* image_mesh_shape = new ImageMeshShape(radius, spacings, center);
		std::pair<ImageMeshShape*, ShapeBooleanOps> shape_and_op(image_mesh_shape, op);
		image_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeImageMesh::addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, Mat3d rotation, ShapeBooleanOps op)
	{
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		throw std::runtime_error("addCylinder is not implemented");
	}
	//=================================================================================================//
	void ComplexShapeImageMesh::addComplexShapeImageMesh(ComplexShapeImageMesh* complex_shape_image_mesh, ShapeBooleanOps op)
	{
		switch (op)
		{
		case ShapeBooleanOps::add:
		{
			for (auto& shape_and_op : complex_shape_image_mesh->image_mesh_shapes_)
			{
				image_mesh_shapes_.push_back(shape_and_op);
			}
			break;
		}
		case ShapeBooleanOps::sub:
		{
			for (auto& shape_and_op : complex_shape_image_mesh->image_mesh_shapes_)
			{
				ImageMeshShape* sp = shape_and_op.first;
				ShapeBooleanOps operation_string 
					= shape_and_op.second == ShapeBooleanOps::add ? ShapeBooleanOps::sub : ShapeBooleanOps::add;
				std::pair<ImageMeshShape*, ShapeBooleanOps> substract_shape_and_op(sp, operation_string);

				image_mesh_shapes_.push_back(substract_shape_and_op);
			}
			break;
		}
		default:
		{
			throw std::runtime_error("unknown operation shape boolean operator");
			break;
		}
		}
	}
	//=================================================================================================//
	bool ComplexShapeImageMesh::checkNotFar(const Vec3d& input_pnt, Real threshold)
	{
		return  checkContain(input_pnt) || checkNearSurface(input_pnt , threshold) ? true : false;
	}
	//=================================================================================================//
	bool ComplexShapeImageMesh::checkNearSurface(const Vec3d& input_pnt, Real threshold)
	{
		return  getMaxAbsoluteElement(input_pnt - findClosestPoint(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	BoundingBox ComplexShapeImageMesh::findBounds()
	{
		//initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);

		for (size_t i = 0; i < image_mesh_shapes_.size(); i++)
		{
			BoundingBox shape_bounds = image_mesh_shapes_[i].first->findBounds();
			for (int j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], shape_bounds.first[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_bounds.second[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
}