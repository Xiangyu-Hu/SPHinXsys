#include "complex_shape_triangle_mesh.h"
#include "geometry.h"

namespace SPH
{
	//=================================================================================================//
	bool ComplexShapeTriangleMesh::checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED)
	{
		bool exist = false;
		bool inside = false;
		
		for (auto& each_shape : triangle_mesh_shapes_)
		{
			TriangleMeshShape* sp = each_shape.first;
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
				std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
				break;
			}
			}
		}
		return exist;
	}
	//=================================================================================================//
	Vec3d ComplexShapeTriangleMesh::findClosestPoint(const Vec3d& input_pnt)
	{
		//a big positive number
		Real large_number(Infinity);
		Real dist_min = large_number;
		Vec3d pnt_closest(0);
		Vec3d pnt_found(0);

		for (auto& each_shape : triangle_mesh_shapes_)
		{
			TriangleMeshShape* sp = each_shape.first;
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
	Real ComplexShapeTriangleMesh::findSignedDistance(const Vec3d& input_pnt)
	{
		Real distance_to_surface = (input_pnt - findClosestPoint(input_pnt)).norm();
		return checkContain(input_pnt) ? -distance_to_surface : distance_to_surface;
	}
	//=================================================================================================//
	Vec3d ComplexShapeTriangleMesh::findNormalDirection(const Vec3d& input_pnt)
	{
		bool is_contain = checkContain(input_pnt);
		Vecd displacement_to_surface = findClosestPoint(input_pnt) - input_pnt;
		while (displacement_to_surface.norm() < Eps) {
			Vecd jittered = input_pnt; //jittering
			for (int l = 0; l != input_pnt.size(); ++l)
				jittered[l] = input_pnt[l] + (((Real)rand() / (RAND_MAX)) - 0.5) * 100.0 * Eps;
			if (checkContain(jittered) == is_contain)
				displacement_to_surface = findClosestPoint(jittered) - jittered;
		}
		Vecd direction_to_surface = displacement_to_surface.normalize();
		return is_contain ? direction_to_surface : -1.0 * direction_to_surface;
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addTriangleMeshShape(TriangleMeshShape* triangle_mesh_shape, ShapeBooleanOps op)
	{
		std::pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addBrick(Vec3d halfsize, int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(halfsize, resolution, translation);
		std::pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addSphere(Real radius, int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(radius, resolution, translation);
		std::pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, ShapeBooleanOps op)
	{
		/** Here SimTK::UnitVec3 give the direction of the cylinder, viz. SimTK::UnitVec3(0,0,1) create a cylinder in z-axis.*/
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(axis, radius, halflength, resolution, translation);
		std::pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addFormSTLFile(std::string file_path_name, Vec3d translation, Real scale_factor, ShapeBooleanOps op)
	{
		TriangleMeshShape* triangle_mesh_shape = new TriangleMeshShape(file_path_name, translation, scale_factor);
		std::pair<TriangleMeshShape*, ShapeBooleanOps> shape_and_op(triangle_mesh_shape, op);
		triangle_mesh_shapes_.push_back(shape_and_op);
	}
	//=================================================================================================//
	void ComplexShapeTriangleMesh::addComplexShapeTriangleMesh(ComplexShapeTriangleMesh* complex_shape_triangle_mesh, ShapeBooleanOps op)
	{
		switch (op)
		{
		case ShapeBooleanOps::add:
		{
			for (auto& shape_and_op : complex_shape_triangle_mesh->triangle_mesh_shapes_)
			{
				triangle_mesh_shapes_.push_back(shape_and_op);
			}
			break;
		}
		case ShapeBooleanOps::sub:
		{
			for (auto& shape_and_op : complex_shape_triangle_mesh->triangle_mesh_shapes_)
			{
				TriangleMeshShape* sp = shape_and_op.first;
				ShapeBooleanOps operation_string 
					= shape_and_op.second == ShapeBooleanOps::add ? ShapeBooleanOps::sub : ShapeBooleanOps::add;
				std::pair<TriangleMeshShape*, ShapeBooleanOps> substract_shape_and_op(sp, operation_string);

				triangle_mesh_shapes_.push_back(substract_shape_and_op);
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
	bool ComplexShapeTriangleMesh::checkNotFar(const Vec3d& input_pnt, Real threshold)
	{
		return  checkContain(input_pnt) || checkNearSurface(input_pnt , threshold) ? true : false;
	}
	//=================================================================================================//
	bool ComplexShapeTriangleMesh::checkNearSurface(const Vec3d& input_pnt, Real threshold)
	{
		return  getMaxAbsoluteElement(input_pnt - findClosestPoint(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	BoundingBox ComplexShapeTriangleMesh::findBounds()
	{
		//initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);

		for (size_t i = 0; i < triangle_mesh_shapes_.size(); i++)
		{
			BoundingBox shape_bounds = triangle_mesh_shapes_[i].first->findBounds();
			for (int j = 0; j != 3; ++j) {
				lower_bound[j] = SMIN(lower_bound[j], shape_bounds.first[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_bounds.second[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
}