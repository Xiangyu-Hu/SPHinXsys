/**
 * @file 	base_geometry.cpp
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#include "base_geometry.h"
namespace SPH
{
	//=================================================================================================//
	bool Shape::checkNotFar(const Vecd &input_pnt, Real threshold)
	{
		return checkContain(input_pnt) || checkNearSurface(input_pnt, threshold) ? true : false;
	}
	//=================================================================================================//
	bool Shape::checkNearSurface(const Vecd &input_pnt, Real threshold)
	{
		return getMaxAbsoluteElement(input_pnt - findClosestPoint(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	Real Shape::findSignedDistance(const Vecd &input_pnt)
	{
		Real distance_to_surface = (input_pnt - findClosestPoint(input_pnt)).norm();
		return checkContain(input_pnt) ? -distance_to_surface : distance_to_surface;
	}
	//=================================================================================================//
	Vecd Shape::findNormalDirection(const Vecd &input_pnt)
	{
		bool is_contain = checkContain(input_pnt);
		Vecd displacement_to_surface = findClosestPoint(input_pnt) - input_pnt;
		while (displacement_to_surface.norm() < Eps)
		{
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
	BoundingBox BinaryShapes::findBounds()
	{
		//initial reference values
		Vecd lower_bound = Vecd(Infinity);
		Vecd upper_bound = Vecd(-Infinity);

		for (auto &shape_and_op : shapes_and_ops_)
		{
			BoundingBox shape_bounds = shape_and_op.first->findBounds();
			for (int j = 0; j != Dimensions; ++j)
			{
				lower_bound[j] = SMIN(lower_bound[j], shape_bounds.first[j]);
				upper_bound[j] = SMAX(upper_bound[j], shape_bounds.second[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
	bool BinaryShapes::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
	{
		bool exist = false;
		bool inside = false;

		for (auto &shape_and_op : shapes_and_ops_)
		{
			Shape *geometry = shape_and_op.first;
			ShapeBooleanOps operation_string = shape_and_op.second;
			switch (operation_string)
			{
			case ShapeBooleanOps::add:
			{
				inside = geometry->checkContain(pnt);
				exist = exist || inside;
				break;
			}
			case ShapeBooleanOps::sub:
			{
				inside = geometry->checkContain(pnt);
				exist = exist && (!inside);
				break;
			}
			default:
			{
				std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				throw;
			}
			}
		}
		return exist;
	}
	//=================================================================================================//
	Vecd BinaryShapes::findClosestPoint(const Vecd &input_pnt)
	{
		//a big positive number
		Real large_number(Infinity);
		Real dist_min = large_number;
		Vecd pnt_closest(0);
		Vecd pnt_found(0);

		for (auto &shape_and_op : shapes_and_ops_)
		{
			Shape *geometry = shape_and_op.first;
			pnt_found = geometry->findClosestPoint(input_pnt);
			Real dist = (input_pnt - pnt_found).norm();

			if (dist <= dist_min)
			{
				dist_min = dist;
				pnt_closest = pnt_found;
			}
		}

		return pnt_closest;
	}
	//=================================================================================================//
	ShapeAndOp *BinaryShapes::getShapeAndOpByName(const std::string &shape_name)
	{
		for (auto &shape_and_op : shapes_and_ops_)
		{
			Shape *shape = shape_and_op.first;
			if (shape->getName() == shape_name)
				return &shape_and_op;
		}
		std::cout << "\n FAILURE: the shape " << shape_name << " has not been created!" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;

		return nullptr;
	}
	//=================================================================================================//
	Shape *BinaryShapes::getShapeByName(const std::string &shape_name)
	{
		return getShapeAndOpByName(shape_name)->first;
	}
	//=================================================================================================//
}