/**
* @file base_geometry.h
* @brief Shape is the base class for all geometries. 
* @details Several pure virtual functions 
* are defined here. (a) closet point on surface: to find the closet point on shape
* surface to a given point. (b) find the lower and upper bounds.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_data_package.h"
#include <string>
using namespace std;

namespace SPH
{
	/**
	 * @class ShapeBooleanOps
	 * @brief Boolian operation for generate complex shapes
	 * @details Note that, for 3D applications, only add and sub boolean operation have been defined right now
	 */
	enum class ShapeBooleanOps { add, sub, sym_diff, intersect };

	/**
	 * @class Shape
	 * @brief Base class for all geometries
	 */
	class Shape
	{
	public:
		Shape(string shape_name) : name_(shape_name) {};
		virtual ~Shape() {};

		string getName() { return name_; };
		virtual Vecd findClosestPoint(Vecd input_pnt) = 0;
		virtual void findBounds(Vecd& lower_bound, Vecd& upper_bound) = 0;

	protected:
		string name_;
	};
}
