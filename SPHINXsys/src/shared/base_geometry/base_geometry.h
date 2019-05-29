/**
* @file shape.h
* @brief Shape is the base class for all geometies. It is a dommain closed by surfaces. 
* @details Several pure virtual functions 
* are defined here. (a) contain: the identify wether a point is conatined in
* the shape. (b) closet point on surface: to find the closet point on shape
* surafce to a given point. (c) find the lower and upper bounds (to be done).
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
	 * @class Shape
	 * @brief Base class for all geometies
	 */
	class Shape
	{
	protected:
		string name;

	public:
		Shape(string shapeName) : name(shapeName) {};
		virtual ~Shape() {};

		string getname() { return name; };

		virtual bool contain(Vecd pnt, bool BOUNDARY_INCLUDED = true) = 0;
		virtual Vecd closestpointonface(Vecd input_pnt) = 0;
		virtual void shapebound(Vecd &lower_bound, Vecd &upper_bound) = 0;
	};
}
