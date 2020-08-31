/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
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
