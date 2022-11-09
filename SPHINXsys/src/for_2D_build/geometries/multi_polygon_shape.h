/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
* @file multi_polygon_shape.h
* @brief Here, we define the 2D geometric algorithms. they are based on the boost library. 
* @details The idea is to define complex geometry based on shapes, usually
* multi-polygon using boost library. we propose only very simple combination
* that the region is composed of shapes without intersection.
* That is, the shapes are those contain each other or without overlap.
* This strict requirement suggests that complex shapes should be finished
* already in modeling using related binary operations before it is included.
* @author	Luhui Han, Chi Zhang and Xiangyu Hu
*/

#ifndef MULTI_POLYGON_SHAPE_H
#define MULTI_POLYGON_SHAPE_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

//boost library
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

#include "base_data_package.h"
#include "base_geometry.h"

#include <iostream>
#include <string>
#include <fstream>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

using namespace boost::geometry;

namespace SPH
{

	/**
	 * @brief preclaimed classes.
	 */
	class Kernel;

	typedef model::polygon<model::d2::point_xy<Real>> boost_poly;
	typedef model::multi_polygon<boost_poly> boost_multi_poly;

	/**
	 * @class MultiPolygon
	 * @brief used to define a closed region
	 */
	class MultiPolygon
	{
	public:
		MultiPolygon(){};
		explicit MultiPolygon(const std::vector<Vecd> &points);
		explicit MultiPolygon(const Vec2d &center, Real radius, int resolution);
		boost_multi_poly &getBoostMultiPoly() { return multi_poly_; };

		BoundingBox findBounds();
		bool checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED = true);
		Vec2d findClosestPoint(const Vec2d &probe_point);

		void addAMultiPolygon(MultiPolygon &multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly &boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(const std::vector<Vecd> &points, ShapeBooleanOps op);
		void addABox(Transform2d transform2d, const Vec2d &halfsize, ShapeBooleanOps op);
		void addACircle(const Vec2d &center, Real radius, int resolution, ShapeBooleanOps op);
		void addAPolygonFromFile(std::string file_path_name, ShapeBooleanOps op, Vec2d translation = Vecd(0), Real scale_factor = 1.0);

	protected:
		boost_multi_poly multi_poly_;
		boost_multi_poly MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
												  boost_multi_poly multi_poly_op,
												  ShapeBooleanOps boolean_op);
	};

	/**
	 * @class MultiPolygonShape
	 * @brief A shape whose geometry is defined by a multi polygon.
	 */
	class MultiPolygonShape : public Shape
	{

	public:
		/** Default constructor. */
		explicit MultiPolygonShape(const std::string &shape_name) : Shape(shape_name){};
		explicit MultiPolygonShape(const MultiPolygon &multi_polygon, const std::string &shape_name = "MultiPolygonShape")
			: Shape(shape_name), multi_polygon_(multi_polygon){};
		virtual ~MultiPolygonShape(){};

		virtual bool isValid() override;
		virtual bool checkContain(const Vec2d &probe_point, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec2d findClosestPoint(const Vec2d &probe_point) override;

	protected:
		MultiPolygon multi_polygon_;

		virtual BoundingBox findBounds() override;
	};
}

#endif //MULTI_POLYGON_SHAPE_H