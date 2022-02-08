/**
* @file multi_polygon_shape.h
* @brief Here, we define the 2D geometric algortihms. they are based on the boost library. 
* @details The idea is to define complex geometry based on shapes, usually
* multi-polygon using boost library. we propose only very simple combinaton
* that the region is composed of shapes without intersection.
* That is, the shapes are those contain each other or without overlap.
* This strict requirement suggests that complex shapes should be finished
* already in modeling using related binary operations before it is included.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
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
		boost_multi_poly &getBoostMultiPoly() { return multi_poly_; };

		BoundingBox findBounds();
		bool checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED = true);
		Vec2d findClosestPoint(const Vec2d &input_pnt);

		void addAMultiPolygon(MultiPolygon &multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly &boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(const std::vector<Vecd> &points, ShapeBooleanOps op);
		void addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op);
		void addAPolygonFromFile(std::string file_path_name, ShapeBooleanOps op, Vec2d translation = Vecd(0), Real scale_factor = 1.0);

	protected:
		boost_multi_poly multi_poly_;
		boost_multi_poly MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
												  boost_multi_poly multi_poly_op,
												  ShapeBooleanOps boolean_op);
	};

	/**
	 * @class MultiPolygonShape
	 * @brief A shape whose geometry is defined by a multipolygen.
	 */
	class MultiPolygonShape : public Shape
	{

	public:
		/** Default constructor. */
		explicit MultiPolygonShape(const MultiPolygon &multi_polygon, const std::string &shape_name = "MultiPolygonShape")
			: Shape(shape_name), multi_polygon_(multi_polygon){};
		virtual ~MultiPolygonShape(){};

		virtual BoundingBox findBounds() override;
		virtual bool checkContain(const Vec2d &input_pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec2d findClosestPoint(const Vec2d &input_pnt) override;
		virtual bool checkNotFar(const Vec2d &input_pnt, Real threshold) override;
		virtual bool checkNearSurface(const Vec2d &input_pnt, Real threshold) override;
		virtual Real findSignedDistance(const Vec2d &input_pnt) override;
		virtual Vec2d findNormalDirection(const Vec2d &input_pnt) override;

	protected:
		MultiPolygon multi_polygon_;
	};
}

#endif //MULTI_POLYGON_SHAPE_H