/**
* @file geometry.h
* @brief Here, we define the 2D geometric algortihms. they are based on the boost library. 
* @details The idea is to define complex geometry based on shapes, usually
* multi-polygon using boost library. we propose only very simple combinaton
* that the region is composed of shapes without intersection.
* That is, the shapes are those contain each other or without overlap.
* This strict requirement suggests that complex shapes should be finished
* already in modeling using related binary operations before it is included.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

//boost library
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

#include "base_data_package.h"
#include "base_geometry.h"

#include <string>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

using namespace boost::geometry;

namespace SPH {

	typedef model::polygon<model::d2::point_xy<Real>> boost_poly;
	typedef model::multi_polygon<boost_poly> boost_multi_poly;

	/**
	 * @class MultiPolygon
	 * @brief used to define a closed region
	 */
	class MultiPolygon : public Shape
	{
	public:
		MultiPolygon() :Shape("MultiPolygon") {};
		boost_multi_poly& getBoostMultiPoly() { return multi_poly_; };
		bool checkContain(Vec2d pnt, bool BOUNDARY_INCLUDED = true);
		virtual Vec2d findClosestPoint(Vec2d input_pnt) override;
		virtual void findBounds(Vec2d &lower_bound, Vec2d &upper_bound) override;

		void addAMultiPolygon(MultiPolygon& multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly& boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(std::vector<Point>& points, ShapeBooleanOps op);
		void addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op);

	protected:
		boost_multi_poly multi_poly_;

		boost_multi_poly MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
						boost_multi_poly multi_poly_op, ShapeBooleanOps boolean_op);
	};

	/**
	 * @class ComplexShape
	 * @brief gives the final geomtrical definition of the SPHBody 
	 */
	class ComplexShape : public Shape
	{
	public:
		/** Default constructor. */
		ComplexShape() : Shape("ComplexShape"), multi_ploygen_() {};
		ComplexShape(string complex_shape_name) : Shape(complex_shape_name), multi_ploygen_() {};
		virtual ~ComplexShape() {};
		virtual Vec2d findClosestPoint(Vec2d input_pnt) override;
		virtual void findBounds(Vec2d& lower_bound, Vec2d& upper_bound) override;
		bool checkContain(Vec2d pnt, bool BOUNDARY_INCLUDED = true);

		void addAMultiPolygon(MultiPolygon& multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly& boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(std::vector<Point>& points, ShapeBooleanOps op);
		void addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op);

	protected:
		MultiPolygon multi_ploygen_;
	};
}

