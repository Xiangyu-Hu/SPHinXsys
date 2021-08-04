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
*/

#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

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
	class MultiPolygon : public Shape
	{
	public:
		MultiPolygon() : Shape("MultiPolygon"){};
		boost_multi_poly &getBoostMultiPoly() { return multi_poly_; };
		bool checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED = true);
		Vec2d findClosestPoint(const Vec2d &input_pnt);
		virtual BoundingBox findBounds() override;

		void addAMultiPolygon(MultiPolygon &multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly &boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(std::vector<Vecd> &points, ShapeBooleanOps op);
		void addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op);

	protected:
		boost_multi_poly multi_poly_;
		boost_multi_poly MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
												  boost_multi_poly multi_poly_op,
												  ShapeBooleanOps boolean_op);
	};

	/**
	 * @class ComplexShape
	 * @brief gives the final geomtrical definition of the SPHBody 
	 */
	class ComplexShape : public Shape
	{
		Vec2d findClosestPoint(const Vec2d &input_pnt);
	public:
		/** Default constructor. */
		ComplexShape() : Shape("ComplexShape"), multi_ploygen_(){};
		ComplexShape(std::string complex_shape_name) : Shape(complex_shape_name), multi_ploygen_(){};
		virtual ~ComplexShape(){};
		virtual BoundingBox findBounds() override;
		void addAMultiPolygon(MultiPolygon &multi_polygon, ShapeBooleanOps op);
		void addABoostMultiPoly(boost_multi_poly &boost_multi_poly, ShapeBooleanOps op);
		void addAPolygon(std::vector<Vecd> &points, ShapeBooleanOps op);
		void addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op);
		void addAPolygonFromFile(std::string file_path_name,
								 ShapeBooleanOps op,
								 Vec2d translation = Vecd(0),
								 Real scale_factor = 1.0);

		virtual bool checkContain(const Vec2d &input_pnt, bool BOUNDARY_INCLUDED = true);
		virtual bool checkNotFar(const Vec2d &input_pnt, Real threshold);
		virtual bool checkNearSurface(const Vec2d &input_pnt, Real threshold);
		/** Signed distance is negative for point within the complex shape. */
		virtual Real findSignedDistance(const Vec2d &input_pnt);
		/** Normal direction point toward outside of the complex shape. */
		virtual Vec2d findNormalDirection(const Vec2d &input_pnt);

	protected:
		MultiPolygon multi_ploygen_;
	};
}

#endif //GEOMETRY_2D_H