/**
* @file geoemtry.h
* @brief Here, we define the 2D geomtric algortihms. they are based on the boost library. 
* @details The idea is to define complex geomerty based on shapes, usually
* multi-polygeon using boost library. we propose only very simple combinaton
* that the region is composed of shapes without intersection.
* That is, the shapes are those contain each other or without overlap.
* This strict requirement suggests that complex shapes should be finished
* already in modeling using related binary oparations before it is included.
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

typedef model::polygon<model::d2::point_xy<Real>> boost_poly;
typedef model::multi_polygon<boost_poly> boost_multi_poly;

namespace SPH {


	/**
	 * @class RegionBooleanOps
	 * @brief Boolian operation for generate complex 2D regions
	 * @details Note that add and sub boolean operation is only defined for particle generator right now
	 * in generating final conbined geometry for level set initialization
	 * add denotes boost::geomerty::union_
	 * sub denotes  boost::geometry::difference
	 * sym_diff denotes boost::geometry::sym_difference
	 * intersect denotes boost::geometry::intersection
	 * also, the final combined geoemtry can be used for particle generator, I will think about this later
	 */
	enum class RegionBooleanOps {add, sub, sym_diff, intersect};

	/**
	 * @class Geometry
	 * @brief One geometric which will be one component of a region
	 */
	class Geometry : public Shape
	{
	protected:
		boost_multi_poly multi_poly;

	public:

		Geometry() :Shape("MLPolygon") {};
		/** a circle */
		Geometry(Vec2d center_, Real radius_, int resolution_);
		/** a multi-polygon */
		Geometry(boost_multi_poly multi_poly_);
		/** a polygon */
		Geometry(std::vector<Point>& points);

		boost_multi_poly get_multi_poly();
		virtual bool contain(Vec2d pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual Vec2d closestpointonface(Vec2d input_pnt) override;
		virtual void shapebound(Vec2d &lower_bound, Vec2d &upper_bound) override;
	};

	/**
	 * @class Region
	 * @brief It will give the final geoemtrical definition of the SPHBody 
	 */
	class Region : public Geometry
	{
	protected:
		/** name of the region */
		std::string region_name_;
		/** geometry container */
		std::vector<Geometry*> geometries;
		/** geometry operation container */
		std::vector<RegionBooleanOps> geometryops;

	public:
		Region(string region_name);
		virtual ~Region() {};
		virtual bool contain(Vec2d pnt, bool BOUNDARY_INCLUDED = true);
		virtual void closestpointonface(Vec2d input_pnt, Vec2d& closest_pnt, Real& phi);
		virtual void regionbound(Vec2d &lower_bound, Vec2d &upper_bound);

		void add_geometry(Geometry *geometry, RegionBooleanOps op);
		/** add a polygen */
		void add_polygon(std::vector<Point>& points, RegionBooleanOps op);
		/** add a circle */
		void add_circle(Vec2d center_, Real radius_, int resolution_, RegionBooleanOps op);
		/** add a multi-polygon */
		void add_boost_multi_polygon(boost_multi_poly multi_poly_, RegionBooleanOps op);
		/** finish the final geometric modeling */
		void done_modeling();
	};
}

