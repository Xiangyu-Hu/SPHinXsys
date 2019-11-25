/**
 * @file 	geometry.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "geometry.h"

using namespace boost::geometry;

namespace SPH {

	Geometry::Geometry(Vec2d center_, Real radius_, int resolution_) : Shape("MLPolygon")
	{
		Vec2d buffer_center = center_;
		Real buffer_radius = radius_;
		int buffer_res = resolution_;

		// Declare the point_circle strategy
		strategy::buffer::join_round join_strategy;
		strategy::buffer::end_round end_strategy;
		strategy::buffer::side_straight side_strategy;
		strategy::buffer::point_circle circle_strategy(buffer_res);
		strategy::buffer::distance_symmetric<double> circle_dist_strategy(buffer_radius);

		// Create the buffer of a multi point
		model::d2::point_xy<Real> circle_center_pnt;

		boost::geometry::set<0>(circle_center_pnt, buffer_center[0]);
		boost::geometry::set<1>(circle_center_pnt, buffer_center[1]);

		buffer(circle_center_pnt, multi_poly,
			circle_dist_strategy, side_strategy,
			join_strategy, end_strategy, circle_strategy);

		if (!is_valid(multi_poly)) {
			std::cout << "\n Error: the multi ploygen is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//===========================================================//
	Geometry::Geometry(boost_multi_poly multi_poly_)
		: Shape("MLPolygon")
	{
		multi_poly = multi_poly_;
	}
	//===========================================================//
	Geometry::Geometry(std::vector<Point>& points)
		: Shape("MLPolygon")
	{
		std::vector<model::d2::point_xy<Real>> pts;
		for (const Point& pnt : points)
		{
			pts.push_back(model::d2::point_xy<Real>(pnt[0], pnt[1]));
		}
		
		boost_poly poly;
		append(poly, pts);
		if (!is_valid(poly)) {
			std::cout << "\n Error: the ploygen is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}

		convert(poly, multi_poly);
	}
	//===========================================================//
	boost_multi_poly  Geometry::get_multi_poly()
	{
		return multi_poly;
	}
	//===========================================================//
	bool Geometry::contain(Vec2d pnt, bool BOUNDARY_INCLUDED /*= true*/)
	{
		if (BOUNDARY_INCLUDED)
		{
			return covered_by(model::d2::point_xy<Real>(pnt[0], pnt[1]), multi_poly);
		}
		else
		{
			return within(model::d2::point_xy<Real>(pnt[0], pnt[1]), multi_poly);
		}
	}
	//===========================================================//
	Vec2d Geometry::closestpointonface(Vec2d input_pnt)
	{
		typedef model::d2::point_xy<Real> pnt_type;
		typedef model::referring_segment<model::d2::point_xy<Real>> seg_type;
		/*
		typedef model::segment<model::d2::point_xy<Real>> seg_type;
		From the documentation on segment and referring_segment, the only difference between the two is that
		referring_segment holds a reference to the points.
		This is what is needed in a for each that modifies the segment since the points modified should be
		reflected in the linestring. In a for each that does not modify the points, it should still take a
		reference (most likely a const reference) since it reduces the amount of copying.
		*/
		pnt_type input_p(input_pnt[0], input_pnt[1]);
		model::segment<model::d2::point_xy<Real>> closest_seg;
		Real closest_dist_2seg = boost::numeric::bounds<Real>::highest();
		std::function<void(seg_type)> findclosestsegment = [&closest_seg, &closest_dist_2seg, &input_p](seg_type seg) {
			Real dist = boost::geometry::distance(input_p, seg);
			if (dist < closest_dist_2seg) {
				closest_dist_2seg = dist;
				//closest_seg.append(seg);
				Real x0 = boost::geometry::get<0, 0>(seg);
				Real y0 = boost::geometry::get<0, 1>(seg);
				Real x1 = boost::geometry::get<1, 0>(seg);
				Real y1 = boost::geometry::get<1, 1>(seg);
				boost::geometry::set<0, 0>(closest_seg, x0);
				boost::geometry::set<0, 1>(closest_seg, y0);
				boost::geometry::set<1, 0>(closest_seg, x1);
				boost::geometry::set<1, 1>(closest_seg, y1);
			}
		};
		boost::geometry::for_each_segment(multi_poly, findclosestsegment);

		Vec2d p_find(0, 0);

		Real x0 = boost::geometry::get<0, 0>(closest_seg);
		Real y0 = boost::geometry::get<0, 1>(closest_seg);
		Real x1 = boost::geometry::get<1, 0>(closest_seg);
		Real y1 = boost::geometry::get<1, 1>(closest_seg);
		Vec2d p_0(x0, y0);
		Vec2d p_1(x1, y1);
		Vec2d vec_v = p_1 - p_0;
		Vec2d vec_w = input_pnt - p_0;

		Real c1 = dot(vec_v, vec_w);
		if (c1 <= 0) {
			p_find = p_0;
		}
		else {
			Real c2 = dot(vec_v, vec_v);
			if (c2 <= c1) {
				p_find = p_1;
			}
			else {
				p_find = p_0 + vec_v * c1 / c2;
			}
		}
		
		return p_find;
	}
	//===========================================================//
	void Geometry::shapebound(Vec2d &lower_bound, Vec2d &upper_bound)
	{
		typedef boost::geometry::model::box<model::d2::point_xy<Real>> box;
		lower_bound[0] = boost::geometry::return_envelope<box>(multi_poly).min_corner().get<0>();
		lower_bound[1] = boost::geometry::return_envelope<box>(multi_poly).min_corner().get<1>();
		upper_bound[0] = boost::geometry::return_envelope<box>(multi_poly).max_corner().get<0>();
		upper_bound[1] = boost::geometry::return_envelope<box>(multi_poly).max_corner().get<1>();
	}
	//===========================================================//
	Region::Region(string region_name)
	{
		region_name_ = region_name;
	}
	//===========================================================//
	bool Region::contain(Vec2d pnt, bool BOUNDARY_INCLUDED /*= true*/)
	{
		return Geometry::contain(pnt);
	}
	//===========================================================//
	void Region::closestpointonface(Vec2d input_pnt, Vec2d& closest_pnt, Real& phi)
	{
		closest_pnt = Geometry::closestpointonface(input_pnt);
		Real phii = (closest_pnt - input_pnt).norm();
		phi = contain(input_pnt) ? phii : -phii;
	}
	//===========================================================//
	void Region::regionbound(Vec2d &lower_bound, Vec2d &upper_bound)
	{
		Geometry::shapebound(lower_bound, upper_bound);
	}
	//===========================================================//
	void Region::add_geometry(Geometry *geometry, RegionBooleanOps op)
	{
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}
	//===========================================================//
	void Region::add_polygon(std::vector<Point>& points, RegionBooleanOps op)
	{
		Geometry * geometry = new Geometry(points);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}
	//===========================================================//
	void Region::add_circle(Vec2d center_, Real radius_, int resolution_, RegionBooleanOps op)
	{
		Geometry * geometry = new Geometry(center_, radius_, resolution_);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}
	//===========================================================//
	void Region::add_boost_multi_polygon(boost_multi_poly multi_poly_, RegionBooleanOps op)
	{
		Geometry * geometry = new Geometry(multi_poly_);
		geometries.push_back(geometry);
		geometryops.push_back(op);
	}
	//===========================================================//
	void Region::done_modeling()
	{
		boost_multi_poly multi_poly_tmp_in = geometries[0]->get_multi_poly();

		for (size_t i = 1; i < geometries.size(); ++i) {
			//out multi-poly need to be emtpy
			//otherwise the operation is not valid
			boost_multi_poly multi_poly_tmp_out;

			RegionBooleanOps boolean_op = geometryops[i];
			switch (boolean_op)
			{
			case RegionBooleanOps::add: {
				boost::geometry::union_(multi_poly_tmp_in,
					geometries[i]->get_multi_poly(), multi_poly_tmp_out);
				break;
			}

			case RegionBooleanOps::sub: {
				boost::geometry::difference(multi_poly_tmp_in,
					geometries[i]->get_multi_poly(), multi_poly_tmp_out);
				break;
			}
			case RegionBooleanOps::sym_diff: {
				boost::geometry::sym_difference(multi_poly_tmp_in,
					geometries[i]->get_multi_poly(), multi_poly_tmp_out);
				break;
			}
			case RegionBooleanOps::intersect: {
				boost::geometry::intersection(multi_poly_tmp_in,
					geometries[i]->get_multi_poly(), multi_poly_tmp_out);
				break;
			}
			default:
			{
				std::cout << "\n FAILURE: the type of boolean operation is undefined!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
				break;
			}
			}
			multi_poly_tmp_in = multi_poly_tmp_out;
		}

		if (!is_valid(multi_poly_tmp_in)) {
			std::cout << "\n FAILURE: the boolean operation of multi ploygen is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		multi_poly = multi_poly_tmp_in;
	}
}