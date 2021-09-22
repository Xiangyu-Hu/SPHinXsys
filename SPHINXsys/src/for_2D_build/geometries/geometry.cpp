/**
 * @file 	geometry.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "geometry.h"

using namespace boost::geometry;

namespace SPH
{
	//=================================================================================================//
	boost_multi_poly MultiPolygon::
		MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
								 boost_multi_poly multi_poly_op, ShapeBooleanOps boolean_op)
	{
		boost_multi_poly multi_poly_tmp_in = multi_poly_in;
		//out multi-poly need to be emtpy
		//otherwise the operation is not valid
		boost_multi_poly multi_poly_tmp_out;

		switch (boolean_op)
		{
		case ShapeBooleanOps::add:
		{
			boost::geometry::union_(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
			break;
		}

		case ShapeBooleanOps::sub:
		{
			boost::geometry::difference(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
			break;
		}
		case ShapeBooleanOps::sym_diff:
		{
			boost::geometry::sym_difference(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
			break;
		}
		case ShapeBooleanOps::intersect:
		{
			boost::geometry::intersection(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
			break;
		}
		default:
		{
			std::cout << "\n FAILURE: the type of boolean operation is undefined!" << std::endl;
			std::cout << "\n Please check the boost libraray reference." << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		}
		return multi_poly_tmp_out;
	}
	//=================================================================================================//
	void MultiPolygon::addAMultiPolygon(MultiPolygon &multi_polygon_op, ShapeBooleanOps op)
	{
		multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_polygon_op.getBoostMultiPoly(), op);
	}
	//=================================================================================================//
	void MultiPolygon::addABoostMultiPoly(boost_multi_poly &boost_multi_poly_op, ShapeBooleanOps op)
	{
		multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, boost_multi_poly_op, op);
	}
	//=================================================================================================//
	void MultiPolygon::addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op)
	{
		Vec2d buffer_center = center;
		Real buffer_radius = radius;
		int buffer_res = resolution;

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

		boost_multi_poly multi_poly_circle;
		buffer(circle_center_pnt, multi_poly_circle,
			   circle_dist_strategy, side_strategy,
			   join_strategy, end_strategy, circle_strategy);

		if (!is_valid(multi_poly_circle))
		{
			std::cout << "\n Error: the multi ploygen is not valid." << std::endl;
			std::cout << "\n The points must be in clockwise. Please check the boost libraray reference." << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}

		multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_circle, op);
	}
	//=================================================================================================//
	void MultiPolygon::addAPolygon(std::vector<Vecd> &points, ShapeBooleanOps op)
	{
		std::vector<model::d2::point_xy<Real>> pts;
		for (const Vecd &pnt : points)
		{
			pts.push_back(model::d2::point_xy<Real>(pnt[0], pnt[1]));
		}

		boost_poly poly;
		append(poly, pts);
		if (!is_valid(poly))
		{
			std::cout << "\n Try to reverse the points to clockwise." << std::endl;
			poly.clear();
			std::vector<model::d2::point_xy<Real>> pts_reverse(pts.rbegin(), pts.rend());
			append(poly, pts_reverse);
			if (!is_valid(poly))
			{
				std::cout << "\n Error: the multi ploygen is still not valid. Please check the boost libraray reference." << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				throw;
			}
		}

		boost_multi_poly multi_poly_polygen;
		convert(poly, multi_poly_polygen);

		multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_polygen, op);
	}
	//=================================================================================================//
	bool MultiPolygon::checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED /*= true*/)
	{
		if (BOUNDARY_INCLUDED)
		{
			return covered_by(model::d2::point_xy<Real>(pnt[0], pnt[1]), multi_poly_);
		}
		else
		{
			return within(model::d2::point_xy<Real>(pnt[0], pnt[1]), multi_poly_);
		}
	}
	//=================================================================================================//
	Vec2d MultiPolygon::findClosestPoint(const Vec2d &input_pnt)
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
		std::function<void(seg_type)> findclosestsegment = [&closest_seg, &closest_dist_2seg, &input_p](seg_type seg)
		{
			Real dist = boost::geometry::distance(input_p, seg);
			if (dist < closest_dist_2seg)
			{
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
		boost::geometry::for_each_segment(multi_poly_, findclosestsegment);

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
		if (c1 <= 0)
		{
			p_find = p_0;
		}
		else
		{
			Real c2 = dot(vec_v, vec_v);
			if (c2 <= c1)
			{
				p_find = p_1;
			}
			else
			{
				p_find = p_0 + vec_v * c1 / c2;
			}
		}

		return p_find;
	}
	//=================================================================================================//
	BoundingBox MultiPolygon::findBounds()
	{
		Vec2d lower_bound(0), upper_bound(0);
		typedef boost::geometry::model::box<model::d2::point_xy<Real>> box;
		lower_bound[0] = boost::geometry::return_envelope<box>(multi_poly_).min_corner().get<0>();
		lower_bound[1] = boost::geometry::return_envelope<box>(multi_poly_).min_corner().get<1>();
		upper_bound[0] = boost::geometry::return_envelope<box>(multi_poly_).max_corner().get<0>();
		upper_bound[1] = boost::geometry::return_envelope<box>(multi_poly_).max_corner().get<1>();
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
	bool ComplexShape::checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED)
	{
		return multi_ploygen_.checkContain(input_pnt, BOUNDARY_INCLUDED);
	}
	//=================================================================================================//
	Vec2d ComplexShape::findClosestPoint(const Vec2d &input_pnt)
	{
		return multi_ploygen_.findClosestPoint(input_pnt);
	}
	//=================================================================================================//
	bool ComplexShape::checkNotFar(const Vec2d &input_pnt, Real threshold)
	{
		return multi_ploygen_.checkContain(input_pnt) || checkNearSurface(input_pnt, threshold) ? true : false;
	}
	//=================================================================================================//
	bool ComplexShape::checkNearSurface(const Vec2d &input_pnt, Real threshold)
	{
		return getMaxAbsoluteElement(input_pnt - multi_ploygen_.findClosestPoint(input_pnt)) < threshold ? true : false;
	}
	//=================================================================================================//
	Real ComplexShape::findSignedDistance(const Vec2d &input_pnt)
	{
		Real distance_to_surface = (findClosestPoint(input_pnt) - input_pnt).norm();
		return checkContain(input_pnt) ? -distance_to_surface : distance_to_surface;
	}
	//=================================================================================================//
	Vec2d ComplexShape::findNormalDirection(const Vec2d &input_pnt)
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
	Vec2d ComplexShape::findNormalDirectionComplexShape(const Vec2d &input_pnt) //function to differentiate from LevelSetComplexShape::findNormalDirection
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
	BoundingBox ComplexShape::findBounds()
	{
		return multi_ploygen_.findBounds();
	}
	//=================================================================================================//
	void ComplexShape::addAMultiPolygon(MultiPolygon &multi_polygon, ShapeBooleanOps op)
	{
		multi_ploygen_.addAMultiPolygon(multi_polygon, op);
	}
	//=================================================================================================//
	void ComplexShape::addABoostMultiPoly(boost_multi_poly &boost_multi_poly, ShapeBooleanOps op)
	{
		multi_ploygen_.addABoostMultiPoly(boost_multi_poly, op);
	}
	//=================================================================================================//
	void ComplexShape::addAPolygon(std::vector<Vecd> &points, ShapeBooleanOps op)
	{
		multi_ploygen_.addAPolygon(points, op);
	}
	//=================================================================================================//
	void ComplexShape::
		addAPolygonFromFile(std::string file_path_name, ShapeBooleanOps op, Vec2d translation, Real scale_factor)
	{
		std::fstream dataFile(file_path_name);
		Vecd temp_point;
		std::vector<Vecd> coordinates;
		double temp1 = 0.0, temp2 = 0.0;
		if (dataFile.fail())
		{
			std::cout << "File can not open.\n"
					  << std::endl;
			;
		}

		while (!dataFile.fail() && !dataFile.eof())
		{
			dataFile >> temp1 >> temp2;
			temp_point[0] = temp1 * scale_factor + translation[0];
			temp_point[1] = temp2 * scale_factor + translation[1];
			coordinates.push_back(temp_point);
		}
		dataFile.close();

		multi_ploygen_.addAPolygon(coordinates, op);
	}
	//=================================================================================================//
	void ComplexShape::addACircle(Vec2d center, Real radius, int resolution, ShapeBooleanOps op)
	{
		multi_ploygen_.addACircle(center, radius, resolution, op);
	}
	//=================================================================================================//
}