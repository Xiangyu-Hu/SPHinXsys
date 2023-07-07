#include "multi_polygon_shape.h"

using namespace boost::geometry;

namespace SPH
{
//=================================================================================================//
MultiPolygon::MultiPolygon(const std::vector<Vecd> &points)
    : MultiPolygon()
{
    addAPolygon(points, ShapeBooleanOps::add);
}
//=================================================================================================//
MultiPolygon::MultiPolygon(const Vecd &center, Real radius, int resolution)
    : MultiPolygon()
{
    addACircle(center, radius, resolution, ShapeBooleanOps::add);
}
//=================================================================================================//
boost_multi_poly MultiPolygon::
    MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
                             boost_multi_poly multi_poly_op, ShapeBooleanOps boolean_op)
{
    boost_multi_poly multi_poly_tmp_in = multi_poly_in;
    /**
     * Out multi-poly need to be emtpy
     * otherwise the operation is not valid.
     */
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
        std::cout << "\n Please check the boost library reference." << std::endl;
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
void MultiPolygon::addABox(Transform transform, const Vecd &halfsize, ShapeBooleanOps op)
{
    Vecd point0 = transform.shiftFrameStationToBase(-halfsize);
    Vecd point1 = transform.shiftFrameStationToBase(Vecd(-halfsize[0], halfsize[1]));
    Vecd point2 = transform.shiftFrameStationToBase(halfsize);
    Vecd point3 = transform.shiftFrameStationToBase(Vecd(halfsize[0], -halfsize[1]));

    std::vector<Vecd> points = {point0, point1, point2, point3, point0};
    addAPolygon(points, op);
}
//=================================================================================================//
void MultiPolygon::addACircle(const Vecd &center, Real radius, int resolution, ShapeBooleanOps op)
{
    Vecd buffer_center = center;
    Real buffer_radius = radius;
    int buffer_res = resolution;

    // Declare the point_circle strategy
    strategy::buffer::join_round join_strategy;
    strategy::buffer::end_round end_strategy;
    strategy::buffer::side_straight side_strategy;
    strategy::buffer::point_circle circle_strategy(buffer_res);
    strategy::buffer::distance_symmetric<Real> circle_dist_strategy(buffer_radius);

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
        std::cout << "\n Error: the multi polygon is not valid." << std::endl;
        std::cout << "\n The points must be in clockwise. Please check the boost library reference." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }

    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_circle, op);
}
//=================================================================================================//
void MultiPolygon::addAPolygon(const std::vector<Vecd> &points, ShapeBooleanOps op)
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
            std::cout << "\n Error: the multi polygon is still not valid. Please check the boost library reference." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            throw;
        }
    }

    boost_multi_poly multi_poly_polygon;
    convert(poly, multi_poly_polygon);

    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_polygon, op);
}
//=================================================================================================//
void MultiPolygon::
    addAPolygonFromFile(std::string file_path_name, ShapeBooleanOps op, Vecd translation, Real scale_factor)
{
    std::fstream dataFile(file_path_name);
    Vecd temp_point;
    std::vector<Vecd> coordinates;
    Real temp1 = 0.0, temp2 = 0.0;
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

    addAPolygon(coordinates, op);
}
//=================================================================================================//
bool MultiPolygon::checkContain(const Vec2d &probe_point, bool BOUNDARY_INCLUDED /*= true*/)
{
    if (BOUNDARY_INCLUDED)
    {
        return covered_by(model::d2::point_xy<Real>(probe_point[0], probe_point[1]), multi_poly_);
    }
    else
    {
        return within(model::d2::point_xy<Real>(probe_point[0], probe_point[1]), multi_poly_);
    }
}
//=================================================================================================//
Vecd MultiPolygon::findClosestPoint(const Vecd &probe_point)
{
    typedef model::d2::point_xy<Real> pnt_type;
    typedef model::referring_segment<model::d2::point_xy<Real>> seg_type;
    /**
     * typedef model::segment<model::d2::point_xy<Real>> seg_type;
     * From the documentation on segment and referring_segment, the only difference between the two is that
     * referring_segment holds a reference to the points.
     * This is what is needed in a for each that modifies the segment since the points modified should be
     * reflected in the line string. In a for each that does not modify the points, it should still take a
     * reference (most likely a const reference) since it reduces the amount of copying.
     */
    pnt_type input_p(probe_point[0], probe_point[1]);
    model::segment<model::d2::point_xy<Real>> closest_seg;
    Real closest_dist_2seg = boost::numeric::bounds<Real>::highest();
    std::function<void(seg_type)> findclosestsegment = [&closest_seg, &closest_dist_2seg, &input_p](seg_type seg)
    {
        Real dist = boost::geometry::distance(input_p, seg);
        if (dist < closest_dist_2seg)
        {
            closest_dist_2seg = dist;
            // closest_seg.append(seg);
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

    Vecd p_find = Vecd::Zero();

    Real x0 = boost::geometry::get<0, 0>(closest_seg);
    Real y0 = boost::geometry::get<0, 1>(closest_seg);
    Real x1 = boost::geometry::get<1, 0>(closest_seg);
    Real y1 = boost::geometry::get<1, 1>(closest_seg);
    Vecd p_0(x0, y0);
    Vecd p_1(x1, y1);
    Vecd vec_v = p_1 - p_0;
    Vecd vec_w = probe_point - p_0;

    Real c1 = vec_v.dot(vec_w);
    if (c1 <= 0)
    {
        p_find = p_0;
    }
    else
    {
        Real c2 = vec_v.dot(vec_v);
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
    Vecd lower_bound = Vecd::Zero();
    Vecd upper_bound = Vecd::Zero();
    typedef boost::geometry::model::box<model::d2::point_xy<Real>> box;
    lower_bound[0] = boost::geometry::return_envelope<box>(multi_poly_).min_corner().get<0>();
    lower_bound[1] = boost::geometry::return_envelope<box>(multi_poly_).min_corner().get<1>();
    upper_bound[0] = boost::geometry::return_envelope<box>(multi_poly_).max_corner().get<0>();
    upper_bound[1] = boost::geometry::return_envelope<box>(multi_poly_).max_corner().get<1>();
    return BoundingBox(lower_bound, upper_bound);
}
//=================================================================================================//
bool MultiPolygonShape::isValid()
{
    return multi_polygon_.getBoostMultiPoly().size() == 0 ? false : true;
}
//=================================================================================================//
bool MultiPolygonShape::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return multi_polygon_.checkContain(probe_point, BOUNDARY_INCLUDED);
}
//=================================================================================================//
Vecd MultiPolygonShape::findClosestPoint(const Vecd &probe_point)
{
    return multi_polygon_.findClosestPoint(probe_point);
}
//=================================================================================================//
BoundingBox MultiPolygonShape::findBounds()
{
    return multi_polygon_.findBounds();
}
//=================================================================================================//
} // namespace SPH